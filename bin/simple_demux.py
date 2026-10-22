import pysam
import argparse
import os
import json
import statistics
import time
import multiprocessing as mp
from itertools import islice
from functools import partial
from collections import Counter, defaultdict

UMI_PATTERN = "ATCAACGCAGAGTNNNNTNNNNTNNNNTCTTGGG"

# DNA Reverse Complement Map
REV_COMP_MAP = str.maketrans("ATGC", "TACG")

def reverse_complement(seq):
    return seq.translate(REV_COMP_MAP)[::-1]

def hamming_distance(s1, s2):
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))

def parse_fasta(fasta_path):
    barcodes = {}
    with pysam.FastxFile(fasta_path) as fh:
        for entry in fh:
            barcodes[entry.name] = entry.sequence.upper()
    return barcodes

def barcode_key(name):
    """'RACE_TSO_5p--bc1001' (lima-style pair) or 'bc1001' -> 'bc1001'."""
    return name.strip().split('--')[-1]

def parse_biosample(csv_path):
    mapping = {}
    with open(csv_path, 'r') as f:
        for line in f:
            if line.startswith("SampleName") or not line.strip():
                continue
            parts = [p.strip() for p in line.strip().split(',')]
            if len(parts) >= 2:
                mapping[barcode_key(parts[1])] = parts[0]
    return mapping

def score_read_scan(target_seq, bc_refs, scan_end=1):
    best_name, best_dist, best_idx = None, 99, -1
    for name, ref in bc_refs.items():
        ref_len = len(ref)
        for i in range(0, scan_end):
            read_slice = target_seq[i : i + ref_len]
            if len(read_slice) < ref_len: break
            dist = hamming_distance(read_slice, ref)
            if dist <= 1: return name, dist, i
            if dist < best_dist:
                best_dist, best_name, best_idx = dist, name, i
            elif dist == best_dist and name != best_name:
                best_name = "AMBIGUOUS"
    return best_name, best_dist, best_idx

def is_umi_valid(seq, umi):
    for i, base in enumerate(umi):
        if base != "N" and seq[i] != base:
            return False
    return True


def process_read_chunk(chunk, tso_seq, sample_bcs, mismatches, max_bc_dist, no_trim, tso_window=1, bc_window=1, bc_trim_len=0):
    results = []
    tso_len = len(list(tso_seq.values())[0])
    for zmw, seq, qual in chunk:
        fwd_name, fwd_dist, fwd_pos = score_read_scan(seq, tso_seq, tso_window)
        rev_seq = reverse_complement(seq)
        rev_name, rev_dist, rev_pos = score_read_scan(rev_seq, tso_seq, tso_window)

        orientation, working_seq = "UNKNOWN", None
        final_seq, final_qual = seq, qual
        tso_dist, tso_pos = 99, -1

        if fwd_dist <= rev_dist and fwd_dist <= mismatches:
            orientation, tso_dist, tso_pos, working_seq = "FORWARD", fwd_dist, fwd_pos, seq
        elif rev_dist <= mismatches:
            orientation, tso_dist, tso_pos, working_seq = "REVERSE", rev_dist, rev_pos, rev_seq
            final_seq, final_qual = rev_seq, qual[::-1]

        bc_name, bc_dist, bc_pos = "N/A", 99, -1
        if working_seq:
            bc_search_seq = reverse_complement(working_seq)
            bc_name, bc_dist, bc_pos = score_read_scan(bc_search_seq, sample_bcs, bc_window)

        umi = working_seq[tso_pos + tso_len:tso_pos + tso_len + len(UMI_PATTERN)] if tso_pos >= 0 else ""

        if len(final_seq) < 500:
            decision = "REJECTED - READ TOO SHORT"
        elif len(final_seq) > 2500:
            decision = "REJECTED - READ TOO LONG"
        elif orientation == "UNKNOWN":
            decision = "REJECTED - TSO NOT FOUND"
        elif bc_name == "AMBIGUOUS":
            decision = "REJECTED - AMBIGUOUS BC"
        elif bc_dist > max_bc_dist:
            decision = "REJECTED - BC DISTANCE TOO HIGH"
        elif not is_umi_valid(umi, UMI_PATTERN):
            decision = "REJECTED - UMI MISMATCH"
        else:
            decision = "ACCEPTED"
            if not no_trim:
                bc_len = len(sample_bcs[bc_name])
                # bc_trim_len > 0: cut only the outer N bases of the matched entry
                # (e.g. the 10-bp sample barcode) and keep the constant-region primer
                if 0 < bc_trim_len < bc_len:
                    bc_len = bc_trim_len
                trim_start = tso_pos + tso_len
                trim_end = len(final_seq) - bc_pos - bc_len
                if trim_start < trim_end:
                    final_seq = final_seq[trim_start:trim_end]
                    final_qual = final_qual[trim_start:trim_end]
            
        log_line = f"{zmw},{umi},{orientation},{len(final_seq)},{tso_dist},{tso_pos},{bc_name},{bc_dist},{bc_pos},{decision}\n"
        results.append((zmw, final_seq, final_qual, bc_name, decision, log_line, tso_pos, bc_pos, bc_dist))
    return results


def make_histogram(lengths, bin_size=50):
    if not lengths:
        return [], []
    max_l = max(lengths)
    bins = list(range(0, max_l + bin_size, bin_size))
    counts = [0] * len(bins)
    for l in lengths:
        idx = l // bin_size
        if idx < len(counts):
            counts[idx] += 1
    return bins, counts


def write_read_length_report(length_data, prefix, suffix):
    """Write read length distribution CSV and HTML from collected lengths."""
    all_assigned = []
    all_unassigned = []

    for name, lengths in length_data.items():
        if name == 'unassigned':
            all_unassigned.extend(lengths)
        else:
            all_assigned.extend(lengths)

    # CSV
    csv_path = f"{prefix}read_length_summary{suffix}.csv"
    with open(csv_path, 'w') as f:
        f.write("category,sample,count,mean_length,median_length,min_length,max_length\n")
        for name, lengths in sorted(length_data.items()):
            if not lengths:
                continue
            s = sorted(lengths)
            n = len(s)
            cat = "unassigned" if name == "unassigned" else "assigned"
            f.write(f"{cat},{name},{n},{sum(s)/n:.1f},{s[n//2]},{min(s)},{max(s)}\n")
        if all_assigned:
            a = sorted(all_assigned)
            n = len(a)
            f.write(f"assigned,ALL_ASSIGNED,{n},{sum(a)/n:.1f},{a[n//2]},{min(a)},{max(a)}\n")
        if all_unassigned:
            u = sorted(all_unassigned)
            n = len(u)
            f.write(f"unassigned,ALL_UNASSIGNED,{n},{sum(u)/n:.1f},{u[n//2]},{min(u)},{max(u)}\n")

    # HTML
    traces = []
    colors = [
        '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
        '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac',
        '#86bcb6', '#d37295', '#8cd17d', '#b6992d', '#499894'
    ]

    bins_a, counts_a = make_histogram(all_assigned)
    bins_u, counts_u = make_histogram(all_unassigned)
    traces.append({
        'x': bins_a, 'y': counts_a,
        'name': f'Assigned (n={len(all_assigned):,})',
        'type': 'bar', 'opacity': 0.7,
        'marker': {'color': '#4e79a7'}
    })
    traces.append({
        'x': bins_u, 'y': counts_u,
        'name': f'Unassigned (n={len(all_unassigned):,})',
        'type': 'bar', 'opacity': 0.7,
        'marker': {'color': '#e15759'}
    })

    color_idx = 0
    for name, lengths in sorted(length_data.items()):
        if name == 'unassigned' or not lengths:
            continue
        bins, counts = make_histogram(lengths)
        traces.append({
            'x': bins, 'y': counts,
            'name': f'{name} (n={len(lengths):,})',
            'type': 'bar', 'opacity': 0.6,
            'visible': 'legendonly',
            'marker': {'color': colors[color_idx % len(colors)]}
        })
        color_idx += 1

    title_str = suffix.replace('_', ' ').strip() if suffix else 'All'
    html_path = f"{prefix}read_length_dist{suffix}.html"
    with open(html_path, 'w') as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
<title>Read Length Distributions - {title_str}</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
    body {{ font-family: Arial, sans-serif; margin: 20px; background: #fafafa; }}
    h1 {{ color: #333; }}
    table {{ border-collapse: collapse; margin: 20px 0; }}
    th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: right; }}
    th {{ background: #4e79a7; color: white; }}
    tr:nth-child(even) {{ background: #f2f2f2; }}
    .assigned {{ color: #4e79a7; font-weight: bold; }}
    .unassigned {{ color: #e15759; font-weight: bold; }}
</style>
</head>
<body>
<h1>Read Length Distributions: {title_str}</h1>

<div id="plot" style="width:100%;max-width:1200px;height:500px;"></div>

<script>
var traces = {json.dumps(traces)};
var layout = {{
    barmode: 'overlay',
    title: 'Read Length Distribution',
    xaxis: {{ title: 'Read Length (bp)', dtick: 250 }},
    yaxis: {{ title: 'Count' }},
    legend: {{ x: 1, xanchor: 'right', y: 1 }},
    template: 'plotly_white'
}};
Plotly.newPlot('plot', traces, layout, {{responsive: true}});
</script>

<h2>Summary Statistics</h2>
<table>
<tr><th>Category</th><th>Sample</th><th>Reads</th><th>Mean Length</th><th>Median Length</th><th>Min</th><th>Max</th></tr>
""")
        for name, lengths in sorted(length_data.items()):
            if not lengths:
                continue
            s = sorted(lengths)
            n = len(s)
            cat = "unassigned" if name == "unassigned" else "assigned"
            f.write(f'<tr><td class="{cat}">{cat}</td><td>{name}</td>'
                    f'<td>{n:,}</td><td>{sum(s)/n:,.1f}</td><td>{s[n//2]:,}</td>'
                    f'<td>{min(s):,}</td><td>{max(s):,}</td></tr>\n')
        if all_assigned:
            a = sorted(all_assigned)
            n = len(a)
            f.write(f'<tr style="font-weight:bold"><td class="assigned">assigned</td><td>TOTAL</td>'
                    f'<td>{n:,}</td><td>{sum(a)/n:,.1f}</td><td>{a[n//2]:,}</td>'
                    f'<td>{min(a):,}</td><td>{max(a):,}</td></tr>\n')
        if all_unassigned:
            u = sorted(all_unassigned)
            n = len(u)
            f.write(f'<tr style="font-weight:bold"><td class="unassigned">unassigned</td><td>TOTAL</td>'
                    f'<td>{n:,}</td><td>{sum(u)/n:,.1f}</td><td>{u[n//2]:,}</td>'
                    f'<td>{min(u):,}</td><td>{max(u):,}</td></tr>\n')
        f.write("</table>\n</body>\n</html>")

    print(f"Read length CSV: {csv_path}")
    print(f"Read length HTML: {html_path}")


def run_demux(args):
    start_time = time.time()
    next_progress_time = start_time + 15

    bc_all = parse_fasta(args.barcodes)
    sample_map = parse_biosample(args.biosample)
    tso_name = "RACE_TSO_5p"
    tso_seq = {tso_name: bc_all[tso_name]}
    print(tso_seq)
    sample_bcs = {k: v for k, v in bc_all.items() if tso_name not in k}
    missing = sorted(set(sample_map) - set(sample_bcs))
    if missing:
        print(f"WARNING: barcodes in --biosample not found in --barcodes FASTA: {missing}")
        print(f"         FASTA barcode names: {sorted(sample_bcs)}")
    print(f"Biosample map: {sample_map}")

    prefix = args.prefix
    suffix = args.suffix
    stats = Counter()
    tso_positions = []
    bc_positions = []
    bc_distances = []
    fq_writers = {}

    # Track read lengths per sample for inline distribution report
    read_lengths = defaultdict(list)
    
    with pysam.AlignmentFile(args.input, "rb", check_sq=False) as reader, \
         open(f"{prefix}audit_log.csv", "w") as log_file, \
         open(f"{prefix}unassigned{suffix}.fastq", "w") as unassigned_fq:
        
        log_file.write("ZMW,umi,Orientation,Read_Length,TSO_Dist,TSO_Pos,BC_Match,BC_Dist,BC_Pos,Decision\n")
        
        def read_chunker():
            fetch_iter = reader.fetch(until_eof=True)
            while True:
                chunk = list(islice(fetch_iter, args.chunk_size))
                if not chunk:
                    break
                yield [(r.query_name, r.query_sequence, "".join([chr(q + 33) for q in r.query_qualities]))
                       for r in chunk if r.query_sequence is not None and r.query_qualities is not None]

        worker_func = partial(
            process_read_chunk, 
            tso_seq=tso_seq, 
            sample_bcs=sample_bcs, 
            mismatches=args.mismatches, 
            max_bc_dist=args.max_bc_dist,
            no_trim=args.no_trim,
            tso_window=args.tso_window,
            bc_window=args.bc_window,
            bc_trim_len=args.bc_trim_len
        )
        
        with mp.Pool(processes=args.cpus) as pool:
            for result_chunk in pool.imap(worker_func, read_chunker()):
                for zmw, final_seq, final_qual, bc_name, decision, log_line, tso_pos, bc_pos, bc_dist in result_chunk:
                    stats['total'] += 1
                    log_file.write(log_line)

                    if decision == "ACCEPTED":
                        s_name = sample_map.get(bc_name, bc_name)
                        out_path = f"{prefix}{s_name}{suffix}.fastq"
                        if out_path not in fq_writers: fq_writers[out_path] = open(out_path, "w")
                        fq_writers[out_path].write(f"@{zmw}\n{final_seq}\n+\n{final_qual}\n")
                        stats[s_name] += 1
                        tso_positions.append(tso_pos)
                        bc_positions.append(bc_pos)
                        bc_distances.append(bc_dist)
                        read_lengths[s_name].append(len(final_seq))
                    else:
                        unassigned_fq.write(f"@{zmw}\n{final_seq}\n+\n{final_qual}\n")
                        stats['unassigned'] += 1
                        read_lengths['unassigned'].append(len(final_seq))
                        if decision == "REJECTED - UMI MISMATCH":
                            stats['umi_mismatch'] += 1

                    now = time.time()
                    if now >= next_progress_time:
                        print(f"[{time.strftime('%H:%M:%S')}] Processed {stats['total']} records...")
                        next_progress_time = now + 15

    for f in fq_writers.values(): f.close()
    
    print("\n" + "="*30)
    print("POSITIONING SUMMARY")
    print("="*30)
    if tso_positions:
        print(f"TSO Median Index: {statistics.median(tso_positions):.1f} (SD: {statistics.pstdev(tso_positions):.2f})")
        print(f"BC Median Index:  {statistics.median(bc_positions):.1f} (SD: {statistics.pstdev(bc_positions):.2f})")
        
        print("\nPosition Frequencies:")
        print(f"{'Position':<10} | {'TSO Count':<12} | {'BC Count':<12}")
        print("-" * 40)
        tso_counts = Counter(tso_positions)
        bc_counts = Counter(bc_positions)
        all_pos = sorted(set(tso_counts.keys()).union(bc_counts.keys()))
        for pos in all_pos:
            print(f"{pos:<10} | {tso_counts.get(pos, 0):<12} | {bc_counts.get(pos, 0):<12}")

    if bc_distances:
        print("\n" + "="*30)
        print("BARCODE DISTANCE SUMMARY (ACCEPTED)")
        print("="*30)
        print(f"{'Distance':<10} | {'Count':<10} | {'% of Accepted':<12}")
        print("-" * 38)
        dist_counts = Counter(bc_distances)
        total_accepted = sum(dist_counts.values())
        for dist in sorted(dist_counts.keys()):
            count = dist_counts[dist]
            pct = (count / total_accepted) * 100
            print(f"{dist:<10} | {count:<10} | {pct:<12.1f}%")
            
    total_seconds = time.time() - start_time
    print(f"Total Time: {total_seconds:.2f} seconds")

    total_reads = stats.get('total', 0)
    unassigned = stats.get('unassigned', 0)
    umi_mismatch = stats.get('umi_mismatch', 0)
    assigned = total_reads - unassigned
    percent_assigned = (assigned / total_reads * 100) if total_reads > 0 else 0.0

    print("\n" + "="*30)
    print("FINAL COUNTS")
    print("="*30)
    print(f"{'Category/Sample':<25} | {'Count':<10} | {'% of Total':<10}")
    print("-" * 52)
    print(f"{'Total Reads':<25} | {total_reads:<10} | {100.0:<10.1f}%")
    print(f"{'Assigned Reads':<25} | {assigned:<10} | {percent_assigned:<10.1f}%")
    print("-" * 52)
    for k, v in sorted(stats.items()):
        if k not in ('total', 'unassigned', 'umi_mismatch'):
            pct = (v / total_reads * 100) if total_reads > 0 else 0.0
            print(f"{k:<25} | {v:<10} | {pct:<10.1f}%")
    pct_umi_mismatch = (umi_mismatch / total_reads * 100) if total_reads > 0 else 0.0
    print(f"{'umi_mismatch':<25} | {umi_mismatch:<10} | {pct_umi_mismatch:<10.1f}%")
    pct_unassigned = (unassigned / total_reads * 100) if total_reads > 0 else 0.0
    print(f"{'unassigned':<25} | {unassigned:<10} | {pct_unassigned:<10.1f}%")

    # Write read length distribution report
    write_read_length_report(read_lengths, prefix, suffix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--barcodes", required=True)
    parser.add_argument("--biosample", required=True)
    parser.add_argument("--mismatches", type=int, default=3)
    parser.add_argument("--max-bc-dist", type=int, default=0)
    parser.add_argument("--prefix", type=str, default="")
    parser.add_argument("--suffix", type=str, default="")
    parser.add_argument("--no-trim", action="store_true", help="Disable trimming of TSO and barcode sequences from accepted reads")
    parser.add_argument("--tso-window", type=int, default=1,
                        help="Search for the TSO in the first N bases of each read (1 = position 0 only)")
    parser.add_argument("--bc-trim-len", type=int, default=0,
                        help="Trim only this many bases of the matched barcode entry from the read end "
                             "(e.g. 10 = sample barcode only, keeping a barcode+primer entry's primer). "
                             "0 = trim the whole entry")
    parser.add_argument("--bc-window", type=int, default=1,
                        help="Search for the barcode in the last N bases of each read (1 = read end only)")
    parser.add_argument("--cpus", type=int, default=mp.cpu_count())
    parser.add_argument("--chunk-size", type=int, default=5000)
    run_demux(parser.parse_args())
