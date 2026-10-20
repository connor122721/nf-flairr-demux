#!/usr/bin/env python
"""
sweep_summary.py - Aggregate audit logs and read length data across parameter sweep.

Parses:
  - audit_{sample}_mm{X}_bc{Y}.csv    (per-read decisions + read lengths)
  - read_length_summary_*_mm{X}_bc{Y}.csv  (pre-computed per-sample stats)

Produces:
  1. sweep_summary.csv          - assignment counts per (mm, bcd, sample)
  2. sweep_read_lengths.csv     - read length stats per (mm, bcd, sample, category)
  3. sweep_summary.html         - interactive report with heatmaps, histograms, boxplots, tables
"""

import argparse
import os
import glob
import re
import json
import math
from collections import defaultdict, Counter


# ============================================================================
# Parsing
# ============================================================================

def parse_audit_log(path):
    """Parse audit log: count decisions AND collect binned read lengths."""
    counts = Counter()
    BIN_SIZE = 25
    # length_bins[category] = Counter of bin_index -> count
    length_bins = defaultdict(Counter)
    # running stats per category
    length_stats = defaultdict(lambda: {
        'n': 0, 'sum': 0, 'sum_sq': 0, 'min': float('inf'), 'max': 0,
        'values': []  # collect for percentiles (capped)
    })
    MAX_VALUES_FOR_PERCENTILES = 500_000  # per category per file

    with open(path, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 10:
                continue
            decision = parts[9]
            counts[decision] += 1
            counts['total'] += 1

            try:
                rlen = int(parts[3])
            except (ValueError, IndexError):
                continue

            cat = 'accepted' if decision == 'ACCEPTED' else 'unassigned'
            length_bins[cat][rlen // BIN_SIZE] += 1

            st = length_stats[cat]
            st['n'] += 1
            st['sum'] += rlen
            st['sum_sq'] += rlen * rlen
            if rlen < st['min']:
                st['min'] = rlen
            if rlen > st['max']:
                st['max'] = rlen
            if st['n'] <= MAX_VALUES_FOR_PERCENTILES:
                st['values'].append(rlen)

    return counts, length_bins, length_stats


def extract_audit_params(filename):
    """Extract mm, bc, sample from audit_{sample}_mm{X}_bc{Y}.csv"""
    m = re.match(r'audit_(.+)_mm(\d+)_bc(\d+)\.csv', filename)
    if m:
        return int(m.group(2)), int(m.group(3)), m.group(1)
    return None, None, None


def parse_biosample(csv_path):
    """Map barcode -> sample name from biosample CSV.

    Columns: SampleName,BarcodeName[,Bam]. With the optional Bam column the
    key is (bam, barcode), so the same barcode can mean different samples in
    different BAMs; without it the key is just the barcode.
    """
    mapping = {}
    if not csv_path or not os.path.exists(csv_path):
        return mapping
    with open(csv_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('SampleName'):
                continue
            parts = [p.strip() for p in line.strip().split(',')]
            if len(parts) < 2:
                continue
            bc = parts[1].split('--')[-1]   # 'RACE_TSO_5p--bc1001' -> 'bc1001'
            if len(parts) >= 3 and parts[2]:
                mapping[(parts[2], bc)] = parts[0]
            else:
                mapping[bc] = parts[0]
    return mapping


def map_sample(sample_map, bam_name, bc):
    """Sample name for barcode `bc` in BAM `bam_name` (falls back to the barcode)."""
    return sample_map.get((bam_name, bc), sample_map.get(bc, bc))


def parse_audit_log_for_umi(path):
    """Yield (umi, bc_name) pairs for ACCEPTED reads."""
    with open(path, 'r') as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip('\n').split(',')
            if len(parts) < 10 or parts[9] != 'ACCEPTED':
                continue
            umi = parts[1]
            bc = parts[6]
            if not umi or not bc or bc in ('N/A', 'AMBIGUOUS'):
                continue
            yield umi, bc


def compute_umi_sharing(audit_files_by_combo, sample_map, output_prefix,
                       max_share_bin=6):
    """For each (mm, bcd) combo, compute UMI-level sample sharing stats and
    emit a per-combo CSV compatible with filter_shared_umis.py.

    Sample identifier is ``<bam>__<biosample>`` so that the same barcode
    re-used across different input BAMs (pools) is treated as a distinct
    sample. Without this, UMIs from bc1001 in pool A and bc1001 in pool B
    would collapse into one sample and look identical.

    Returns: dict[(mm, bcd)] -> {
        'total_unique_umis', 'sample_count_dist', 'shared_umis',
        'per_sample': {sample -> {'total', 'shared', 'unique'}},
        'samples': sorted list,
        'top_pairs': [(sample_a, sample_b, count), ...]
    }
    """
    out = {}
    for key, files in sorted(audit_files_by_combo.items()):
        mm, bcd = key
        # umi -> Counter(sample -> read_count)
        umi_to_samples = defaultdict(Counter)
        for fp in files:
            fname = os.path.basename(fp)
            _, _, bam_name = extract_audit_params(fname)
            for umi, bc in parse_audit_log_for_umi(fp):
                mapped = map_sample(sample_map, bam_name, bc)
                # <bam>__<biosample>: BAM names contain single underscores but
                # not double underscores, so this is splittable by FILTER_SHARED_UMIS.
                sample = f"{bam_name}__{mapped}"
                umi_to_samples[umi][sample] += 1

        total_umis = len(umi_to_samples)
        sample_count_dist = Counter()
        samples_seen = set()
        per_sample_total = Counter()
        per_sample_shared = Counter()
        pair_counts = Counter()

        for umi, counts in umi_to_samples.items():
            n_samples = len(counts)
            sample_count_dist[n_samples] += 1
            samples_seen.update(counts.keys())
            for s in counts:
                per_sample_total[s] += 1
                if n_samples >= 2:
                    per_sample_shared[s] += 1
            if n_samples >= 2:
                samples_in_umi = sorted(counts.keys())
                for i in range(len(samples_in_umi)):
                    for j in range(i + 1, len(samples_in_umi)):
                        pair_counts[(samples_in_umi[i], samples_in_umi[j])] += 1

        # ---- per-combo CSV: only UMIs shared by >=2 samples ----
        sample_list = sorted(samples_seen)
        csv_path = f"barcode_sharing_summary_mm{mm}_bc{bcd}.csv"
        with open(csv_path, 'w') as f:
            f.write("umi_combo," +
                    ",".join(f"count_{s}" for s in sample_list) +
                    ",max_count,next_max_count\n")
            for umi, counts in umi_to_samples.items():
                if len(counts) < 2:
                    continue
                vals = [counts.get(s, 0) for s in sample_list]
                sorted_v = sorted(vals, reverse=True)
                mc = sorted_v[0]
                nmc = sorted_v[1] if len(sorted_v) > 1 else 0
                f.write(f"{umi}," + ",".join(str(v) for v in vals) +
                        f",{mc},{nmc}\n")

        per_sample = {}
        for s in sample_list:
            tot = per_sample_total[s]
            shr = per_sample_shared[s]
            per_sample[s] = {
                'total': tot,
                'shared': shr,
                'unique': tot - shr,
                'pct_shared': (shr / tot * 100) if tot > 0 else 0.0,
            }

        top_pairs = sorted(pair_counts.items(), key=lambda kv: -kv[1])[:25]

        out[key] = {
            'total_unique_umis': total_umis,
            'sample_count_dist': dict(sample_count_dist),
            'shared_umis': sum(v for k, v in sample_count_dist.items() if k >= 2),
            'samples': sample_list,
            'per_sample': per_sample,
            'top_pairs': [(a, b, c) for (a, b), c in top_pairs],
            'csv_path': csv_path,
        }
        # Free memory before processing next combo
        del umi_to_samples

    # ---- aggregate CSV across all combos ----
    agg_csv = f"{output_prefix}_umi_sharing.csv"
    with open(agg_csv, 'w') as f:
        f.write("mismatches,max_bc_dist,n_samples,total_unique_umis,"
                "umis_unique_to_one_sample,umis_shared_2plus,pct_shared\n")
        for key in sorted(out.keys()):
            mm, bcd = key
            s = out[key]
            n_one = s['sample_count_dist'].get(1, 0)
            shared = s['shared_umis']
            total = s['total_unique_umis']
            pct = (shared / total * 100) if total > 0 else 0.0
            f.write(f"{mm},{bcd},{len(s['samples'])},{total},"
                    f"{n_one},{shared},{pct:.2f}\n")

    return out, agg_csv


def extract_rl_params(filename):
    """Extract mm, bc, sample from read_length_summary_{sample}_mm{X}_bc{Y}.csv"""
    m = re.match(r'read_length_summary_(.+)_mm(\d+)_bc(\d+)\.csv', filename)
    if m:
        return int(m.group(2)), int(m.group(3)), m.group(1)
    return None, None, None


def parse_read_length_csv(path):
    """Parse a read_length_summary CSV into list of row dicts."""
    rows = []
    with open(path, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 7:
                continue
            rows.append({
                'category': parts[0],
                'sample_bc': parts[1],
                'count': int(parts[2]),
                'mean': float(parts[3]),
                'median': int(parts[4]),
                'min_len': int(parts[5]),
                'max_len': int(parts[6]),
            })
    return rows


def compute_percentiles(values, percentiles=[25, 50, 75]):
    """Compute percentiles from a list of values."""
    if not values:
        return {p: 0 for p in percentiles}
    s = sorted(values)
    n = len(s)
    result = {}
    for p in percentiles:
        k = (n - 1) * p / 100
        f = int(k)
        c = f + 1 if f + 1 < n else f
        result[p] = s[f] + (k - f) * (s[c] - s[f])
    return result


def percentiles_from_bins(bin_counts, bin_size=25, percentiles=[25, 50, 75]):
    """Approximate percentiles from binned counts."""
    if not bin_counts:
        return {p: 0 for p in percentiles}
    total = sum(bin_counts.values())
    if total == 0:
        return {p: 0 for p in percentiles}
    sorted_bins = sorted(bin_counts.items())
    result = {}
    for p in percentiles:
        target = total * p / 100
        cumulative = 0
        for bin_idx, count in sorted_bins:
            cumulative += count
            if cumulative >= target:
                result[p] = bin_idx * bin_size + bin_size // 2
                break
        else:
            result[p] = sorted_bins[-1][0] * bin_size + bin_size // 2
    return result


# ============================================================================
# Main
# ============================================================================

def main(args):
    BIN_SIZE = 25

    # ---- Find files ----
    audit_files = glob.glob(os.path.join(args.input_dir, 'audit_*_mm*_bc*.csv'))
    rl_files = glob.glob(os.path.join(args.input_dir, 'read_length_summary_*_mm*_bc*.csv'))

    if not audit_files:
        print("No audit log files found!")
        open(f"{args.output_prefix}.csv", 'w').write(
            "mismatches,max_bc_dist,sample,total,accepted,pct_accepted\n")
        open(f"{args.output_prefix}.html", 'w').write(
            "<html><body>No data found</body></html>")
        return

    print(f"Found {len(audit_files)} audit logs, {len(rl_files)} read length summaries")

    # ---- Parse audit logs ----
    audit_rows = []
    # Aggregate length data per (mm, bcd, category)
    combo_bins = defaultdict(lambda: defaultdict(Counter))      # (mm,bcd) -> cat -> bins
    combo_stats = defaultdict(lambda: defaultdict(lambda: {
        'n': 0, 'sum': 0, 'min': float('inf'), 'max': 0, 'values': []
    }))
    MAX_VALS = 200_000  # per combo per category

    for csv_path in sorted(audit_files):
        fname = os.path.basename(csv_path)
        mm, bcd, sample = extract_audit_params(fname)
        if mm is None:
            continue

        print(f"  Parsing {fname}...")
        counts, length_bins, length_stats_file = parse_audit_log(csv_path)
        total = counts.get('total', 0)
        accepted = counts.get('ACCEPTED', 0)
        pct = (accepted / total * 100) if total > 0 else 0.0

        audit_rows.append({
            'mm': mm, 'bcd': bcd, 'sample': sample,
            'total': total, 'accepted': accepted, 'pct': pct,
            'rej_tso':     counts.get('REJECTED - TSO NOT FOUND', 0),
            'rej_bc_dist': counts.get('REJECTED - BC DISTANCE TOO HIGH', 0),
            'rej_ambig':   counts.get('REJECTED - AMBIGUOUS BC', 0),
            'rej_umi':     counts.get('REJECTED - UMI MISMATCH', 0),
            'rej_short':   counts.get('REJECTED - READ TOO SHORT', 0),
            'rej_long':    counts.get('REJECTED - READ TOO LONG', 0),
        })

        # Merge into combo-level aggregates
        key = (mm, bcd)
        for cat in ['accepted', 'unassigned']:
            for bin_idx, cnt in length_bins[cat].items():
                combo_bins[key][cat][bin_idx] += cnt
            st_file = length_stats_file[cat]
            st_combo = combo_stats[key][cat]
            st_combo['n'] += st_file['n']
            st_combo['sum'] += st_file['sum']
            if st_file['min'] < st_combo['min']:
                st_combo['min'] = st_file['min']
            if st_file['max'] > st_combo['max']:
                st_combo['max'] = st_file['max']
            if len(st_combo['values']) < MAX_VALS:
                remaining = MAX_VALS - len(st_combo['values'])
                st_combo['values'].extend(st_file['values'][:remaining])

    # ---- Parse read length summaries ----
    rl_rows = []
    for csv_path in sorted(rl_files):
        fname = os.path.basename(csv_path)
        mm, bcd, sample = extract_rl_params(fname)
        if mm is None:
            continue
        for r in parse_read_length_csv(csv_path):
            if r['sample_bc'] in ('ALL_ASSIGNED', 'ALL_UNASSIGNED'):
                continue
            r['mm'] = mm
            r['bcd'] = bcd
            r['sample'] = sample
            rl_rows.append(r)

    # ---- UMI sharing across samples (per parameter combo) ----
    sample_map = parse_biosample(args.biosample)
    if sample_map:
        print(f"Loaded biosample map: {len(sample_map)} barcode->sample entries")
    else:
        print("No biosample mapping provided; using barcode names as sample IDs")

    audit_files_by_combo = defaultdict(list)
    for csv_path in audit_files:
        fname = os.path.basename(csv_path)
        mm, bcd, _ = extract_audit_params(fname)
        if mm is not None:
            audit_files_by_combo[(mm, bcd)].append(csv_path)

    print(f"Computing UMI sharing for {len(audit_files_by_combo)} parameter combos...")
    umi_sharing, umi_csv_path = compute_umi_sharing(
        audit_files_by_combo, sample_map, args.output_prefix
    )

    # ---- Derived data ----
    all_mm = sorted(set(r['mm'] for r in audit_rows))
    all_bcd = sorted(set(r['bcd'] for r in audit_rows))
    all_samples = sorted(set(r['sample'] for r in audit_rows))
    colors6 = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f', '#edc948']
    colors15 = [
        '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
        '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac',
        '#86bcb6', '#d37295', '#8cd17d', '#b6992d', '#499894'
    ]

    # Assignment aggregates
    agg = defaultdict(lambda: {'total': 0, 'accepted': 0})
    for r in audit_rows:
        agg[(r['mm'], r['bcd'])]['total'] += r['total']
        agg[(r['mm'], r['bcd'])]['accepted'] += r['accepted']

    # Compute combo-level length stats
    combo_length_summary = {}
    for key in combo_bins:
        combo_length_summary[key] = {}
        for cat in ['accepted', 'unassigned']:
            st = combo_stats[key][cat]
            n = st['n']
            if n == 0:
                combo_length_summary[key][cat] = {
                    'n': 0, 'mean': 0, 'median': 0, 'q1': 0, 'q3': 0,
                    'min': 0, 'max': 0, 'iqr': 0
                }
                continue
            mean = st['sum'] / n
            pcts = compute_percentiles(st['values'], [25, 50, 75]) if st['values'] else \
                   percentiles_from_bins(combo_bins[key][cat], BIN_SIZE, [25, 50, 75])
            iqr = pcts[75] - pcts[25]
            combo_length_summary[key][cat] = {
                'n': n, 'mean': round(mean, 1),
                'median': round(pcts[50]), 'q1': round(pcts[25]), 'q3': round(pcts[75]),
                'min': st['min'], 'max': st['max'], 'iqr': round(iqr),
            }

    # ---- Write CSVs ----
    csv_path = f"{args.output_prefix}.csv"
    with open(csv_path, 'w') as f:
        f.write("mismatches,max_bc_dist,sample,total,accepted,pct_accepted,"
                "rej_tso,rej_bc_dist,rej_ambiguous,rej_umi,rej_short,rej_long\n")
        for r in sorted(audit_rows, key=lambda x: (x['mm'], x['bcd'], x['sample'])):
            f.write(f"{r['mm']},{r['bcd']},{r['sample']},{r['total']},{r['accepted']},{r['pct']:.2f},"
                    f"{r['rej_tso']},{r['rej_bc_dist']},{r['rej_ambig']},{r['rej_umi']},"
                    f"{r['rej_short']},{r['rej_long']}\n")

    rl_csv_path = f"{args.output_prefix}_read_lengths.csv"
    with open(rl_csv_path, 'w') as f:
        f.write("mismatches,max_bc_dist,category,"
                "n,mean,q1,median,q3,min,max,iqr\n")
        for key in sorted(combo_length_summary.keys()):
            mm, bcd = key
            for cat in ['accepted', 'unassigned']:
                s = combo_length_summary[key][cat]
                if s['n'] == 0:
                    continue
                f.write(f"{mm},{bcd},{cat},"
                        f"{s['n']},{s['mean']},{s['q1']},{s['median']},{s['q3']},"
                        f"{s['min']},{s['max']},{s['iqr']}\n")

    # Per-sample read length table
    rl_detail_csv = f"{args.output_prefix}_read_lengths_detail.csv"
    with open(rl_detail_csv, 'w') as f:
        f.write("mismatches,max_bc_dist,sample,category,barcode,count,mean,median,min,max\n")
        for r in sorted(rl_rows, key=lambda x: (x['mm'], x['bcd'], x['sample'], x['category'])):
            f.write(f"{r['mm']},{r['bcd']},{r['sample']},{r['category']},{r['sample_bc']},"
                    f"{r['count']},{r['mean']:.1f},{r['median']},{r['min_len']},{r['max_len']}\n")

    # ---- Write HTML ----
    html_path = f"{args.output_prefix}.html"
    with open(html_path, 'w') as f:

        # === Header ===
        f.write(f"""<!DOCTYPE html>
<html>
<head>
<title>Parameter Sweep Summary</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; padding: 20px 30px; background: #f8f9fa; color: #333; }}
    h1 {{ color: #1a1a2e; border-bottom: 3px solid #4e79a7; padding-bottom: 8px; }}
    h2 {{ color: #16213e; margin-top: 40px; }}
    h3 {{ color: #4e79a7; }}
    .plot-container {{ width: 100%; max-width: 1100px; height: 520px; margin: 15px 0 30px 0; }}
    .plot-wide {{ height: 600px; }}
    .info {{ background: #e8f0fe; padding: 14px 18px; border-radius: 8px; margin: 10px 0 25px 0;
             border-left: 4px solid #4e79a7; }}
    .nav {{ background: #1a1a2e; padding: 12px 20px; border-radius: 8px; margin-bottom: 25px;
            position: sticky; top: 0; z-index: 100; }}
    .nav a {{ color: #a8d8ea; text-decoration: none; margin-right: 18px; font-size: 14px; }}
    .nav a:hover {{ color: #fff; text-decoration: underline; }}
    table {{ border-collapse: collapse; margin: 15px 0 30px 0; width: 100%; max-width: 1100px; font-size: 13px; }}
    th, td {{ border: 1px solid #dee2e6; padding: 7px 10px; text-align: right; }}
    th {{ background: #4e79a7; color: white; position: sticky; top: 44px; }}
    tr:nth-child(even) {{ background: #f2f4f8; }}
    tr:hover {{ background: #dce6f0; }}
    td:first-child, td:nth-child(2), td:nth-child(3) {{ text-align: left; }}
    .section {{ margin-bottom: 50px; }}
    .tab-container {{ display: flex; gap: 8px; margin-bottom: -1px; }}
    .tab {{ padding: 8px 16px; cursor: pointer; background: #dee2e6; border: 1px solid #ccc;
            border-bottom: none; border-radius: 6px 6px 0 0; font-size: 13px; }}
    .tab.active {{ background: #fff; border-bottom: 1px solid #fff; font-weight: bold; }}
    .tab-content {{ display: none; border: 1px solid #ccc; padding: 15px; background: #fff; border-radius: 0 6px 6px 6px; }}
    .tab-content.active {{ display: block; }}
</style>
</head>
<body>

<h1>Parameter Sweep Report</h1>

<div class="nav">
    <a href="#sec-assignment">Assignment Rates</a>
    <a href="#sec-rejections">Rejections</a>
    <a href="#sec-readlen-heatmaps">Read Length Heatmaps</a>
    <a href="#sec-histograms">Histograms</a>
    <a href="#sec-boxplots">Boxplots</a>
    <a href="#sec-umi-sharing">UMI Sharing</a>
    <a href="#sec-tables">Data Tables</a>
</div>

<div class="info">
    <strong>Sweep grid:</strong> mismatches [{min(all_mm)}&ndash;{max(all_mm)}] &times;
    max-bc-dist [{min(all_bcd)}&ndash;{max(all_bcd)}] = {len(all_mm) * len(all_bcd)} combos &times;
    {len(all_samples)} samples = {len(audit_rows)} total runs
</div>
""")

        # =================================================================
        # SECTION 1: Assignment Rate
        # =================================================================
        z_pct, z_text = [], []
        for bcd in all_bcd:
            row_p, row_t = [], []
            for mm in all_mm:
                d = agg.get((mm, bcd), {'total': 0, 'accepted': 0})
                pct = (d['accepted'] / d['total'] * 100) if d['total'] > 0 else 0
                row_p.append(round(pct, 2))
                row_t.append(f"{d['accepted']:,} / {d['total']:,}<br>{pct:.1f}%")
            z_pct.append(row_p)
            z_text.append(row_t)

        f.write(f"""
<div class="section" id="sec-assignment">
<h2>Assignment Rate (% Accepted)</h2>
<div id="heatmap_accept" class="plot-container"></div>
<script>
Plotly.newPlot('heatmap_accept', [{{
    z: {json.dumps(z_pct)}, x: {json.dumps(all_mm)},
    y: {json.dumps(['bc ' + str(b) for b in all_bcd])},
    text: {json.dumps(z_text)}, hoverinfo: 'text',
    type: 'heatmap', colorscale: 'YlGnBu',
    colorbar: {{ title: '% Accepted' }}, zmin: 0, zmax: 100
}}], {{
    title: 'Assignment Rate by Parameter Combination',
    xaxis: {{ title: 'Mismatches (TSO)', dtick: 1 }},
    yaxis: {{ title: 'Max Barcode Distance' }},
    template: 'plotly_white'
}}, {{responsive: true}});
</script>

<h3>Assignment Rate vs Mismatches (by Max BC Distance)</h3>
<div id="line_chart" class="plot-container"></div>
<script>
var lineTraces = [
""")
        for i, bcd in enumerate(all_bcd):
            pcts = []
            for mm in all_mm:
                d = agg.get((mm, bcd), {'total': 0, 'accepted': 0})
                pct = (d['accepted'] / d['total'] * 100) if d['total'] > 0 else 0
                pcts.append(round(pct, 2))
            f.write(f"  {{ x: {json.dumps(all_mm)}, y: {json.dumps(pcts)}, "
                    f"name: 'bc_dist={bcd}', type: 'scatter', mode: 'lines+markers', "
                    f"marker: {{ color: '{colors6[i % len(colors6)]}' }} }},\n")

        f.write(f"""];
Plotly.newPlot('line_chart', lineTraces, {{
    xaxis: {{ title: 'Mismatches (TSO)', dtick: 1 }},
    yaxis: {{ title: '% Accepted', range: [0, 100] }},
    template: 'plotly_white',
    legend: {{ title: {{ text: 'Max BC Distance' }} }}
}}, {{responsive: true}});
</script>
</div>
""")

        # =================================================================
        # SECTION 2: Rejection heatmaps
        # =================================================================
        rej_types = {
            'rej_tso': 'TSO Not Found',
            'rej_bc_dist': 'BC Distance Too High',
            'rej_ambig': 'Ambiguous BC',
            'rej_umi': 'UMI Mismatch',
        }

        f.write('<div class="section" id="sec-rejections">\n<h2>Rejection Breakdown</h2>\n')

        for idx, (rej_key, rej_label) in enumerate(rej_types.items()):
            rej_agg = defaultdict(int)
            total_agg = defaultdict(int)
            for r in audit_rows:
                k = (r['mm'], r['bcd'])
                rej_agg[k] += r[rej_key]
                total_agg[k] += r['total']

            z_rej, z_rej_t = [], []
            for bcd in all_bcd:
                row, row_t = [], []
                for mm in all_mm:
                    k = (mm, bcd)
                    cnt = rej_agg.get(k, 0)
                    tot = total_agg.get(k, 1)
                    pct = cnt / tot * 100 if tot > 0 else 0
                    row.append(round(pct, 2))
                    row_t.append(f"{cnt:,}<br>{pct:.1f}%")
                z_rej.append(row)
                z_rej_t.append(row_t)

            div_id = f"heatmap_rej_{idx}"
            f.write(f"""
<h3>Rejection: {rej_label}</h3>
<div id="{div_id}" class="plot-container"></div>
<script>
Plotly.newPlot('{div_id}', [{{
    z: {json.dumps(z_rej)}, x: {json.dumps(all_mm)},
    y: {json.dumps(['bc ' + str(b) for b in all_bcd])},
    text: {json.dumps(z_rej_t)}, hoverinfo: 'text',
    type: 'heatmap', colorscale: 'OrRd',
    colorbar: {{ title: '% Rejected' }}
}}], {{
    title: '{rej_label} Rate',
    xaxis: {{ title: 'Mismatches (TSO)', dtick: 1 }},
    yaxis: {{ title: 'Max Barcode Distance' }},
    template: 'plotly_white'
}}, {{responsive: true}});
</script>
""")
        f.write('</div>\n')

        # =================================================================
        # SECTION 3: Read Length Heatmaps (mean & median for accepted)
        # =================================================================
        f.write('<div class="section" id="sec-readlen-heatmaps">\n'
                '<h2>Read Length Heatmaps (Accepted Reads)</h2>\n')

        for metric, label in [('mean', 'Mean Read Length'), ('median', 'Median Read Length')]:
            z_len, z_len_t = [], []
            for bcd in all_bcd:
                row, row_t = [], []
                for mm in all_mm:
                    s = combo_length_summary.get((mm, bcd), {}).get('accepted', {})
                    val = s.get(metric, 0)
                    n = s.get('n', 0)
                    row.append(round(val, 1) if val else 0)
                    row_t.append(f"n={n:,}<br>{metric}={val:,.0f} bp")
                z_len.append(row)
                z_len_t.append(row_t)

            div_id = f"heatmap_rl_{metric}"
            f.write(f"""
<h3>{label} (Accepted)</h3>
<div id="{div_id}" class="plot-container"></div>
<script>
Plotly.newPlot('{div_id}', [{{
    z: {json.dumps(z_len)}, x: {json.dumps(all_mm)},
    y: {json.dumps(['bc ' + str(b) for b in all_bcd])},
    text: {json.dumps(z_len_t)}, hoverinfo: 'text',
    type: 'heatmap', colorscale: 'Viridis',
    colorbar: {{ title: 'bp' }}
}}], {{
    title: '{label} by Parameter Combination',
    xaxis: {{ title: 'Mismatches (TSO)', dtick: 1 }},
    yaxis: {{ title: 'Max Barcode Distance' }},
    template: 'plotly_white'
}}, {{responsive: true}});
</script>
""")
        f.write('</div>\n')

        # =================================================================
        # SECTION 4: Read Length Histograms
        # =================================================================
        f.write('<div class="section" id="sec-histograms">\n'
                '<h2>Read Length Histograms</h2>\n'
                '<p>Average read length distributions across all samples. '
                'Use the dropdown to select parameter combinations.</p>\n')

        # Build histogram traces for a selection of combos + the "all combos averaged"
        # Aggregate all accepted bins across all combos
        total_accepted_bins = Counter()
        total_unassigned_bins = Counter()
        for key in combo_bins:
            for bin_idx, cnt in combo_bins[key].get('accepted', {}).items():
                total_accepted_bins[bin_idx] += cnt
            for bin_idx, cnt in combo_bins[key].get('unassigned', {}).items():
                total_unassigned_bins[bin_idx] += cnt

        def bins_to_xy(bin_counts, bin_size=BIN_SIZE):
            if not bin_counts:
                return [], []
            x = [b * bin_size + bin_size // 2 for b in sorted(bin_counts.keys())]
            y = [bin_counts[b] for b in sorted(bin_counts.keys())]
            return x, y

        # Global histogram
        x_acc, y_acc = bins_to_xy(total_accepted_bins)
        x_un, y_un = bins_to_xy(total_unassigned_bins)

        f.write(f"""
<h3>All Combos Aggregated</h3>
<div id="hist_global" class="plot-container"></div>
<script>
Plotly.newPlot('hist_global', [
    {{ x: {json.dumps(x_acc)}, y: {json.dumps(y_acc)}, name: 'Accepted', type: 'bar',
       marker: {{ color: '#4e79a7' }}, opacity: 0.7 }},
    {{ x: {json.dumps(x_un)}, y: {json.dumps(y_un)}, name: 'Unassigned', type: 'bar',
       marker: {{ color: '#e15759' }}, opacity: 0.7 }}
], {{
    barmode: 'overlay',
    title: 'Read Length Distribution (All Combos)',
    xaxis: {{ title: 'Read Length (bp)' }},
    yaxis: {{ title: 'Count' }},
    template: 'plotly_white'
}}, {{responsive: true}});
</script>
""")

        # Per-bc_dist histograms (accepted only, one line per bc_dist, aggregated across mm)
        f.write("""
<h3>Accepted Read Lengths by Max BC Distance (aggregated across mismatches)</h3>
<div id="hist_by_bcd" class="plot-container"></div>
<script>
var histBcdTraces = [
""")
        for i, bcd in enumerate(all_bcd):
            bcd_bins = Counter()
            for mm in all_mm:
                for bin_idx, cnt in combo_bins.get((mm, bcd), {}).get('accepted', {}).items():
                    bcd_bins[bin_idx] += cnt
            x, y = bins_to_xy(bcd_bins)
            # Normalize to density for comparison
            total_count = sum(y) if y else 1
            y_norm = [round(v / total_count * 100, 3) for v in y]
            f.write(f"  {{ x: {json.dumps(x)}, y: {json.dumps(y_norm)}, "
                    f"name: 'bc_dist={bcd} (n={total_count:,})', type: 'scatter', "
                    f"mode: 'lines', fill: 'tozeroy', opacity: 0.5, "
                    f"line: {{ color: '{colors6[i % len(colors6)]}' }} }},\n")

        f.write(f"""];
Plotly.newPlot('hist_by_bcd', histBcdTraces, {{
    title: 'Accepted Read Length Density by Max BC Distance',
    xaxis: {{ title: 'Read Length (bp)' }},
    yaxis: {{ title: '% of Reads' }},
    template: 'plotly_white',
    legend: {{ title: {{ text: 'Max BC Distance' }} }}
}}, {{responsive: true}});
</script>
""")

        # Per-mismatches histograms (accepted, one line per mm, aggregated across bc_dist)
        f.write("""
<h3>Accepted Read Lengths by Mismatches (aggregated across bc_dist)</h3>
<div id="hist_by_mm" class="plot-container plot-wide"></div>
<script>
var histMmTraces = [
""")
        for i, mm in enumerate(all_mm):
            mm_bins = Counter()
            for bcd in all_bcd:
                for bin_idx, cnt in combo_bins.get((mm, bcd), {}).get('accepted', {}).items():
                    mm_bins[bin_idx] += cnt
            x, y = bins_to_xy(mm_bins)
            total_count = sum(y) if y else 1
            y_norm = [round(v / total_count * 100, 3) for v in y]
            vis = 'true' if mm <= 5 else "'legendonly'"
            f.write(f"  {{ x: {json.dumps(x)}, y: {json.dumps(y_norm)}, "
                    f"name: 'mm={mm} (n={total_count:,})', type: 'scatter', "
                    f"mode: 'lines', visible: {vis}, "
                    f"line: {{ color: '{colors15[i % len(colors15)]}' }} }},\n")

        f.write(f"""];
Plotly.newPlot('hist_by_mm', histMmTraces, {{
    title: 'Accepted Read Length Density by Mismatches',
    xaxis: {{ title: 'Read Length (bp)' }},
    yaxis: {{ title: '% of Reads' }},
    template: 'plotly_white',
    legend: {{ title: {{ text: 'Mismatches' }} }}
}}, {{responsive: true}});
</script>
""")
        f.write('</div>\n')

        # =================================================================
        # SECTION 5: Boxplots
        # =================================================================
        f.write('<div class="section" id="sec-boxplots">\n'
                '<h2>Read Length Boxplots</h2>\n')

        # Boxplot: accepted read lengths, one box per (mm, bcd) combo
        # Use Plotly's pre-computed stat boxplots
        f.write("""
<h3>Accepted Read Lengths by Parameter Combination</h3>
<p>Each box represents one (mismatches, bc_dist) combination aggregated across all samples.
   Use the dropdown to group by mismatches or by bc_dist.</p>
<div id="boxplot_combo" class="plot-container plot-wide"></div>
<script>
var boxTraces = [
""")
        for i, bcd in enumerate(all_bcd):
            x_labels = []
            q1s, medians, q3s, mins, maxs, means = [], [], [], [], [], []
            for mm in all_mm:
                s = combo_length_summary.get((mm, bcd), {}).get('accepted', {})
                if not s or s.get('n', 0) == 0:
                    continue
                x_labels.append(f"mm{mm}")
                q1s.append(s['q1'])
                medians.append(s['median'])
                q3s.append(s['q3'])
                mins.append(s['min'])
                maxs.append(s['max'])
                means.append(s['mean'])

            f.write(f"""  {{
    type: 'box', name: 'bc_dist={bcd}',
    x: {json.dumps(x_labels)},
    lowerfence: {json.dumps(mins)}, q1: {json.dumps(q1s)},
    median: {json.dumps(medians)}, q3: {json.dumps(q3s)},
    upperfence: {json.dumps(maxs)}, mean: {json.dumps(means)},
    marker: {{ color: '{colors6[i % len(colors6)]}' }},
    boxmean: true
  }},
""")

        f.write(f"""];
Plotly.newPlot('boxplot_combo', boxTraces, {{
    title: 'Accepted Read Length Distribution',
    xaxis: {{ title: 'Mismatches' }},
    yaxis: {{ title: 'Read Length (bp)' }},
    boxmode: 'group',
    template: 'plotly_white',
    legend: {{ title: {{ text: 'Max BC Distance' }} }}
}}, {{responsive: true}});
</script>
""")

        # Boxplot: accepted vs unassigned side-by-side, grouped by mismatches (bc_dist=0 as default)
        f.write("""
<h3>Accepted vs Unassigned Read Lengths (bc_dist=0)</h3>
<div id="boxplot_vs" class="plot-container plot-wide"></div>
<script>
var boxVsTraces = [
""")
        for cat, color, label in [('accepted', '#4e79a7', 'Accepted'),
                                   ('unassigned', '#e15759', 'Unassigned')]:
            x_labels = []
            q1s, medians, q3s, mins, maxs, means = [], [], [], [], [], []
            for mm in all_mm:
                s = combo_length_summary.get((mm, 0), {}).get(cat, {})
                if not s or s.get('n', 0) == 0:
                    continue
                x_labels.append(f"mm{mm}")
                q1s.append(s['q1'])
                medians.append(s['median'])
                q3s.append(s['q3'])
                mins.append(s['min'])
                maxs.append(s['max'])
                means.append(s['mean'])

            f.write(f"""  {{
    type: 'box', name: '{label}',
    x: {json.dumps(x_labels)},
    lowerfence: {json.dumps(mins)}, q1: {json.dumps(q1s)},
    median: {json.dumps(medians)}, q3: {json.dumps(q3s)},
    upperfence: {json.dumps(maxs)}, mean: {json.dumps(means)},
    marker: {{ color: '{color}' }},
    boxmean: true
  }},
""")

        f.write(f"""];
Plotly.newPlot('boxplot_vs', boxVsTraces, {{
    title: 'Accepted vs Unassigned Read Lengths (bc_dist=0)',
    xaxis: {{ title: 'Mismatches' }},
    yaxis: {{ title: 'Read Length (bp)' }},
    boxmode: 'group',
    template: 'plotly_white'
}}, {{responsive: true}});
</script>
</div>
""")

        # =================================================================
        # SECTION 6: UMI Sharing Across Samples
        # =================================================================
        f.write('<div class="section" id="sec-umi-sharing">\n'
                '<h2>UMI Sharing Across Samples</h2>\n'
                '<p>For each parameter combination, the number of unique UMIs observed '
                'in <strong>accepted</strong> reads is broken down by how many samples each '
                'UMI is seen in. UMIs appearing in &ge; 2 samples may indicate '
                'cross-sample contamination, barcode hopping, or PCR jackpotting. '
                'Per-combo CSVs (<code>barcode_sharing_summary_mm{X}_bc{Y}.csv</code>) '
                'are emitted in the same format consumed by '
                '<code>bin/filter_shared_umis.py</code>.</p>\n')

        # --- Heatmap: % of UMIs shared across 2+ samples per (mm, bcd) ---
        z_share, z_share_t = [], []
        for bcd in all_bcd:
            row, row_t = [], []
            for mm in all_mm:
                s = umi_sharing.get((mm, bcd))
                if not s or s['total_unique_umis'] == 0:
                    row.append(0)
                    row_t.append('no data')
                    continue
                pct = s['shared_umis'] / s['total_unique_umis'] * 100
                row.append(round(pct, 2))
                row_t.append(f"{s['shared_umis']:,} shared / "
                             f"{s['total_unique_umis']:,} total<br>{pct:.1f}%")
            z_share.append(row)
            z_share_t.append(row_t)

        f.write(f"""
<h3>% UMIs Shared Across 2+ Samples</h3>
<div id="heatmap_umi_share" class="plot-container"></div>
<script>
Plotly.newPlot('heatmap_umi_share', [{{
    z: {json.dumps(z_share)}, x: {json.dumps(all_mm)},
    y: {json.dumps(['bc ' + str(b) for b in all_bcd])},
    text: {json.dumps(z_share_t)}, hoverinfo: 'text',
    type: 'heatmap', colorscale: 'OrRd',
    colorbar: {{ title: '% Shared' }}
}}], {{
    title: 'UMI Cross-Sample Sharing Rate',
    xaxis: {{ title: 'Mismatches (TSO)', dtick: 1 }},
    yaxis: {{ title: 'Max Barcode Distance' }},
    template: 'plotly_white'
}}, {{responsive: true}});
</script>
""")

        # --- Distribution: # of UMIs by sample-count ---
        max_share_seen = 1
        for s in umi_sharing.values():
            if s['sample_count_dist']:
                max_share_seen = max(max_share_seen, max(s['sample_count_dist'].keys()))
        # cap bin labels at 6+ for display
        BIN_CAP = 6
        share_x_labels = [str(i) for i in range(1, BIN_CAP)] + [f"{BIN_CAP}+"]

        f.write("""
<h3>Distribution of UMIs by Number of Samples Observed In</h3>
<p>Stacked bar per parameter combo. UMIs appearing in 6 or more samples are pooled.</p>
<div id="bar_umi_dist" class="plot-container plot-wide"></div>
<script>
var umiDistTraces = [
""")
        # one trace per bin (1, 2, 3, 4, 5, 6+)
        stack_colors = ['#4e79a7', '#76b7b2', '#f28e2b', '#e15759', '#b07aa1', '#59a14f']
        combo_labels = [f"mm{mm}_bc{bcd}" for (mm, bcd) in sorted(umi_sharing.keys())]
        for b, label in enumerate(share_x_labels):
            counts_for_bin = []
            for key in sorted(umi_sharing.keys()):
                dist = umi_sharing[key]['sample_count_dist']
                if b < BIN_CAP - 1:
                    counts_for_bin.append(dist.get(b + 1, 0))
                else:
                    counts_for_bin.append(sum(v for k, v in dist.items() if k >= BIN_CAP))
            color = stack_colors[b % len(stack_colors)]
            name = (f"{label} sample" if b == 0
                    else f"{label} samples")
            f.write(f"  {{ x: {json.dumps(combo_labels)}, y: {json.dumps(counts_for_bin)}, "
                    f"name: '{name}', type: 'bar', "
                    f"marker: {{ color: '{color}' }} }},\n")
        f.write(f"""];
Plotly.newPlot('bar_umi_dist', umiDistTraces, {{
    barmode: 'stack',
    title: 'Unique UMIs per Parameter Combo by Sample Multiplicity',
    xaxis: {{ title: 'Parameter Combination', tickangle: -45 }},
    yaxis: {{ title: 'Unique UMIs (count)' }},
    template: 'plotly_white',
    legend: {{ title: {{ text: '# samples per UMI' }} }}
}}, {{responsive: true}});
</script>
""")

        # --- Per-(mm, bcd) summary table ---
        f.write("""
<h3>Per-Combo UMI Sharing Summary</h3>
<table>
<tr><th>Mismatches</th><th>BC Dist</th><th>Samples</th><th>Total Unique UMIs</th>
<th>Unique to 1</th><th>In 2</th><th>In 3</th><th>In 4</th><th>In 5</th>
<th>In 6+</th><th>% Shared (&ge;2)</th></tr>
""")
        for key in sorted(umi_sharing.keys()):
            mm, bcd = key
            s = umi_sharing[key]
            dist = s['sample_count_dist']
            total = s['total_unique_umis']
            shared = s['shared_umis']
            pct = (shared / total * 100) if total > 0 else 0.0
            in_6plus = sum(v for k, v in dist.items() if k >= 6)
            f.write(f"<tr><td>{mm}</td><td>{bcd}</td><td>{len(s['samples'])}</td>"
                    f"<td>{total:,}</td>"
                    f"<td>{dist.get(1, 0):,}</td><td>{dist.get(2, 0):,}</td>"
                    f"<td>{dist.get(3, 0):,}</td><td>{dist.get(4, 0):,}</td>"
                    f"<td>{dist.get(5, 0):,}</td><td>{in_6plus:,}</td>"
                    f"<td>{pct:.1f}%</td></tr>\n")
        f.write("</table>\n")

        # --- Per-sample sharing (default combo = first one) ---
        default_key = sorted(umi_sharing.keys())[0] if umi_sharing else None
        if default_key:
            mm0, bcd0 = default_key
            per_sample = umi_sharing[default_key]['per_sample']
            f.write(f"""
<h3>Per-Sample UMI Sharing (mm={mm0}, bc_dist={bcd0})</h3>
<p>Counts of UMIs assigned to each sample, broken down by whether the UMI is
   unique to that sample or also observed in at least one other sample.</p>
<div style="max-height: 500px; overflow-y: auto;">
<table>
<tr><th>Sample</th><th>UMIs (Total)</th><th>UMIs Unique</th>
<th>UMIs Shared</th><th>% Shared</th></tr>
""")
            for s in sorted(per_sample.keys()):
                d = per_sample[s]
                f.write(f"<tr><td>{s}</td><td>{d['total']:,}</td>"
                        f"<td>{d['unique']:,}</td><td>{d['shared']:,}</td>"
                        f"<td>{d['pct_shared']:.1f}%</td></tr>\n")
            f.write("</table>\n</div>\n")

            # --- Top shared sample pairs ---
            top_pairs = umi_sharing[default_key]['top_pairs']
            if top_pairs:
                f.write(f"""
<h3>Top Cross-Sample UMI Overlaps (mm={mm0}, bc_dist={bcd0})</h3>
<p>Sample pairs with the largest number of co-observed UMIs (top 25).</p>
<table>
<tr><th>Sample A</th><th>Sample B</th><th>Shared UMIs</th></tr>
""")
                for a, b, c in top_pairs:
                    f.write(f"<tr><td>{a}</td><td>{b}</td><td>{c:,}</td></tr>\n")
                f.write("</table>\n")

        f.write('</div>\n')

        # =================================================================
        # SECTION 7: Data Tables
        # =================================================================
        f.write('<div class="section" id="sec-tables">\n<h2>Data Tables</h2>\n')

        # --- Table 1: Assignment summary per combo ---
        f.write("""
<h3>Assignment Summary per Parameter Combination</h3>
<table>
<tr><th>Mismatches</th><th>BC Dist</th><th>Total</th><th>Accepted</th><th>% Accepted</th>
<th>Rej TSO</th><th>Rej BC Dist</th><th>Rej Ambig</th><th>Rej UMI</th>
<th>Rej Short</th><th>Rej Long</th></tr>
""")
        # Aggregate per combo
        combo_agg = defaultdict(lambda: {
            'total': 0, 'accepted': 0,
            'rej_tso': 0, 'rej_bc_dist': 0, 'rej_ambig': 0,
            'rej_umi': 0, 'rej_short': 0, 'rej_long': 0,
        })
        for r in audit_rows:
            k = (r['mm'], r['bcd'])
            for field in combo_agg[k]:
                combo_agg[k][field] += r[field]

        for k in sorted(combo_agg.keys()):
            mm, bcd = k
            d = combo_agg[k]
            pct = (d['accepted'] / d['total'] * 100) if d['total'] > 0 else 0
            f.write(f"<tr><td>{mm}</td><td>{bcd}</td><td>{d['total']:,}</td>"
                    f"<td>{d['accepted']:,}</td><td>{pct:.1f}%</td>"
                    f"<td>{d['rej_tso']:,}</td><td>{d['rej_bc_dist']:,}</td>"
                    f"<td>{d['rej_ambig']:,}</td><td>{d['rej_umi']:,}</td>"
                    f"<td>{d['rej_short']:,}</td><td>{d['rej_long']:,}</td></tr>\n")
        f.write("</table>\n")

        # --- Table 2: Read length stats per combo ---
        f.write("""
<h3>Read Length Statistics per Parameter Combination (Accepted Reads)</h3>
<table>
<tr><th>Mismatches</th><th>BC Dist</th><th>N</th><th>Mean</th>
<th>Q1</th><th>Median</th><th>Q3</th><th>IQR</th>
<th>Min</th><th>Max</th></tr>
""")
        for k in sorted(combo_length_summary.keys()):
            mm, bcd = k
            s = combo_length_summary[k].get('accepted', {})
            if not s or s.get('n', 0) == 0:
                continue
            f.write(f"<tr><td>{mm}</td><td>{bcd}</td><td>{s['n']:,}</td>"
                    f"<td>{s['mean']:,.1f}</td>"
                    f"<td>{s['q1']:,}</td><td>{s['median']:,}</td><td>{s['q3']:,}</td>"
                    f"<td>{s['iqr']:,}</td>"
                    f"<td>{s['min']:,}</td><td>{s['max']:,}</td></tr>\n")
        f.write("</table>\n")

        # --- Table 3: Read length stats per combo (Unassigned) ---
        f.write("""
<h3>Read Length Statistics per Parameter Combination (Unassigned Reads)</h3>
<table>
<tr><th>Mismatches</th><th>BC Dist</th><th>N</th><th>Mean</th>
<th>Q1</th><th>Median</th><th>Q3</th><th>IQR</th>
<th>Min</th><th>Max</th></tr>
""")
        for k in sorted(combo_length_summary.keys()):
            mm, bcd = k
            s = combo_length_summary[k].get('unassigned', {})
            if not s or s.get('n', 0) == 0:
                continue
            f.write(f"<tr><td>{mm}</td><td>{bcd}</td><td>{s['n']:,}</td>"
                    f"<td>{s['mean']:,.1f}</td>"
                    f"<td>{s['q1']:,}</td><td>{s['median']:,}</td><td>{s['q3']:,}</td>"
                    f"<td>{s['iqr']:,}</td>"
                    f"<td>{s['min']:,}</td><td>{s['max']:,}</td></tr>\n")
        f.write("</table>\n")

        # --- Table 4: Per-sample detail (from read_length_summary CSVs) ---
        if rl_rows:
            f.write("""
<h3>Per-Sample Read Length Detail</h3>
<p>Showing per-barcode stats for each sample × parameter combination.
   Full data available in <code>sweep_summary_read_lengths_detail.csv</code>.</p>
<div style="max-height: 600px; overflow-y: auto;">
<table>
<tr><th>MM</th><th>BC Dist</th><th>Sample</th><th>Category</th><th>Barcode</th>
<th>Count</th><th>Mean</th><th>Median</th><th>Min</th><th>Max</th></tr>
""")
            for r in sorted(rl_rows, key=lambda x: (x['mm'], x['bcd'], x['sample'],
                                                     x['category'], x['sample_bc'])):
                f.write(f"<tr><td>{r['mm']}</td><td>{r['bcd']}</td>"
                        f"<td>{r['sample']}</td><td>{r['category']}</td>"
                        f"<td>{r['sample_bc']}</td>"
                        f"<td>{r['count']:,}</td><td>{r['mean']:,.1f}</td>"
                        f"<td>{r['median']:,}</td><td>{r['min_len']:,}</td>"
                        f"<td>{r['max_len']:,}</td></tr>\n")
            f.write("</table>\n</div>\n")

        f.write("</div>\n</body>\n</html>")

    print(f"\nOutputs:")
    print(f"  Assignment CSV:       {csv_path}")
    print(f"  Read length CSV:      {rl_csv_path}")
    print(f"  Read length detail:   {rl_detail_csv}")
    print(f"  UMI sharing CSV:      {umi_csv_path}")
    print(f"  Per-combo UMI CSVs:   barcode_sharing_summary_mm*_bc*.csv "
          f"({len(umi_sharing)} files)")
    print(f"  HTML report:          {html_path}")
    print(f"  Parsed {len(audit_files)} audit logs, {len(rl_files)} read length files")
    print(f"  {len(all_samples)} samples, {len(all_mm) * len(all_bcd)} parameter combos")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', required=True)
    parser.add_argument('--output_prefix', default='sweep_summary')
    parser.add_argument('--biosample', default=None,
                        help='Optional biosample CSV (SampleName,BarcodeName) '
                             'to remap barcodes to sample names in UMI sharing stats')
    args = parser.parse_args()
    main(args)
