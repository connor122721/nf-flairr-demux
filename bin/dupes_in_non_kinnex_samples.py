"""UMI sharing heatmaps from a sweep_summary barcode_sharing_summary CSV.

Consumes one `barcode_sharing_summary_mm{X}_bc{Y}.csv` produced by
`bin/sweep_summary.py` (one row per UMI that appears in 2+ samples, with one
`count_<sample>` column per sample) and emits two heatmaps showing pairwise
shared-UMI counts:

  1. Pool x barcode view  — labels are the full "<bam>__<biosample>" sample
     IDs from the CSV. Reveals same-barcode-across-pools contamination.
  2. Biosample view       — labels are collapsed to the part after the last
     "__". Reveals cross-barcode contamination independent of pool.

Both views are emitted every run. No new dependencies — stdlib + Plotly via
the same CDN URL already used in `bin/sweep_summary.py`.
"""

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter
from itertools import combinations


PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.27.0.min.js"


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def parse_sharing_csv(path):
    """Yield (active_samples_full, header_sample_ids) for every row of the CSV.

    `active_samples_full` is the list of `<bam>__<biosample>` IDs whose
    `count_<id>` column is > 0 on that row. `header_sample_ids` is the full
    ordered list extracted from the header (returned on every row so the
    caller can snapshot it once).
    """
    with open(path, 'r', newline='') as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"Empty CSV: {path}")
        count_cols = [c for c in reader.fieldnames if c.startswith('count_')]
        if not count_cols:
            raise ValueError(
                f"No count_* columns found in {path}; not a barcode_sharing_summary CSV.")
        sample_ids = [c[len('count_'):] for c in count_cols]

        for row in reader:
            active = []
            for sid, col in zip(sample_ids, count_cols):
                raw = row.get(col, '0') or '0'
                try:
                    if int(raw) > 0:
                        active.append(sid)
                except ValueError:
                    continue
            yield active, sample_ids


def parse_sample_mapping(path):
    """Optional label remap. CSV with two columns: sample_id,label."""
    mapping = {}
    if path is None:
        return mapping
    if not os.path.exists(path):
        print(f"WARN: --sample_mapping {path} not found; ignoring", file=sys.stderr)
        return mapping
    with open(path, 'r', newline='') as fh:
        reader = csv.reader(fh)
        for row in reader:
            if not row or row[0].lower() == 'sample_id':
                continue
            if len(row) < 2:
                continue
            mapping[row[0]] = row[1]
    return mapping


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(sharing_csv_path):
    """Single pass: build full and biosample-level pair/diag counters."""
    pair_full = Counter()
    diag_full = Counter()
    pair_bio = Counter()
    diag_bio = Counter()
    all_full_seen = set()
    all_bio_seen = set()
    n_rows = 0

    for active_full, header_sample_ids in parse_sharing_csv(sharing_csv_path):
        if not all_full_seen:
            all_full_seen.update(header_sample_ids)
            all_bio_seen.update(s.rsplit('__', 1)[-1] for s in header_sample_ids)
        n_rows += 1
        for a in active_full:
            diag_full[a] += 1
        for a, b in combinations(sorted(active_full), 2):
            pair_full[(a, b)] += 1
        # Biosample view: dedup so identical biosamples across pools don't
        # create artificial self-pairs.
        active_bio = sorted(set(s.rsplit('__', 1)[-1] for s in active_full))
        for a in active_bio:
            diag_bio[a] += 1
        for a, b in combinations(active_bio, 2):
            pair_bio[(a, b)] += 1

    return {
        'n_rows': n_rows,
        'pair_full': pair_full, 'diag_full': diag_full,
        'all_full': sorted(all_full_seen),
        'pair_bio': pair_bio, 'diag_bio': diag_bio,
        'all_bio': sorted(all_bio_seen),
    }


def select_top_n(labels, diag_counts, top_n):
    """Pick the N labels with the largest diagonals; preserve a stable order."""
    if top_n is None or top_n >= len(labels):
        return list(labels)
    ranked = sorted(labels, key=lambda s: (-diag_counts.get(s, 0), s))
    return sorted(ranked[:top_n])


def remap_labels(labels, mapping):
    return [mapping.get(s, s) for s in labels]


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_pair_csv(path, labels, pair_counts, diag_counts, display_for):
    """Long-format CSV: sample_a,sample_b,shared_umis (diagonal included)."""
    label_set = set(labels)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['sample_a', 'sample_b', 'shared_umis'])
        # diagonal first (label participation count)
        for s in labels:
            n = diag_counts.get(s, 0)
            if n > 0:
                w.writerow([display_for(s), display_for(s), n])
        for (a, b), n in sorted(pair_counts.items(), key=lambda kv: -kv[1]):
            if a in label_set and b in label_set and n > 0:
                w.writerow([display_for(a), display_for(b), n])


def build_matrix(labels, pair_counts, diag_counts):
    """Return a symmetric NxN matrix M[i][j] = pair_counts.get((s_i, s_j))."""
    n = len(labels)
    idx = {s: i for i, s in enumerate(labels)}
    M = [[0] * n for _ in range(n)]
    for s, c in diag_counts.items():
        if s in idx:
            i = idx[s]
            M[i][i] = c
    for (a, b), c in pair_counts.items():
        if a in idx and b in idx:
            i, j = idx[a], idx[b]
            if i == j:
                continue
            M[i][j] = c
            M[j][i] = c
    return M


def write_heatmap_html(path, title, labels, matrix):
    """Render a Plotly heatmap to a self-contained HTML file."""
    z = matrix
    text = [
        [f"{labels[i]}<br>{labels[j]}<br>shared UMIs: {matrix[i][j]}"
         for j in range(len(labels))]
        for i in range(len(labels))
    ]
    with open(path, 'w') as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
<title>{title}</title>
<script src="{PLOTLY_CDN}"></script>
<style>
    body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 0; padding: 20px 30px; background: #f8f9fa; color: #333; }}
    h1 {{ color: #1a1a2e; border-bottom: 3px solid #4e79a7; padding-bottom: 8px; }}
    .info {{ background: #e8f0fe; padding: 10px 14px; border-radius: 6px;
             border-left: 4px solid #4e79a7; margin-bottom: 20px; }}
    .plot {{ width: 100%; max-width: 1400px; height: 800px; }}
</style>
</head>
<body>
<h1>{title}</h1>
<div class="info">
  Diagonal cells = number of UMI-sharing rows in which the sample appears at all.<br>
  Off-diagonal cells = number of rows in which both samples appear with count &gt; 0
  (= pairwise shared-UMI count).
</div>
<div id="hm" class="plot"></div>
<script>
var z = {json.dumps(z)};
var labels = {json.dumps(labels)};
var text = {json.dumps(text)};
Plotly.newPlot('hm', [{{
    z: z, x: labels, y: labels,
    text: text, hoverinfo: 'text',
    type: 'heatmap', colorscale: 'OrRd',
    colorbar: {{ title: 'shared UMIs' }}
}}], {{
    title: {json.dumps(title)},
    xaxis: {{ tickangle: -45, automargin: true }},
    yaxis: {{ automargin: true }},
    template: 'plotly_white'
}}, {{responsive: true}});
</script>
</body>
</html>
""")


# ---------------------------------------------------------------------------
# Stdout summary
# ---------------------------------------------------------------------------

def print_top_pairs(label, pair_counts, display_for, k=10):
    pairs = sorted(pair_counts.items(), key=lambda kv: -kv[1])[:k]
    print(f"\nTop {k} pairs ({label}):")
    for (a, b), c in pairs:
        print(f"  {display_for(a)}  <->  {display_for(b)}  :  {c}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def derive_default_prefix(input_path):
    base = os.path.basename(input_path)
    m = re.match(r'barcode_sharing_summary_(mm\d+_bc\d+)\.csv$', base)
    if m:
        return f"umi_sharing_{m.group(1)}"
    stem = os.path.splitext(base)[0]
    return f"umi_sharing_{stem}"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input_csv",
                   help="Path to a barcode_sharing_summary_mm*_bc*.csv emitted "
                        "by sweep_summary.py.")
    p.add_argument("--sample_mapping", default=None,
                   help="Optional 2-column CSV (sample_id,label) for label "
                        "remapping in the heatmap.")
    p.add_argument("--output_prefix", default=None,
                   help="Output filename prefix. Default: umi_sharing_mm{X}_bc{Y} "
                        "derived from input filename.")
    p.add_argument("--top_n", type=int, default=None,
                   help="Keep only the N samples with the most rows of "
                        "participation (per heatmap view). Default: no cap.")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    if not os.path.exists(args.input_csv):
        sys.exit(f"ERROR: input CSV not found: {args.input_csv}")

    prefix = args.output_prefix or derive_default_prefix(args.input_csv)
    mapping = parse_sample_mapping(args.sample_mapping)
    display_for = lambda s: mapping.get(s, s)

    print(f"Reading {args.input_csv} ...")
    agg = aggregate(args.input_csv)
    print(f"  Rows (shared UMIs):                 {agg['n_rows']:,}")
    print(f"  Distinct samples (pool x barcode):  {len(agg['all_full'])}")
    print(f"  Distinct biosamples:                {len(agg['all_bio'])}")

    # ---- pool x barcode view ----
    full_labels = select_top_n(agg['all_full'], agg['diag_full'], args.top_n)
    full_matrix = build_matrix(full_labels, agg['pair_full'], agg['diag_full'])
    full_display = remap_labels(full_labels, mapping)

    full_csv = f"{prefix}_pairwise_pool_barcode.csv"
    full_html = f"{prefix}_heatmap_pool_barcode.html"
    write_pair_csv(full_csv, full_labels, agg['pair_full'], agg['diag_full'], display_for)
    write_heatmap_html(full_html,
                       f"UMI sharing - pool x barcode  ({len(full_labels)} samples)",
                       full_display, full_matrix)
    print(f"\n  -> {full_csv}")
    print(f"  -> {full_html}")
    print_top_pairs("pool x barcode", agg['pair_full'], display_for)

    # ---- biosample view ----
    bio_labels = select_top_n(agg['all_bio'], agg['diag_bio'], args.top_n)
    bio_matrix = build_matrix(bio_labels, agg['pair_bio'], agg['diag_bio'])
    bio_display = remap_labels(bio_labels, mapping)

    bio_csv = f"{prefix}_pairwise_biosample.csv"
    bio_html = f"{prefix}_heatmap_biosample.html"
    write_pair_csv(bio_csv, bio_labels, agg['pair_bio'], agg['diag_bio'], display_for)
    write_heatmap_html(bio_html,
                       f"UMI sharing - biosample  ({len(bio_labels)} samples)",
                       bio_display, bio_matrix)
    print(f"\n  -> {bio_csv}")
    print(f"  -> {bio_html}")
    print_top_pairs("biosample", agg['pair_bio'], display_for)


if __name__ == "__main__":
    main()
