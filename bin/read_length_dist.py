#!/usr/bin/env python
"""
read_length_dist.py - Read length distributions for assigned vs unassigned FASTQs.
Generates a summary CSV and an interactive HTML plot.
"""

import argparse
import os
import glob
import json
from collections import defaultdict

def count_read_lengths(fastq_path):
    """Parse FASTQ and return list of read lengths."""
    lengths = []
    with open(fastq_path, 'r') as f:
        line_num = 0
        for line in f:
            line_num += 1
            if line_num % 4 == 2:  # sequence line
                lengths.append(len(line.strip()))
    return lengths

def main(args):
    fastq_files = glob.glob(os.path.join(args.input_dir, '*.fastq'))
    if not fastq_files:
        print("No FASTQ files found!")
        return

    # Collect lengths per file and categorize
    file_data = {}
    all_assigned = []
    all_unassigned = []

    for fq in sorted(fastq_files):
        name = os.path.basename(fq).replace('.fastq', '')
        lengths = count_read_lengths(fq)
        if not lengths:
            continue

        file_data[name] = lengths
        if name == 'unassigned':
            all_unassigned.extend(lengths)
        else:
            all_assigned.extend(lengths)

    # Write summary CSV
    csv_path = f"{args.output_prefix}_summary.csv"
    with open(csv_path, 'w') as f:
        f.write("category,sample,count,mean_length,median_length,min_length,max_length\n")
        for name, lengths in sorted(file_data.items()):
            lengths_sorted = sorted(lengths)
            n = len(lengths_sorted)
            mean_len = sum(lengths_sorted) / n
            median_len = lengths_sorted[n // 2]
            category = "unassigned" if name == "unassigned" else "assigned"
            f.write(f"{category},{name},{n},{mean_len:.1f},{median_len},{min(lengths_sorted)},{max(lengths_sorted)}\n")

        # Totals
        if all_assigned:
            a_sorted = sorted(all_assigned)
            n = len(a_sorted)
            f.write(f"assigned,ALL_ASSIGNED,{n},{sum(a_sorted)/n:.1f},{a_sorted[n//2]},{min(a_sorted)},{max(a_sorted)}\n")
        if all_unassigned:
            u_sorted = sorted(all_unassigned)
            n = len(u_sorted)
            f.write(f"unassigned,ALL_UNASSIGNED,{n},{sum(u_sorted)/n:.1f},{u_sorted[n//2]},{min(u_sorted)},{max(u_sorted)}\n")

    # Build histogram data for HTML plot
    def make_histogram(lengths, bin_size=50):
        if not lengths:
            return [], []
        min_l = min(lengths)
        max_l = max(lengths)
        bins = list(range(0, max_l + bin_size, bin_size))
        counts = [0] * len(bins)
        for l in lengths:
            idx = l // bin_size
            if idx < len(counts):
                counts[idx] += 1
        return bins, counts

    # Per-barcode histograms + assigned/unassigned aggregates
    traces = []
    colors = [
        '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
        '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac',
        '#86bcb6', '#d37295', '#8cd17d', '#b6992d', '#499894'
    ]

    # Aggregate traces first
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

    # Per-barcode traces (hidden by default)
    color_idx = 0
    for name, lengths in sorted(file_data.items()):
        if name == 'unassigned':
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

    html_path = f"{args.output_prefix}_dist.html"
    with open(html_path, 'w') as f:
        f.write(f"""<!DOCTYPE html>
<html>
<head>
<title>Read Length Distributions</title>
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
<h1>Read Length Distributions: Assigned vs Unassigned</h1>

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
        for name, lengths in sorted(file_data.items()):
            lengths_sorted = sorted(lengths)
            n = len(lengths_sorted)
            mean_len = sum(lengths_sorted) / n
            median_len = lengths_sorted[n // 2]
            cat = "unassigned" if name == "unassigned" else "assigned"
            css_class = cat
            f.write(f'<tr><td class="{css_class}">{cat}</td><td>{name}</td>'
                    f'<td>{n:,}</td><td>{mean_len:,.1f}</td><td>{median_len:,}</td>'
                    f'<td>{min(lengths_sorted):,}</td><td>{max(lengths_sorted):,}</td></tr>\n')

        # Totals
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

    print(f"Summary CSV: {csv_path}")
    print(f"HTML report: {html_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Read length distributions for demuxed FASTQs')
    parser.add_argument('--input_dir', required=True, help='Directory containing FASTQ files')
    parser.add_argument('--output_prefix', default='read_length', help='Output file prefix')
    args = parser.parse_args()
    main(args)
