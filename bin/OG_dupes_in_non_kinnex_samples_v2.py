# Check for duplicated sequence/UMI in segmented reads.
# Accepts both .fastq and .fastq.gz inputs; the path to each file is given
# in the sample mapping CSV (no directory scanning).
#
# Output streams in real time:
#   - Log lines are timestamped and force-flushed.
#   - barcode_sharing_summary.csv and sharing_details.csv are rewritten
#     atomically after each fastq is processed (configurable via
#     --checkpoint_every), so you can peek at intermediate results in R /
#     pandas while the script is still running. Atomic replace (tmp + rename)
#     means readers never see a half-written file.
#
# Important: intermediate CSVs are partial views. A UMI combo only appears
# once it's been seen in >=2 samples, so combos can appear as later files
# are processed, and barcode_count / max_count / samples_in_combo can only
# grow with each checkpoint. Use mid-run snapshots for format/sanity checks;
# trust the final write for quantitative analysis.

import argparse
import csv
import gzip
import os
import sys
from datetime import datetime

from Bio import SeqIO
from Bio.Seq import Seq

SUMMARY_CSV = "barcode_sharing_summary_05_13_2026.csv"
DETAILS_CSV = "sharing_details_05_13_2026.csv"

def log(msg):
    """Timestamped print, force-flushed so it appears in the log immediately."""
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def open_fastq(path):
    """Open a fastq file in text mode, transparently handling gzip."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def write_outputs(umi_combos, sample_order, sample_to_pool):
    """Rewrite both output CSVs atomically. Returns (n_summary_rows, n_detail_rows).

    Each CSV is written to a `.tmp` sibling, then `os.replace`d into place so
    a concurrent reader either sees the previous complete file or the new
    complete file -- never a half-written one.
    """
    # ---- summary ----
    n_summary_rows = 0
    summary_tmp = SUMMARY_CSV + ".tmp"
    with open(summary_tmp, "w", newline="") as csvfile:
        fieldnames = (
            ["umi_combo", "barcode_count"]
            + [f"count_{bc}" for bc in sample_order]
            + ["total_count", "max_count", "next_max_count", "pools"]
        )
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for umi_combo, barcodes in umi_combos.items():
            if len(barcodes) <= 1:
                continue
            rec = {
                "umi_combo": umi_combo,
                "barcode_count": len(barcodes),
                **{f"count_{bc}": 0 for bc in sample_order},
            }
            for barcode_name, recs in barcodes.items():
                rec[f"count_{barcode_name}"] = len(recs)
            rec["total_count"] = sum(rec[f"count_{bc}"] for bc in sample_order)
            rec["max_count"] = max(rec[f"count_{bc}"] for bc in sample_order)
            rec["next_max_count"] = (
                sorted([rec[f"count_{bc}"] for bc in sample_order], reverse=True)[1]
                if len(sample_order) > 1
                else 0
            )
            rec["pools"] = "|".join(
                sorted(
                    set(
                        sample_to_pool[bc]
                        for bc in sample_order
                        if rec[f"count_{bc}"] > 0
                    )
                )
            )
            writer.writerow(rec)
            n_summary_rows += 1
    os.replace(summary_tmp, SUMMARY_CSV)

    # ---- details ----
    n_detail_rows = 0
    details_tmp = DETAILS_CSV + ".tmp"
    with open(details_tmp, "w", newline="") as csvfile:
        fieldnames = ["umi", "sample", "samples_in_combo", "reads_in_combo", "zmw", "umi_combo"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for umi_combo, barcodes in umi_combos.items():
            if len(barcodes) <= 1:
                continue
            for barcode_name, recs in barcodes.items():
                for rec in recs:
                    writer.writerow({
                        "umi": umi_combo[:30],
                        "sample": barcode_name,
                        "samples_in_combo": len(barcodes.keys()),
                        "reads_in_combo": len(recs),
                        "zmw": rec,
                        "umi_combo": umi_combo,
                    })
                    n_detail_rows += 1
    os.replace(details_tmp, DETAILS_CSV)

    return n_summary_rows, n_detail_rows


def main():
    # Belt-and-suspenders: even without `python -u`, force line-buffered
    # stdout/stderr so anything that escapes log() still flushes promptly.
    sys.stdout.reconfigure(line_buffering=True)
    sys.stderr.reconfigure(line_buffering=True)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample_mapping",
        required=True,
        help="CSV with columns: sample_id, pool, fastq_path "
             "(fastq_path may be .fastq or .fastq.gz)",
    )
    parser.add_argument("-s", "--sample", type=int, help="Limit to first n reads (global)")
    parser.add_argument(
        "--progress_every",
        type=int,
        default=100000,
        help="Emit a progress line every N reads (default: 100000)",
    )
    parser.add_argument(
        "--checkpoint_every",
        type=int,
        default=1,
        help="Rewrite output CSVs after every N processed files (default: 1, "
             "i.e. after every file). Set to 0 to disable intermediate writes "
             "(only write at the end).",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    umi_prefix = "ATCAACGCAGAGT"

    # Read sample mapping. Each row = one fastq.
    # pool is treated as an opaque run/batch label (any string).
    sample_order = []
    sample_to_pool = {}
    file_to_sample = {}  # absolute fastq path -> sample_id

    log(f"Reading sample mapping: {args.sample_mapping}")
    with open(args.sample_mapping, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        required_cols = {"sample_id", "pool", "fastq_path"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"sample_mapping is missing required columns: {sorted(missing)}. "
                f"Found: {reader.fieldnames}"
            )

        for row in reader:
            sample = row["sample_id"].replace(".fastq.gz", "").replace(".fastq", "")
            fastq_path = os.path.abspath(os.path.expanduser(row["fastq_path"].strip()))

            if not os.path.isfile(fastq_path):
                raise FileNotFoundError(f"fastq file not found: {fastq_path}")
            if sample in sample_to_pool:
                raise ValueError(f"Duplicate sample_id in sample_mapping: {sample}")

            sample_order.append(sample)
            sample_to_pool[sample] = row["pool"].strip()
            file_to_sample[fastq_path] = sample

    n_files = len(file_to_sample)
    n_pools = len(set(sample_to_pool.values()))
    log(f"Loaded {n_files} fastq entries across {n_pools} pool(s)")
    log(f"Checkpoint mode: "
        + (f"rewrite CSVs every {args.checkpoint_every} file(s)"
           if args.checkpoint_every > 0 else "final write only"))

    read_count = 0
    valid_read_count = 0
    no_index_count = 0

    umi_combos = {}
    total_barcode_counts = {s: 0 for s in sample_order}

    # Process each fastq / fastq.gz
    for file_idx, (fastq_path, sample_name) in enumerate(file_to_sample.items(), start=1):
        log(f"[{file_idx}/{n_files}] Processing: {fastq_path} (sample={sample_name})")
        file_reads_before = read_count
        file_valid_before = valid_read_count

        with open_fastq(fastq_path) as f:
            for record in SeqIO.parse(f, "fastq"):
                if args.sample and read_count >= args.sample:
                    break

                if read_count > 0 and read_count % args.progress_every == 0:
                    log(f"  ...processed {read_count:,} reads total ({valid_read_count:,} valid)")
                    report_stats(umi_combos, valid_read_count, total_barcode_counts)

                raw_seq = str(record.seq).upper()
                read_count += 1

                if not raw_seq.startswith(umi_prefix):
                    raw_seq = str(Seq(raw_seq).reverse_complement())

                if not raw_seq.startswith(umi_prefix):
                    if args.verbose:
                        log(f"Read {record.id} does not start with UMI prefix, skipping")
                    no_index_count += 1
                    continue

                valid_read_count += 1
                total_barcode_counts[sample_name] += 1
                umi_combo = raw_seq

                if umi_combo not in umi_combos:
                    umi_combos[umi_combo] = {}
                if sample_name not in umi_combos[umi_combo]:
                    umi_combos[umi_combo][sample_name] = []
                umi_combos[umi_combo][sample_name].append(record.id)

        file_reads = read_count - file_reads_before
        file_valid = valid_read_count - file_valid_before
        log(f"  -> finished {sample_name}: {file_reads:,} reads ({file_valid:,} valid)")

        # Incremental CSV checkpoint
        if args.checkpoint_every > 0 and file_idx % args.checkpoint_every == 0:
            n_sum, n_det = write_outputs(umi_combos, sample_order, sample_to_pool)
            log(f"  -> checkpoint [{file_idx}/{n_files}]: "
                f"{n_sum:,} summary rows, {n_det:,} detail rows written")

        if args.sample and read_count >= args.sample:
            log(f"Reached --sample limit ({args.sample:,}); stopping early.")
            break

    log(f"Total reads in files: {read_count:,}")
    log(f"Reads without UMI prefix: {no_index_count:,}")
    log(f"Valid reads with UMI prefix: {valid_read_count:,}")
    log(f"Unique UMI combos found: {len(umi_combos):,}")

    report_stats(umi_combos, valid_read_count, total_barcode_counts)

    # Final write -- always do it, even if checkpoints were disabled or the
    # last file didn't land on a checkpoint boundary.
    log("Final write of output CSVs ...")
    n_sum, n_det = write_outputs(umi_combos, sample_order, sample_to_pool)
    log(f"  -> {SUMMARY_CSV}: {n_sum:,} rows")
    log(f"  -> {DETAILS_CSV}: {n_det:,} rows")
    log("Done.")


def report_stats(umi_combos, total_reads, total_barcode_counts):
    dupes = 0
    max_shared = 0
    duplicated_reads = 0

    for umi_combo, barcodes in umi_combos.items():
        if len(barcodes) > 1:
            dupes += 1
            if len(barcodes) > max_shared:
                max_shared = len(barcodes)
            duplicated_reads += sum(len(recs) for recs in barcodes.values())

    perc_shared = (
        round(dupes / len(umi_combos) * 100, 2) if len(umi_combos) > 0 else 0.00
    )
    perc_duplicated_reads = (
        round(duplicated_reads / total_reads * 100, 2) if total_reads > 0 else 0.00
    )
    log(f"  Shared UMI combos: {dupes:,} ({perc_shared}% of unique combos)")
    log(f"  Reads in shared combos: {duplicated_reads:,} ({perc_duplicated_reads}% of valid reads)")
    log(f"  Max barcodes sharing a single combo: {max_shared}")
    log(
        "  Per-barcode totals: "
        + ", ".join(f"{bc}={n:,}" for bc, n in total_barcode_counts.items())
    )


if __name__ == "__main__":
    main()