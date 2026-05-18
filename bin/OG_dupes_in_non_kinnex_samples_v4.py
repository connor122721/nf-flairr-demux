# Check for duplicated sequence/UMI in segmented reads.
# Accepts both .fastq and .fastq.gz inputs; the path to each file is given
# in the sample mapping CSV (no directory scanning).
#
# Parallelization:
#   Each fastq is processed in its own worker process via stdlib
#   concurrent.futures.ProcessPoolExecutor. Workers return a local
#   {sequence: [read_ids]} dict; the main process merges results into the
#   global umi_combos as workers complete, and rewrites the two output CSVs
#   atomically after each merge (configurable via --checkpoint_every).
#
#   No third-party dependencies are added; only stdlib modules. SeqIO/Seq
#   from BioPython are already required by the original script.
#
# Output streams in real time: log lines are timestamped and force-flushed;
# output CSVs are written atomically (tmp + rename), so a concurrent
# reader in R or pandas never sees a half-written file.
#
# Note on partial snapshots:
#   Intermediate CSVs are partial views. A UMI combo only appears once it
#   has been seen in >=2 samples, so combos can appear in the CSV as more
#   files complete, and barcode_count / max_count / samples_in_combo can
#   only grow. Spot-check format/sanity mid-run; trust the final write for
#   quantitative analysis.
#
# Note on --sample in parallel mode:
#   With workers > 1, --sample is enforced PER FILE (each worker stops at N
#   reads from its own file), not globally. For a strict global limit, use
#   --workers 1.

import argparse
import csv
import gzip
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

from Bio import SeqIO
from Bio.Seq import Seq


SUMMARY_CSV = "barcode_sharing_summary_05152026.csv"
DETAILS_CSV = "sharing_details_05152026.csv"


def log(msg):
    """Timestamped print, force-flushed. Safe to call from worker processes:
    on Linux, child stdout is inherited from the parent and per-line writes
    under PIPE_BUF (4096) are atomic, so interleaving stays line-clean."""
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def open_fastq(path):
    """Open a fastq file in text mode, transparently handling gzip."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def default_worker_count():
    """Honor SLURM/cgroup CPU affinity when present, else fall back to cpu_count."""
    if hasattr(os, "sched_getaffinity"):
        try:
            return len(os.sched_getaffinity(0))
        except OSError:
            pass
    return os.cpu_count() or 1


# ---------------------------------------------------------------------------
# Worker: runs in a child process. Must be module-level for pickling.
# ---------------------------------------------------------------------------

def process_one_fastq(fastq_path, sample_name, umi_prefix, progress_every,
                      sample_limit, verbose):
    """Parse one fastq and return UMI-grouped read IDs for this sample.

    Returns a dict with keys:
        fastq_path, sample_name, umis ({seq: [read_ids]}),
        read_count, valid_count, no_index_count
    """
    pid = os.getpid()
    log(f"  [pid={pid}] start: {sample_name}  ({fastq_path})")

    local_umis = {}
    read_count = 0
    valid_count = 0
    no_index_count = 0

    with open_fastq(fastq_path) as f:
        for record in SeqIO.parse(f, "fastq"):
            if sample_limit and read_count >= sample_limit:
                break

            if read_count > 0 and read_count % progress_every == 0:
                log(f"  [pid={pid} {sample_name}] {read_count:,} reads "
                    f"({valid_count:,} valid)")

            raw_seq = str(record.seq).upper()
            read_count += 1

            if not raw_seq.startswith(umi_prefix):
                raw_seq = str(Seq(raw_seq).reverse_complement())

            if not raw_seq.startswith(umi_prefix):
                if verbose:
                    log(f"  [pid={pid} {sample_name}] {record.id}: no UMI prefix")
                no_index_count += 1
                continue

            valid_count += 1
            if raw_seq not in local_umis:
                local_umis[raw_seq] = []
            local_umis[raw_seq].append(record.id)

    log(f"  [pid={pid}] done:  {sample_name}  "
        f"({read_count:,} reads, {valid_count:,} valid)")

    return {
        "fastq_path": fastq_path,
        "sample_name": sample_name,
        "umis": local_umis,
        "read_count": read_count,
        "valid_count": valid_count,
        "no_index_count": no_index_count,
    }


# ---------------------------------------------------------------------------
# CSV writer (atomic, main process only).
# ---------------------------------------------------------------------------

def write_outputs(umi_combos, sample_order, sample_to_pool):
    """Rewrite both output CSVs atomically. Returns (n_summary, n_detail)."""
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

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
    parser.add_argument(
        "-s", "--sample",
        type=int,
        help="Limit each file (worker) to first N reads. "
             "With --workers > 1, this is per file, not global.",
    )
    parser.add_argument(
        "--progress_every",
        type=int,
        default=100000,
        help="Each worker emits a progress line every N reads (default: 100000)",
    )
    parser.add_argument(
        "--checkpoint_every",
        type=int,
        default=1,
        help="Rewrite output CSVs after every N completed files (default: 1). "
             "Set to 0 to disable intermediate writes (final write only).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=default_worker_count(),
        help=f"Number of parallel worker processes "
             f"(default: detected = {default_worker_count()}). "
             f"Use --workers 1 for fully serial behavior.",
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

    n_workers = max(1, min(args.workers, n_files))
    log(f"Worker pool size: {n_workers}")
    log(f"Checkpoint mode: "
        + (f"rewrite CSVs every {args.checkpoint_every} completed file(s)"
           if args.checkpoint_every > 0 else "final write only"))

    # Global accumulators (only the main process touches these)
    read_count = 0
    valid_read_count = 0
    no_index_count = 0
    umi_combos = {}  # seq -> {sample_name: [read_ids]}
    total_barcode_counts = {s: 0 for s in sample_order}

    # Submit every fastq to the pool; merge as workers finish.
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                process_one_fastq,
                path,
                sample,
                umi_prefix,
                args.progress_every,
                args.sample,
                args.verbose,
            ): (path, sample)
            for path, sample in file_to_sample.items()
        }
        log(f"Submitted {len(futures)} fastq job(s); processing...")

        completed = 0
        for future in as_completed(futures):
            path, sample = futures[future]
            try:
                result = future.result()
            except Exception as e:
                log(f"  [!] {sample} FAILED ({path}): {e}")
                raise

            # Merge worker output into the global umi_combos. Each sample is
            # unique to one file (enforced at CSV-load time), so the inner
            # assignment never collides with another worker's contribution.
            for seq, read_ids in result["umis"].items():
                if seq not in umi_combos:
                    umi_combos[seq] = {}
                umi_combos[seq][sample] = read_ids

            read_count += result["read_count"]
            valid_read_count += result["valid_count"]
            no_index_count += result["no_index_count"]
            total_barcode_counts[sample] = result["valid_count"]

            completed += 1
            log(f"  [{completed}/{n_files}] merged {sample}: "
                f"{result['read_count']:,} reads, "
                f"{result['valid_count']:,} valid, "
                f"{len(result['umis']):,} unique UMI combos in this file")

            if args.checkpoint_every > 0 and completed % args.checkpoint_every == 0:
                n_sum, n_det = write_outputs(umi_combos, sample_order, sample_to_pool)
                log(f"  -> checkpoint [{completed}/{n_files}]: "
                    f"{n_sum:,} summary rows, {n_det:,} detail rows written")

    log(f"Total reads in files: {read_count:,}")
    log(f"Reads without UMI prefix: {no_index_count:,}")
    log(f"Valid reads with UMI prefix: {valid_read_count:,}")
    log(f"Unique UMI combos found: {len(umi_combos):,}")

    report_stats(umi_combos, valid_read_count, total_barcode_counts)

    # Final write -- always, even if checkpoints were disabled or the last
    # completion didn't fall on a checkpoint_every boundary.
    log("Final write of output CSVs ...")
    n_sum, n_det = write_outputs(umi_combos, sample_order, sample_to_pool)
    log(f"  -> {SUMMARY_CSV}: {n_sum:,} rows")
    log(f"  -> {DETAILS_CSV}: {n_det:,} rows")
    log("Done.")


if __name__ == "__main__":
    main()
