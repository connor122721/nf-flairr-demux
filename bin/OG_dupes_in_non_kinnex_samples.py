# Check for duplicated sequence/UMI in segmented reads.
# Accepts both .fastq and .fastq.gz inputs; the path to each file is given
# in the sample mapping CSV (no directory scanning).

import argparse
import csv
import gzip
import os

from Bio import SeqIO
from Bio.Seq import Seq


def open_fastq(path):
    """Open a fastq file in text mode, transparently handling gzip."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample_mapping",
        required=True,
        help="CSV with columns: sample_id, pool, fastq_path "
             "(fastq_path may be .fastq or .fastq.gz)",
    )
    parser.add_argument("-s", "--sample", type=int, help="Limit to first n reads (global)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")
    args = parser.parse_args()

    umi_prefix = "ATCAACGCAGAGT"

    # Read sample mapping. Each row = one fastq.
    # pool is treated as an opaque run/batch label (any string).
    sample_order = []
    sample_to_pool = {}
    file_to_sample = {}  # absolute fastq path -> sample_id

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

    read_count = 0
    valid_read_count = 0
    no_index_count = 0

    umi_combos = {}
    total_barcode_counts = {s: 0 for s in sample_order}

    # Process each fastq / fastq.gz
    for fastq_path, sample_name in file_to_sample.items():
        print(f"Processing file: {fastq_path}  (sample={sample_name})")
        with open_fastq(fastq_path) as f:
            for record in SeqIO.parse(f, "fastq"):
                if args.sample and read_count >= args.sample:
                    break

                if read_count > 0 and read_count % 100000 == 0:
                    print(f"Processed {read_count} reads...")
                    report_stats(umi_combos, valid_read_count, total_barcode_counts)

                raw_seq = str(record.seq).upper()
                read_count += 1

                if not raw_seq.startswith(umi_prefix):
                    raw_seq = str(Seq(raw_seq).reverse_complement())

                if not raw_seq.startswith(umi_prefix):
                    if args.verbose:
                        print(f"Read {record.id} does not start with UMI prefix, skipping")
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

        if args.sample and read_count >= args.sample:
            break  # also exit the outer loop, not just the inner one

    print(f"Total reads in files: {read_count}")
    print(f"Reads without UMI prefix: {no_index_count}")
    print(f"Valid reads with UMI prefix: {valid_read_count}")
    print(f"Unique UMI combos found: {len(umi_combos)}")

    report_stats(umi_combos, valid_read_count, total_barcode_counts)

    with open("barcode_sharing_summary.csv", "w", newline="") as csvfile:
        fieldnames = (
            ["umi_combo", "barcode_count"]
            + [f"count_{bc}" for bc in sample_order]
            + ["total_count", "max_count", "next_max_count", "pools"]
        )
        print(fieldnames)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for umi_combo, barcodes in umi_combos.items():
            if len(barcodes) > 1:
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

    with open("sharing_details.csv", "w", newline="") as csvfile:
        fieldnames = ["umi", "sample", "samples_in_combo", "reads_in_combo", "zmw", "umi_combo"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for umi_combo, barcodes in umi_combos.items():
            if len(barcodes) > 1:
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
        str(round(dupes / len(umi_combos) * 100, 2)) if len(umi_combos) > 0 else "0.00"
    )
    perc_duplicated_reads = (
        str(round(duplicated_reads / total_reads * 100, 2)) if total_reads > 0 else "0.00"
    )
    print(f"UMI combos shared across multiple barcodes: {dupes}  ({perc_shared}% of unique combos)")
    print(f"Reads within duplicated UMI combos: {duplicated_reads} ({perc_duplicated_reads}% of valid reads)")
    print(f"Maximum number of barcodes sharing a UMI combo: {max_shared}")
    print(
        "Total count for each barcode: "
        + ", ".join(f"{bc}: {n}" for bc, n in total_barcode_counts.items())
    )
    print("\n")


if __name__ == "__main__":
    main()