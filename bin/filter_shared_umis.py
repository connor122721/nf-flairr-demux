import argparse
import os

from receptor_utils import simple_bio_seq as simple
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description="Filter shared UMIs from a list of UMIs.")
    parser.add_argument("barcode_sharing_summary", type=str, help="Path to barcode_sharing_summary.csv file.")
    parser.add_argument("fastq_suffix", type=str, help="Suffix of the FASTQ files containing UMIs (eg _IgM).")
    parser.add_argument("input_path", type=str, help="Path to input fastq files.")    
    parser.add_argument("output_path", type=str, help="Path to output fastq files.")
    parser.add_argument("--next_max_thresh", type=int, default=1, help="reads will be retained in the max_count sample if next_max_count is <= this value (default 1).")
    parser.add_argument("--multiple_thresh", type=int, default=10, help="reads will be retained in the max_count sample if max_count/next_max_count is above this value (default 10).")
    args = parser.parse_args()

    umi_prefix = "ATCAACGCAGAGT"
    
    filter_log = {}

    shared_recs = simple.read_csv(args.barcode_sharing_summary)
    remove_combos = {}
    sample_names = []

    for rec in shared_recs:
        seq = rec['umi_combo']
        filter_log[seq] = 0

    for h in shared_recs[0].keys():
        if 'count_' in h:
            sample_names.append(h.replace('count_', ''))

    print('Samples:', sample_names)

    for sample in sample_names:
        remove_combos[sample] = list()

    for rec in shared_recs:
        combo = rec['umi_combo']
        max_count = float(rec['max_count'])
        next_max_count = float(rec['next_max_count'])
        retain_max = False
        if next_max_count <= args.next_max_thresh and (max_count / next_max_count) > float(args.multiple_thresh):
            retain_max = True

        max_count = int(max_count)
        for sample in sample_names:
            count = int(rec['count_' + sample])
            if count > 0 and (count != max_count or not retain_max):
                remove_combos[sample].append(combo)

    for sample in sample_names:
        input_file = f"{args.input_path}/{sample}{args.fastq_suffix}.fastq"

        if not os.path.exists(input_file):
            print(f"Input file {input_file} does not exist. Skipping sample {sample}.")
            continue
    
        print(input_file)
        print(f"Removing {len(remove_combos[sample])} UMI combos")

        output_file = f"{args.output_path}/{sample}{args.fastq_suffix}.fastq"
        with open(output_file, "wt", encoding="ascii") as out_fh:
            removed = 0
            for rec in SeqIO.parse(input_file, "fastq"):
                seq = str(rec.seq)

                if not seq.startswith(umi_prefix):
                    seq = simple.reverse_complement(seq)

                if not seq.startswith(umi_prefix):
                    print(f"Read {rec.id} does not have the UMI prefix.")
                    SeqIO.write(rec, out_fh, "fastq")
                    continue

                if seq not in remove_combos[sample]:
                    SeqIO.write(rec, out_fh, "fastq")
                else:
                    removed += 1

                    if seq not in filter_log:
                        print(f"seq {seq} not found in filter log. This should not happen.")
                    else:
                        filter_log[seq] += 1

            print(f"{removed} reads removed.")

    print("Filter log:")
    print("seq\tCount")
    for seq, count in filter_log.items():
        print(f"{seq}\t{count}")

if __name__ == "__main__":
    main()