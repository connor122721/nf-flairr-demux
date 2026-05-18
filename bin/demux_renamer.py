#!/usr/bin/env python3
"""
Rename and move demux files based on barcode-to-sample mapping.

Reads a CSV mapping pool/barcode to sample_id, finds demux.RACE_TSO_5p--bcXXXX.{bam,bam.pbi}
files in pool directories, and renames them to <sample_id>.{bam,bai}, moving to demux/ directory.
"""

import argparse
import csv
import sys
from pathlib import Path


def load_mapping(csv_path):
    """Load barcode to sample_id mapping from CSV."""
    mapping = {}
    try:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pool = row['pool'].strip()
                barcode = row['barcode'].strip()
                sample_id = row['sample_id'].strip()
                key = (pool, barcode)
                mapping[key] = sample_id
        return mapping
    except FileNotFoundError:
        print(f"Error: CSV file not found: {csv_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        sys.exit(1)


def find_barcode_files(pool_dir, barcode):
    """Find fastq files for a given barcode in pool directory."""
    files = []
    pattern = f"*{barcode}*.fastq"
    
    matches = list(pool_dir.glob(pattern))
    
    if matches and len(matches) == 1:
        files.append(matches[0])
    
    return files


def get_output_name(input_file, sample_id, suffix=''):
    """Convert input filename to output filename."""
    if input_file.name.endswith('.fastq'):
        return f"{sample_id}{suffix}.fastq"
    else:
        raise ValueError(f"Unexpected file extension: {input_file.name}")


def main():
    parser = argparse.ArgumentParser(
        description='Rename demux files based on barcode-to-sample mapping.'
    )
    parser.add_argument(
        'csv_file',
        help='Path to barcodes_to_samples.csv'
    )
    parser.add_argument(
        'base_dir',
        nargs='?',
        default='.',
        help='Base directory containing bc2041, bc2042, bc2043 subdirs (default: current dir)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without making changes'
    )
    parser.add_argument(
        '--suffix',
        default='',
        help='Suffix to add to sample_id in output filenames (e.g., "_IgM" → CRC_7183_DS_IgM.fastq)'
    )
    
    args = parser.parse_args()
    
    csv_path = Path(args.csv_file)
    base_dir = Path(args.base_dir)
    
    if not csv_path.exists():
        print(f"Error: CSV file not found: {csv_path}", file=sys.stderr)
        sys.exit(1)
    
    if not base_dir.exists():
        print(f"Error: Base directory not found: {base_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Load mapping
    mapping = load_mapping(csv_path)
    
    # Create or verify demux output directory
    output_dir = base_dir / 'demux'
    if not output_dir.exists():
        if args.dry_run:
            print(f"[DRY RUN] Would create: {output_dir}")
        else:
            output_dir.mkdir(parents=True, exist_ok=True)
            print(f"Created output directory: {output_dir}")
    
    # Collect all operations
    operations = []
    
    for (pool, barcode), sample_id in sorted(mapping.items()):
        pool_dir = base_dir / pool
        
        if not pool_dir.exists():
            print(f"Warning: Pool directory not found: {pool_dir}", file=sys.stderr)
            continue
        
        files = find_barcode_files(pool_dir, barcode)
        
        if not files:
            print(f"Warning: No files found for {pool}/{barcode} ({sample_id})", file=sys.stderr)
            continue
        
        for input_file in files:
            output_name = get_output_name(input_file, sample_id, args.suffix)
            output_file = output_dir / output_name
            operations.append({
                'input': input_file,
                'output': output_file,
                'sample_id': sample_id,
                'pool': pool,
                'barcode': barcode
            })
    
    # Display operations
    if not operations:
        print("No files found to process.", file=sys.stderr)
        sys.exit(1)
    
    print(f"\n{'=' * 80}")
    print(f"Mapping and moving {len(operations)} file(s):")
    print(f"{'=' * 80}\n")
    
    for op in operations:
        marker = "[DRY RUN]" if args.dry_run else "[EXEC]"
        print(f"{marker} {op['pool']} / {op['barcode']} → {op['sample_id']}")
        print(f"        {op['input'].name}")
        print(f"        ↓")
        print(f"        {op['output'].name}")
        print()
    
    # Execute operations
    if not args.dry_run:
        print(f"\n{'=' * 80}")
        print("Executing moves...")
        print(f"{'=' * 80}\n")
        
        for i, op in enumerate(operations, 1):
            try:
                op['input'].rename(op['output'])
                print(f"[{i}/{len(operations)}] ✓ Moved: {op['input'].name} → {op['output'].name}")
            except Exception as e:
                print(f"[{i}/{len(operations)}] ✗ Error: {e}", file=sys.stderr)
                sys.exit(1)
        
        print(f"\n{'=' * 80}")
        print(f"Success! Moved {len(operations)} file(s) to {output_dir}")
        print(f"{'=' * 80}\n")
    else:
        print(f"\n{'=' * 80}")
        print("Dry run completed. Use without --dry-run to execute.")
        print(f"{'=' * 80}\n")


if __name__ == '__main__':
    main()
