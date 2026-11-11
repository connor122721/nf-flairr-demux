#!/usr/bin/env bash
# Launch flairr-dsl2 (preprocess) for one locus with HyperQueue on SLURM.
# Copy, edit the paths in the block below, and run from a login node inside tmux/screen.
set -euo pipefail

# ---- edit these ----------------------------------------------------------
FLAIRR_DIR=/path/to/flairr_dsl2            # flairr-dsl2 checkout (has preprocess/ and annotate/)
RUN_DIR=/path/to/flairr_runs/RUN01         # where work/ and results/ go
SHEET=/path/to/metadata/RUN01_flairr_dsl2_TRA.tsv
LOCUS=TRA                                  # TRA or TRB
OUTPUT_LOCUS=TRA
REF_DIR=/path/to/germline/ogr-ref/
PARTITION=cpu384g
# --------------------------------------------------------------------------

module load hyperqueue nextflow

# Some clusters' login scripts `cd` into your work folder for every bash script
# (via BASH_ENV). That breaks Nextflow and HyperQueue; clear it for this run.
export BASH_ENV=

# HyperQueue server + SLURM worker allocation (one allocation queue only)
if ! hq server info >/dev/null 2>&1; then
  nohup hq server start > "$HOME/hq-server.log" 2>&1 &
  hq server wait
fi
if ! hq alloc list | grep -q "$PARTITION"; then
  hq alloc add slurm --time-limit 8h --max-worker-count 8 --backlog 2 --idle-timeout 5m \
    -- --partition="$PARTITION"
fi

mkdir -p "$RUN_DIR" && cd "$RUN_DIR"
[[ "$(bash -c pwd)" == "$PWD" ]] || { echo "bash still changes folder at startup; see README Troubleshooting" >&2; exit 1; }

nextflow run "$FLAIRR_DIR/preprocess/main.nf" -offline -resume \
  -profile apptainer,hyperqueue \
  -w "$RUN_DIR/work" \
  --input "$SHEET" \
  --locus "$LOCUS" --output_locus "$OUTPUT_LOCUS" --constant_region TR \
  --species Homo_sapiens \
  --germline_ref_dir "$REF_DIR" \
  --outdir "$RUN_DIR/results"
