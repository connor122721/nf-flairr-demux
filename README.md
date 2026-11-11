# nf-flairr-demux

Demultiplexing and running FLAIRR-seq on Kinnex (PacBio) TCR libraries, from the
sequencing core's HiFi BAMs (or FASTQs) to per-sample, per-chain clone tables.

```
HiFi BAM (one per SMRT cell / Kinnex adaptor)      or      deconcatenated FASTQ
   │                                                            │
   ├─ 1. skera split            deconcatenate Kinnex arrays     │  (skipped for FASTQ input)
   │                                                            │
   ├─ 2. lima  (only if needed) split chains pooled under the same barcodes
   ├─ 3. simple_demux.py   ◄────────────────────────────────────┘
   │                            TSO + sample barcode + UMI check  →  one FASTQ per sample
   ├─ 4. shared-UMI filter      drop reads whose UMI is shared across samples
   │
   └─ 5. flairr-dsl2 (Nextflow + HyperQueue)   preprocess → annotate → Change-O clone tables
```

All file names, patients and barcodes below are dummy examples (`PT01`, `PT02`,
`RUN01`, ...). Example input files are in [`example/`](example/).

---

## 0. Which path do I need?

Look at what the core delivered for each SMRT cell:

| You have | What to do |
|---|---|
| Per-sample, per-chain FASTQs already, e.g. `RUN01.RACE_TSO_5p--bc1001_TRA.fastq` (the core ran lima) | Skip to **[5. Run flairr-dsl2](#5-run-flairr-dsl2-with-hyperqueue--nextflow)**. |
| One HiFi BAM per Kinnex adaptor, with **one chain per adaptor** (e.g. TRA on `bcM0001`, TRBC1 on `bcM0002`) | Run the **Nextflow pipeline** ([section 3](#3-demultiplex-samples-nextflow-pipeline)). It does skera, demux and filtering. |
| FASTQs (or already-deconcatenated BAMs) **with the TSO and sample barcodes still on the reads**, one chain per file: e.g. `samtools fastq` of a skera BAM, or core lima FASTQs that were not clipped | Run the **Nextflow pipeline** with a samplesheet ([section 3](#3-demultiplex-samples-nextflow-pipeline)). FASTQs skip skera; for BAMs pass `--kinnex false`. |
| One HiFi BAM holding **several chains with the same sample barcodes** (e.g. TRA, TRBC1 and TRBC2 all on `bc1001`) | Split the chains with **lima** first ([section 2](#2-optional-split-pooled-chains-with-lima)), then run `simple_demux.py` per chain. |

The sample demux matches reads on the sample barcode only. If two chains share a
barcode inside one BAM or FASTQ, it cannot tell them apart, which is why the lima step exists.

**FASTQ input must keep the TSO and sample barcodes.** The demux finds samples by
looking for the `RACE_TSO_5p` sequence at one end of each read and the sample barcode at
the other. Lima trims those off by default, so standard lima output (including most
core-delivered `*.RACE_TSO_5p--bcXXXX_<chain>.fastq` files) is already demultiplexed and
**cannot** be demuxed again. Feed those straight to flairr-dsl2 (section 5). To check a
FASTQ, run `simple_demux.py` on it (section 2). If almost every read in `audit_log.csv`
says `REJECTED - TSO NOT FOUND`, the adapters have been trimmed.

---

## 1. Setup

### Software
- Nextflow ≥ 23.04 (tested with 25.10) and Java 17 or 21
- Apptainer or Singularity
- SLURM (and HyperQueue for the flairr-dsl2 step)

### Build the container (once)
The container holds `skera`, `lima`, `pysam`, `biopython`, `pigz` and `receptor_utils`.
```bash
apptainer build demux_py3.sif containers/demux_py3.def
```
Pass its path with `--demux_container /path/to/demux_py3.sif` (the default in
`main_v1.nf` points to a personal folder).

To turn a skera BAM into FASTQ (e.g. to share it, or to demux it later with `--input_fastqs`):
```bash
apptainer exec demux_py3.sif skera split movie.hifi_reads.bam mas_primers.fasta movie.skera.bam
samtools fastq movie.skera.bam | gzip > movie.skera.fastq.gz     # or: bam2fastq -o movie.skera movie.skera.bam
```

### Reference files you need
| File | What it is |
|---|---|
| `mas_primers.fasta` | Kinnex MAS adapter sequences for `skera` (from PacBio, matching your Kinnex kit). Only needed for BAM input. |
| `flairr_barcodes.fasta` | One entry named **`RACE_TSO_5p`** (the 5′ RACE template-switch oligo) plus one entry per sample barcode (`bc1001`, `bc1002`, ...). See [`example/flairr_barcodes.example.fasta`](example/flairr_barcodes.example.fasta); replace the placeholder sequences with your real ones. The TSO entry must be named exactly `RACE_TSO_5p`. |
| `biosample.csv` | Sample barcode → sample name map (below) |

### Samplesheet or biosample sheet
The easiest input is a **samplesheet** with one row per sample that names its FASTQ or
BAM; see [section 3](#3-demultiplex-samples-nextflow-pipeline). Each file gets its own
barcode → sample map, so the same barcode can mean different samples in different files.

The older way is a folder of inputs plus one **biosample sheet**: header
`SampleName,BarcodeName`, one row per sample barcode
([`example/biosample.example.csv`](example/biosample.example.csv)):
```csv
SampleName,BarcodeName
PT01_D0,bc1001
PT01_D14,bc1002
PT01_D28,bc1003
PT01_TIL,bc1004
PT02_D0,bc1005
PT02_D14,bc1006
```
That map is **global**: the same barcode gets the same name in every input file. If
you demux several BAMs in one run and they reuse barcodes for different samples,
use the barcode itself as the name
([`example/biosample.shared_barcodes.example.csv`](example/biosample.shared_barcodes.example.csv)).
Then attach real sample names in the flairr-dsl2 sheet (section 4). Output files
always carry the input BAM / FASTQ name, so nothing collides.

---

## 2. (Optional) Split pooled chains with lima

Only needed when several chains share sample barcodes inside one BAM. Use the
container, so `skera` and `lima` match the pipeline's versions.

```bash
SIF=/path/to/demux_py3.sif
BAM=m00000_260101_000000_s1.hifi_reads.bam

# 1. deconcatenate the Kinnex arrays
apptainer exec $SIF skera split $BAM mas_primers.fasta RUN01.skera.bam

# 2. split by chain. chain_barcodes.fasta pairs RACE_TSO_5p with the chain-specific
#    3' (C-region) barcodes, so output pairs are named like RACE_TSO_5p--bc1001_TRA
apptainer exec $SIF lima RUN01.skera.bam chain_barcodes.fasta RUN01.lima.bam \
    --hifi-preset ASYMMETRIC --split-named --num-threads 16

# 3. check the split
cat RUN01.lima.summary
cat RUN01.lima.counts
```
Match the lima settings to how your libraries were built. If the core has already
run lima once, copy its barcode FASTA and settings (its `*.lima.summary` and
`*.lima.counts` list what it used). Many cores already deliver the per-pair FASTQs,
in which case you can skip straight to section 5.

Then run the sample demux on each chain's BAM with the script directly.
`simple_demux.py` reads BAM or FASTQ (`.fastq`, `.fq`, optionally `.gz`), picked from the
file extension:
```bash
for CHAIN in TRA TRBC1 TRBC2; do
  mkdir -p demux_$CHAIN && cd demux_$CHAIN
  apptainer exec $SIF python /path/to/nf-flairr-demux/bin/simple_demux.py \
      --input ../RUN01.lima.<pair_for_${CHAIN}>.bam \
      --barcodes ../flairr_barcodes.fasta \
      --biosample ../biosample.csv \
      --mismatches 0 \
      --max-bc-dist 1 \
      --suffix _${CHAIN} \
      --cpus 16
  cd ..
done
```
Each folder gets one `<SampleName>_<CHAIN>.fastq` per sample, plus `unassigned_<CHAIN>.fastq`,
`audit_log.csv` (one line per read with its decision) and
`read_length_summary_<CHAIN>.csv` / `.html`.

### What `simple_demux.py` checks
For each read (either orientation):
1. Finds the `RACE_TSO_5p` sequence within `--mismatches` errors, starting anywhere in the
   first `--tso-window` bases of the read (default 1 = it must start at base 0).
2. Finds the sample barcode at the other end within `--max-bc-dist` errors, within the last
   `--bc-window` bases (default 1). Ties → `AMBIGUOUS`.
3. Checks the UMI after the TSO matches `ATCAACGCAGAGTNNNNTNNNNTNNNNTCTTGGG`.
4. Keeps reads 500–2,500 bp.
5. Trims TSO and barcode from accepted reads (turn off with `--no-trim`).

Rejected reads are written to `unassigned*.fastq`, and the reason is in `audit_log.csv`
(`TSO NOT FOUND`, `BC DISTANCE TOO HIGH`, `AMBIGUOUS BC`, `UMI MISMATCH`, `READ TOO SHORT/LONG`).

---

## 3. Demultiplex samples (Nextflow pipeline)

The pipeline demuxes each input over a grid of TSO-mismatch × barcode-distance
settings so you can pick the best combination. It also runs FastQC/MultiQC, a sweep
summary, and the shared-UMI filter.

### Input: a samplesheet (recommended)
One row per sample, CSV or TSV, with a header. Columns are found by name; any other
columns are ignored, so a full master metadata sheet works as it is.

| Needed | Accepted column names |
|---|---|
| sample name | `SampleName`, `sample_name`, `sample` |
| sample barcode | `BarcodeName`, `primer_barcode`, `barcode` |
| reads file (FASTQ or BAM) | `fastq`, `fastq_path`, `reads`, `input`, `bam`, `bam_path` (first one found wins) |
| *optional* chain | `locus`, e.g. `TRA`, `TRBC1`. If the barcode FASTA names barcodes per chain (`bc1004_TRA`), the pipeline uses `<barcode>_<locus>`; otherwise plain `<barcode>`. |
| *optional* keep / skip | `sequenced` or `include`. Rows with `no`, `n`, `false` or `0` are skipped. |

Minimal sheet ([`example/samplesheet.example.csv`](example/samplesheet.example.csv)):
```csv
SampleName,BarcodeName,locus,fastq
PT01_D0_TRA,bc1001,TRA,/path/to/core_delivery/hifi_reads/RUN01.RACE_TSO_5p--bc1001_TRA.fastq
PT01_D0_TRBC1,bc1001,TRBC1,/path/to/core_delivery/hifi_reads/RUN01.RACE_TSO_5p--bc1001_TRBC1.fastq
PT02_D0_TRBC2,bc1001,TRBC2,/path/to/RUN02/m00000_260101_000000_s2.hifi_reads.bam
PT02_D14_TRBC2,bc1002,TRBC2,/path/to/RUN02/m00000_260101_000000_s2.hifi_reads.bam
```
A master metadata TSV works too ([`example/master_metadata.example.tsv`](example/master_metadata.example.tsv),
using `sample_name`, `primer_barcode`, `locus`, `sequenced`, `fastq_path`).

- **Rows that point to the same file are demuxed together.** In the example, PT02's two samples share one BAM.
- **The barcode map is per file.** `bc1001` can be PT01 in one file and PT02 in another.
- **Rows with the same barcode in the same file must be the same sample.** Otherwise the pipeline stops.

```bash
module load nextflow
export BASH_ENV=      # see Troubleshooting
nextflow -bg run main_v1.nf -profile slurm,apptainer -resume \
    --samplesheet /path/to/metadata/RUN01_samplesheet.csv \
    --kinnex false \
    --barcodes /path/to/barcodes/FLAIRR_bc_loci_barcodes.fasta \
    --tso_window 40 --bc_window 40 \
    --bc_trim_len 16 \
    --mismatches_range 0,2 \
    --max_bc_dist_range 0,1,2 \
    --outdir ./RUN01_FLAIRR
```

- **`--kinnex`:** whether BAM rows are raw Kinnex arrays.
  - `true` (default): skera runs first, and `--mas_primers` is needed.
  - `false`: BAMs are demuxed as they are.
  - FASTQs never go through skera.
- **`--bc_trim_len 16`:** set this when the barcode FASTA has chain-specific entries that contain the C-region primer (`bc1004_TRA`). Otherwise the demux trims the primer along with the barcode, and flairr-dsl2 fails at primer masking.
- **`--tso_window` / `--bc_window`:** how far into each end of the read to look for the TSO and the barcode. The default of 1 means "must start at the first base". Raise it, e.g. to 40, if the reads carry extra bases before the TSO or after the barcode. Without that, most reads come out `TSO NOT FOUND`.
- **`-bg`:** runs Nextflow in the background, so follow it with `tail -f .nextflow.log`. Use absolute paths in the samplesheet.

### Input: folders of BAMs or FASTQs + biosample sheet
For many files that all use the same barcode → sample map:

| Input | Flags | What happens |
|---|---|---|
| Raw Kinnex HiFi BAMs (`*.hifi_reads.bam`) | `--input_bams` + `--mas_primers` + `--biosample` | skera → demux → filter |
| BAMs that don't need skera | `--input_bams` + `--kinnex false` + `--biosample` | demux → filter |
| FASTQs (`.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`) | `--input_fastqs` + `--biosample` | demux → filter |

Each flag takes one file, a glob in quotes (`'/path/*.fastq.gz'`), or a folder. A folder
means every matching file under it, at any depth. Don't send BAMs that have already been
through skera or lima back through skera (use `--kinnex false`).

```bash
nextflow run main_v1.nf -profile slurm,apptainer -resume \
    --input_bams /path/to/RUN01/hifi_bams/ \
    --mas_primers /path/to/mas_primers.fasta \
    --barcodes    /path/to/flairr_barcodes.fasta \
    --biosample   /path/to/biosample.csv \
    --mismatches_range 0 \
    --max_bc_dist_range 0,1,2,3 \
    --partition compute \
    --outdir /path/to/RUN01/demux_results
```

### Input names
Each input file gets an ID that appears in every output file name:
- BAM through skera: `<bam name>.skera`
- BAM with `--kinnex false`: `<bam name>`
- FASTQ: the file name minus `.fastq(.gz)` / `.fq(.gz)`, e.g. `RUN01.RACE_TSO_5p--bc1004_TRA.fastq` → `RUN01.RACE_TSO_5p--bc1004_TRA`

IDs must be unique and must not contain a double underscore `__`; the pipeline stops
with an error otherwise. Rename or symlink files that break either rule.

### Parameters
| Parameter | Default | Meaning |
|---|---|---|
| `--samplesheet` | – | per-sample CSV/TSV (use this **or** `--input_bams` / `--input_fastqs`) |
| `--input_bams` | – | HiFi BAM(s) or folder (+ `--biosample`) |
| `--input_fastqs` | – | FASTQ(s) or folder (+ `--biosample`) |
| `--biosample` | – | `SampleName,BarcodeName` CSV (not used with `--samplesheet`) |
| `--kinnex` | `true` | run skera on BAM inputs first |
| `--mas_primers` | – | Kinnex MAS adapters for skera (BAMs with `--kinnex true`) |
| `--barcodes` | – | FASTA with `RACE_TSO_5p` + sample barcodes |
| `--tso_window`, `--bc_window` | `1`, `1` | bases from each read end searched for the TSO / barcode |
| `--bc_trim_len` | `0` | trim only this many barcode bases from the read end (`0` = the whole barcode entry). Set it to the sample-barcode length (e.g. `16`) when barcode entries also contain the C-region primer (`bc1004_TRA`-style), so flairr-dsl2 can still find the primer. |
| `--mismatches_range` | `0` | TSO mismatches to try (comma list) |
| `--max_bc_dist_range` | `0,1,2,3,5` | sample-barcode errors to try (comma list) |
| `--demux_container` | personal path | demux container (`demux_py3.sif`) |
| `--outdir` | `results` | output folder |
| `--no_trim` | `false` | keep TSO / barcode on reads |
| `--run_filter_shared_umis` | `true` | run the shared-UMI filter |
| `--filter_next_max_thresh` | `1` | keep a shared UMI in its top sample if the runner-up has ≤ this many reads |
| `--filter_multiple_thresh` | `10` | ... or if the top sample has > this many times the runner-up's reads |
| `--cpus`, `--chunk_size` | `12`, `5000` | demux threads, reads per chunk |
| `--partition`, `--threads`, `--memory` | `compute`, `12`, `50 GB` | SLURM settings (`slurm` profile) |

### Outputs
```
demux_results/
├── skera/                                    deconcatenated BAMs (Kinnex BAM input only)
├── demux/mm<X>_bc<Y>/<id>/                   <SampleName>_<id>_mm<X>_bc<Y>.fastq.gz  + audit / read-length files
├── fastqc/  multiqc/mm<X>_bc<Y>/             read QC per setting (sample FASTQs; unassigned_* skipped)
├── sweep_summary/sweep_summary.html          assigned reads per sample for every setting  ← start here
├── filtered_demux/mm<X>_bc<Y>/output/        <id>__<SampleName>_mm<X>_bc<Y>.fastq.gz  (shared UMIs removed)
└── reports/                                  Nextflow execution report + timeline
```

### Choosing the settings
Open `sweep_summary/sweep_summary.html`:
- Pick the lowest barcode distance after which assigned reads stop climbing. Higher
  distances rescue a few reads but also mis-assign more.
- Check the UMI-sharing table: many UMIs shared between samples means barcode bleed.
  The `filtered_demux/` FASTQs have those removed.
- Use the `filtered_demux/mm<X>_bc<Y>/output/` FASTQs for that setting in section 4.

---

## 4. Build the flairr-dsl2 input sheets

One tab-separated file per chain group, **no header**, two columns: sample name, FASTQ
path. TRA in one file; TRBC1 and TRBC2 together in the TRB file. Use a consistent
name pattern such as `<patient>_<timepoint>_<chain>`, so the downstream analysis
can parse patient, timepoint and chain from it.

[`example/flairr_dsl2_TRA.example.tsv`](example/flairr_dsl2_TRA.example.tsv) (from this pipeline):
```
PT01_D0_TRA	/path/to/demux/results/filtered_demux/mm0_bc1/output/m00000_260101_000000_s1.hifi_reads.bcM0001.skera__PT01_D0_mm0_bc1.fastq.gz
PT01_D14_TRA	/path/to/demux/results/filtered_demux/mm0_bc1/output/m00000_260101_000000_s1.hifi_reads.bcM0001.skera__PT01_D14_mm0_bc1.fastq.gz
```
[`example/flairr_dsl2_TRB.example.tsv`](example/flairr_dsl2_TRB.example.tsv) (straight from core-delivered lima FASTQs):
```
PT01_D0_TRBC1	/path/to/core_delivery/hifi_reads/RUN01.RACE_TSO_5p--bc1001_TRBC1.fastq
PT01_D0_TRBC2	/path/to/core_delivery/hifi_reads/RUN01.RACE_TSO_5p--bc1001_TRBC2.fastq
```
Tips:
- **Absolute paths only.**
- **Unix line endings.** Sheets edited on Windows get `\r` at line ends; fix with `sed -i 's/\r$//' sheet.tsv`.
- **FASTQ-input runs** are named the same way, with the FASTQ ID in place of the BAM ID: `filtered_demux/mm0_bc1/output/RUN01_TRA__PT01_D0_mm0_bc1.fastq.gz`.
- **Drop empty barcodes.** lima and the demux write a FASTQ for any barcode with reads, including stray barcodes that were never in the pool and samples that failed library prep (a few KB). Only list the samples you sequenced.
- **Keep a master sheet** per run (sample, patient, timepoint, chain, barcode, adaptor, FASTQ path, notes) next to these files.

---

## 5. Run flairr-dsl2 with HyperQueue + Nextflow

flairr-dsl2 runs many small jobs. HyperQueue packs them onto a few SLURM worker jobs
instead of submitting thousands of tiny SLURM jobs.

[`example/run_flairr_hq.example.sh`](example/run_flairr_hq.example.sh) does all of this.
Edit the paths at the top and run it inside `tmux` or `screen` on a login node.

### 5.1 Start HyperQueue (once per session)
```bash
module load hyperqueue nextflow
export BASH_ENV=          # see Troubleshooting: stops login scripts cd-ing every job

nohup hq server start > ~/hq-server.log 2>&1 &   # survives logout; or run inside tmux
hq server wait

# one SLURM allocation queue: workers are started on demand, stop after 5 min idle
hq alloc list | grep -q cpu384g || \
hq alloc add slurm --time-limit 8h --max-worker-count 8 --backlog 2 --idle-timeout 5m \
    -- --partition=cpu384g
```
Check it: `hq server info`, `hq worker list`, `hq alloc list`, `squeue -u $USER`.
Running `hq alloc add` again adds another queue; check `hq alloc list` first.

### 5.2 Preprocess, then annotate, per chain
```bash
cd /path/to/flairr_runs/RUN01          # work/ and .nextflow.log go here

nextflow run /path/to/flairr_dsl2/preprocess/main.nf -offline -resume \
  -profile apptainer,hyperqueue \
  -w /path/to/flairr_runs/RUN01/work \
  --input /path/to/metadata/RUN01_flairr_dsl2_TRA.tsv \
  --locus TRA --output_locus TRA --constant_region TR \
  --species Homo_sapiens \
  --germline_ref_dir /path/to/germline/ogr-ref/ \
  --outdir /path/to/flairr_runs/RUN01/results
```
- **TRB:** use the TRB sheet with `--locus TRB` and the TRB output settings from the flairr-dsl2 README.
- **Annotate:** when preprocess has finished, run `annotate/main.nf` with the same sheet and locus to get the Change-O / AIRR clone tables (`*clone_pass.tsv`).
- **Paths:** keep every path absolute (`main.nf`, `-w`, `--input`, `--outdir`).
- **Monitor:** watch with `tail -f .nextflow.log` and `hq job list`.

### 5.3 Results layout for downstream analysis
Put the clone tables in one folder per chain, one sub-folder per sample:
```
flairr_results/
├── TRA/PT01_D0_TRA/.../clones/*clone_pass.tsv
├── TRBC1/PT01_D0_TRBC1/.../clones/*clone_pass.tsv
└── TRBC2/PT01_D0_TRBC2/.../clones/*clone_pass.tsv
```

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Cannot find script file: main.nf`, or `Unknown project nextflow-io/...`, while `main.nf` is right there | A system login script (e.g. `/etc/profile.d/*.sh`, loaded for every bash script via `BASH_ENV`) runs `cd $WORK`. The Nextflow launcher is a bash script, so it starts in the wrong folder. Test: `bash -c pwd` prints a different folder than `pwd`. | `export BASH_ENV=` before starting `hq` and `nextflow`. Ask the admins to guard the line with `[[ $- == *i* ]]`. |
| `work/` appears a few folders up | Same as above | Same, and pass `-w /absolute/path/work` |
| HyperQueue tasks fail; `.command.log/*.hqs` says `bash: .command.run: No such file or directory` | Same as above: the task's bash `cd`s away from its work folder | Same. Restart the HQ server **after** `export BASH_ENV=` so workers inherit it. |
| Run "hangs", nothing in `hq job list` | HQ server died (started with plain `&` and you logged out), or no workers | `hq server info`, `hq worker list`, `hq alloc list`, `squeue -u $USER`; restart with `nohup` / tmux |
| Jobs stay `WAITING` with workers running | Task asks for more CPUs / memory than a worker has | Lower the process resources, or use a bigger partition in `hq alloc add` |
| `.hqs` log files look like binary | HQ stream format | `hq job cat <id> stderr` (or `stdout`) |
| Few reads assigned in the demux | Wrong `RACE_TSO_5p` / barcode sequences, or settings too strict | Check `audit_log.csv` decisions and `sweep_summary.html`; try a higher `--max_bc_dist_range` |
| FASTQ input: nearly all reads `REJECTED - TSO NOT FOUND` | The FASTQ was trimmed by lima (already demultiplexed) | Don't re-demux it; use it directly in flairr-dsl2 |
| `input name '...' contains '__'` / `input files share a name` | An input file name has a double underscore, or two inputs have the same name | Rename or symlink the files |
| `give one of --samplesheet, --input_bams or --input_fastqs` / `give only one of ...` | No input, or more than one kind | Pass exactly one |
| `--samplesheet needs columns for ...` | Header names not recognised | Rename columns to one of the accepted names (section 3) |
| `barcode ... not found in <barcodes>` | Samplesheet barcode (or `<barcode>_<locus>`) isn't a FASTA header | Match names to the barcode FASTA, or add / fix the `locus` column |
| Most reads `TSO NOT FOUND` though the TSO is in the reads | TSO doesn't start at base 0 | Raise `--tso_window` (and `--bc_window`) |
| flairr-dsl2 MaskPrimers: almost every read `primers-fail` (C primer not found) | Barcode entries include the C-region primer (`bc1004_TRA` = barcode + primer) and the demux trimmed the whole entry | Re-run the demux with `--bc_trim_len 16` (your sample-barcode length) |
| Same barcode used by different samples in one run | Biosample map is global | Use barcode names in `biosample.csv` and rename in the flairr-dsl2 sheet |
| Chains mixed in one sample's FASTQ | Chains pooled under the same barcodes without lima | Split chains with lima first (section 2) |
| Bad paths / trailing `\r` in sample names | Sheet edited on Windows | `sed -i 's/\r$//' sheet.tsv` |
