#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// nf-flairr-demux: FLAIRR-seq parameter sweep (Kinnex or standard HiFi)
// ============================================================================
// Input, either:
//   --samplesheet sheet.csv   CSV with columns bam,sample_name,sample_barcode
//                             (sample_barcode: bc1001 or RACE_TSO_5p--bc1001)
//                             (extra columns are ignored). Each BAM gets its own
//                             barcode->sample map, so barcodes can be reused
//                             across BAMs.
// or (legacy):
//   --input_bams <dir|bam> --biosample <SampleName,BarcodeName csv>
//                             one barcode->sample map shared by every BAM.
//
// --kinnex true  (default) run skera split (needs --mas_primers) before demux
// --kinnex false            standard (non-Kinnex) HiFi reads: demux directly
params.samplesheet = null
params.input_bams = null
params.mas_primers = null
params.barcodes = null
params.biosample = null
params.kinnex = true
params.outdir = 'results'
params.no_trim = false
params.cpus = 12
params.chunk_size = 5000

// Parameter sweep ranges
params.mismatches_range  = '0'
params.max_bc_dist_range = '0,1,2,3,5'

// Consolidated container: built from containers/demux_py.def via
//   apptainer build demux_py.sif containers/demux_py.def
// Override on the command line: --demux_container /path/to/demux_py.sif
// Includes pbskera, lima, pysam, biopython, pigz, receptor_utils.
//params.demux_container = "${projectDir}/containers/demux_py.sif"
params.demux_container = "/home/c0murr09/sif/demux_py3.sif"

// filter_shared_umis.py thresholds (see bin/filter_shared_umis.py)
params.run_filter_shared_umis = true
params.filter_next_max_thresh = 1
params.filter_multiple_thresh = 10

// ============================================================================
// Validate
// ============================================================================
def is_kinnex = params.kinnex.toString().toLowerCase() in ['true', 'yes', '1']

if (params.samplesheet && (params.input_bams || params.biosample)) {
    error "ERROR: use either --samplesheet OR --input_bams/--biosample, not both"
}
if (!params.samplesheet && !(params.input_bams && params.biosample)) {
    error "ERROR: provide --samplesheet, or both --input_bams and --biosample"
}
if (is_kinnex && !params.mas_primers) { error "ERROR: --mas_primers is required when --kinnex true" }
if (!params.barcodes) { error "ERROR: --barcodes is required" }

// BAM id used in output names: file name minus .bam (e.g. movie_s1.hifi_reads)
def bamId(f) { f.name.replaceAll(/\.bam$/, '') }

log.info "nf-flairr-demux | kinnex: ${is_kinnex} | input: ${params.samplesheet ?: params.input_bams}"

// ============================================================================
// Processes
// ============================================================================
process SKERA_SPLIT {
    tag "${id}"
    label 'process_medium'
    container "${params.demux_container}"

    publishDir "${params.outdir}/skera", mode: 'copy', pattern: '*.skera.bam*'

    input:
    tuple val(id), path(bam), path(biosample)
    path mas_primers

    output:
    tuple val(id), path("${id}.skera.bam"), path(biosample), emit: skera_bam
    path "${id}.skera.bam.pbi", emit: skera_pbi, optional: true

    script:
    """
    skera split \\
        ${bam} \\
        ${mas_primers} \\
        ${id}.skera.bam
    """
}

process SIMPLE_DEMUX {
    tag "${id}_mm${mismatches}_bc${max_bc_dist}"
    label 'process_medium'
    container "${params.demux_container}"

    publishDir "${params.outdir}/demux/mm${mismatches}_bc${max_bc_dist}/${id}", mode: 'copy'

    input:
    tuple val(id), path(reads_bam), path(biosample), val(mismatches), val(max_bc_dist)
    path barcodes

    output:
    tuple val(mismatches), val(max_bc_dist), path("*.fastq.gz"),            emit: fastqs
    tuple val(mismatches), val(max_bc_dist), path("audit_*.csv"),           emit: audit_log
    tuple val(mismatches), val(max_bc_dist), path("read_length_*.csv"),     emit: rl_csv
    tuple val(mismatches), val(max_bc_dist), path("read_length_*.html"),    emit: rl_html
    path "demux.out"

    script:
    def trim_flag = params.no_trim ? '--no-trim' : ''
    def unique_suffix = "_${id}_mm${mismatches}_bc${max_bc_dist}"
    """
    python ${projectDir}/bin/simple_demux.py \\
        --input ${reads_bam} \\
        --barcodes ${barcodes} \\
        --biosample ${biosample} \\
        --mismatches ${mismatches} \\
        --max-bc-dist ${max_bc_dist} \\
        --cpus ${task.cpus} \\
        --chunk-size ${params.chunk_size} \\
        --suffix '${unique_suffix}' \\
        ${trim_flag}

    cp .command.out demux.out

    # Rename audit log with unique suffix
    mv audit_log.csv audit${unique_suffix}.csv

    # Gzip all FASTQs
    for fq in *.fastq; do
        pigz -p ${task.cpus} "\${fq}" 2>/dev/null || gzip "\${fq}"
    done
    """
}

process FASTQC {
    tag "${fastq.simpleName}"
    label 'process_low'
    container "docker://biocontainers/fastqc:v0.11.9_cv8"

    publishDir "${params.outdir}/fastqc/mm${mismatches}_bc${max_bc_dist}", mode: 'copy'

    input:
    tuple val(mismatches), val(max_bc_dist), path(fastq)

    output:
    tuple val(mismatches), val(max_bc_dist), path("*.html"), emit: html
    tuple val(mismatches), val(max_bc_dist), path("*.zip"),  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${fastq}
    """
}

process MULTIQC {
    tag "multiqc_mm${mismatches}_bc${max_bc_dist}"
    label 'process_low'
    container "docker://multiqc/multiqc:pdf-latest"

    publishDir "${params.outdir}/multiqc/mm${mismatches}_bc${max_bc_dist}", mode: 'copy'

    input:
    tuple val(mismatches), val(max_bc_dist), path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc . --force
    """
}

process SWEEP_SUMMARY {
    tag "sweep_summary"
    label 'process_medium'
    container "${params.demux_container}"

    publishDir "${params.outdir}/sweep_summary", mode: 'copy'

    // Audit logs and read-length summaries are passed as two SEPARATE
    // path inputs rather than mixed into one channel. Mixing two collected
    // channels caused Nextflow to count each file twice in the staging map
    // (input file name collision), even though the filenames are unique.
    input:
    path audit_csvs
    path rl_csvs
    path biosample

    output:
    path "sweep_summary.csv",                  emit: csv
    path "sweep_summary_read_lengths.csv",     emit: rl_csv
    path "sweep_summary_read_lengths_detail.csv", emit: rl_detail_csv
    path "sweep_summary_umi_sharing.csv",      emit: umi_csv
    path "barcode_sharing_summary_mm*_bc*.csv", emit: per_combo_umi_csv, optional: true
    path "sweep_summary.html",                 emit: html

    script:
    """
    python ${projectDir}/bin/sweep_summary.py \\
        --input_dir . \\
        --biosample ${biosample} \\
        --output_prefix sweep_summary
    """
}

process FILTER_SHARED_UMIS {
    tag "filter_mm${mm}_bc${bcd}"
    label 'process_low'
    container "${params.demux_container}"

    publishDir "${params.outdir}/filtered_demux/mm${mm}_bc${bcd}", mode: 'copy'

    input:
    tuple val(mm), val(bcd), path(sharing_csv), path(fastqs)

    output:
    tuple val(mm), val(bcd), path("output/*.fastq.gz"), emit: filtered_fastqs, optional: true
    tuple val(mm), val(bcd), path("filter_log_mm${mm}_bc${bcd}.txt"), emit: log

    script:
    """
    set -euo pipefail

    # Extract sample names from the sharing CSV header (count_<sample> columns).
    # Sample names from sweep_summary.py are of the form: <bam>__<biosample>
    SAMPLES=\$(head -1 ${sharing_csv} | tr ',' '\\n' | grep '^count_' | sed 's/^count_//')

    mkdir -p input output

    # Each <bam>__<biosample> sample maps back to one SIMPLE_DEMUX output file
    # named <biosample>_<bam>_mm${mm}_bc${bcd}.fastq.gz. Split on the LAST '__'
    # (BAM names contain single underscores but not double underscores).
    for sample in \$SAMPLES; do
        bam_part="\${sample%__*}"     # before last __  → BAM name
        bs_part="\${sample##*__}"     # after  last __  → biosample
        src="\${bs_part}_\${bam_part}_mm${mm}_bc${bcd}.fastq.gz"
        if [ -f "\$src" ]; then
            zcat "\$src" > "input/\${sample}_mm${mm}_bc${bcd}.fastq"
        else
            echo "WARN: missing fastq for sample '\${sample}' (looked for \${src})" >&2
        fi
    done

    python ${projectDir}/bin/filter_shared_umis.py \\
        ${sharing_csv} \\
        _mm${mm}_bc${bcd} \\
        input \\
        output \\
        --next_max_thresh ${params.filter_next_max_thresh} \\
        --multiple_thresh ${params.filter_multiple_thresh} \\
        > filter_log_mm${mm}_bc${bcd}.txt 2>&1

    # gzip filtered outputs
    for fq in output/*.fastq; do
        [ -f "\$fq" ] || continue
        pigz -p ${task.cpus} "\$fq" 2>/dev/null || gzip "\$fq"
    done
    """
}

// ============================================================================
// Workflow
// ============================================================================
workflow {
    // Input channels
    //ch_bams = Channel.fromPath("${params.input_bams}/**/*.bam", checkIfExists: false)
    //            .filter { it.exists() }
    //            .unique()
    
    ch_barcodes = Channel.fromPath(params.barcodes, checkIfExists: true).first()

    // ch_inputs: tuple(id, bam, per-BAM biosample csv)
    // ch_biosample: one CSV for SWEEP_SUMMARY (3rd column Bam = id in samplesheet mode)
    if (params.samplesheet) {
        def sheet = file(params.samplesheet, checkIfExists: true)
        def rows  = sheet.splitCsv(header: true).collect { r ->
            r.collectEntries { k, v -> [(k?.trim()): v?.trim()] }
        }
        ['bam', 'sample_name', 'sample_barcode'].each { c ->
            if (!rows || !rows[0].containsKey(c)) {
                error "ERROR: samplesheet ${sheet} is missing column '${c}' (needs bam,sample_name,sample_barcode)"
            }
        }
        rows = rows.findAll { it.bam }
        // Accept 'RACE_TSO_5p--bc1001' (lima-style) or 'bc1001'; the demux matches on 'bc1001'
        rows.each { r -> r.sample_barcode = r.sample_barcode.tokenize('-').last() }

        def bs_dir = file("${workflow.workDir}/biosample_sheets")
        bs_dir.mkdirs()
        def inputs = []
        def combined = ['SampleName,BarcodeName,Bam']
        rows.groupBy { it.bam }.each { bam_path, rs ->
            def bam = file(bam_path, checkIfExists: true)
            def id  = bamId(bam)
            def map = [:]
            rs.each { r ->
                if (map.containsKey(r.sample_barcode) && map[r.sample_barcode] != r.sample_name) {
                    error "ERROR: ${r.sample_barcode} maps to both '${map[r.sample_barcode]}' and '${r.sample_name}' in ${bam_path}"
                }
                map[r.sample_barcode] = r.sample_name
            }
            def bs = bs_dir.resolve("${id}.biosample.csv")
            bs.text = (['SampleName,BarcodeName'] + map.collect { bc, s -> "${s},${bc}" }).join('\n') + '\n'
            map.each { bc, s -> combined << "${s},${bc},${id}" }
            inputs << tuple(id, bam, bs)
            log.info "  ${id}: ${map.size()} samples"
        }
        def ids = inputs.collect { it[0] }
        if (ids.size() != ids.unique(false).size()) {
            error "ERROR: two BAMs in the samplesheet share the same file name: ${ids}"
        }
        def all_bs = bs_dir.resolve('all_samples.biosample.csv')
        all_bs.text = combined.join('\n') + '\n'

        ch_inputs    = Channel.fromList(inputs)
        ch_biosample = Channel.value(all_bs)
    } else {
        def bam_input = params.input_bams
        def bam_glob = file(bam_input).isDirectory() ? "${bam_input}/**/*.bam" : bam_input
        ch_bams = Channel.fromPath(bam_glob, checkIfExists: false)
        ch_biosample = Channel.fromPath(params.biosample, checkIfExists: true).first()
        ch_inputs    = ch_bams.unique().combine(ch_biosample).map { bam, bs -> tuple(bamId(bam), bam, bs) }
    }

    // Step 1: Deconcatenate Kinnex arrays (once per BAM), or use HiFi reads as-is
    if (is_kinnex) {
        ch_mas_primers = Channel.fromPath(params.mas_primers, checkIfExists: true).first()
        SKERA_SPLIT(ch_inputs, ch_mas_primers)
        ch_reads = SKERA_SPLIT.out.skera_bam
    } else {
        ch_reads = ch_inputs
    }

    // Step 2: Build parameter sweep grid.
    // Dedupe at the Groovy list level so a CLI like `--max_bc_dist_range 0,0`
    // doesn't fan the sweep out into two identical (mm, bcd) combos — when
    // cached, both branches resolve to the same SIMPLE_DEMUX work-dir, and
    // SIMPLE_DEMUX.out.audit_log then emits the same file path twice, blowing
    // up SWEEP_SUMMARY with an input-staging collision.
    ch_mismatches  = Channel.from(params.mismatches_range .tokenize(',').collect { it.trim().toInteger() }.unique())
    ch_max_bc_dist = Channel.from(params.max_bc_dist_range.tokenize(',').collect { it.trim().toInteger() }.unique())

    // Combine: every input BAM × every (mismatches, max_bc_dist) pair
    ch_sweep = ch_reads
        .combine(ch_mismatches)
        .combine(ch_max_bc_dist)
        .map { id, bam, bs, mm, bcd -> tuple(id, bam, bs, mm, bcd) }

    // Step 3: Demux for each parameter combo (also produces read length reports)
    SIMPLE_DEMUX(ch_sweep, ch_barcodes)

    // Step 4: FastQC per individual FASTQ - names are now globally unique
    ch_fq_flat = SIMPLE_DEMUX.out.fastqs
        .flatMap { mm, bcd, fqs ->
            fqs instanceof List
                ? fqs.collect { fq -> tuple(mm, bcd, fq) }
                : [ tuple(mm, bcd, fqs) ]
        }
    FASTQC(ch_fq_flat)

    // Step 5: MultiQC per parameter combo
    ch_mqc = FASTQC.out.zip
        .map { mm, bcd, zip -> tuple("${mm}_${bcd}", mm, bcd, zip) }
        .groupTuple(by: 0)
        .map { key, mms, bcds, zips -> tuple(mms[0], bcds[0], zips.flatten()) }
    MULTIQC(ch_mqc)

    // Step 6: Sweep summary across all audit logs + read length summaries.
    // Pass audit CSVs and read-length CSVs as two separate `path` inputs.
    // `.unique()` before `.collect()` guards against duplicate path emissions
    // (e.g. when a user accidentally writes `--max_bc_dist_range 0,0`, both
    // branches resolve to the same cached work-dir under `-resume`).
    ch_audit_csvs = SIMPLE_DEMUX.out.audit_log
        .map { mm, bcd, csv -> csv }
        .unique()
        .collect()
    ch_rl_csvs    = SIMPLE_DEMUX.out.rl_csv
        .map { mm, bcd, csv -> csv }
        .unique()
        .collect()

    SWEEP_SUMMARY(ch_audit_csvs, ch_rl_csvs, ch_biosample)

    // Step 7: filter UMIs shared across samples per parameter combo.
    // Uses the per-combo barcode_sharing_summary CSVs emitted by SWEEP_SUMMARY.
    if (params.run_filter_shared_umis) {
        // Per-combo sharing CSVs: barcode_sharing_summary_mm{X}_bc{Y}.csv
        ch_per_combo_csv = SWEEP_SUMMARY.out.per_combo_umi_csv
            .flatten()
            .map { csv ->
                def m = (csv.name =~ /barcode_sharing_summary_mm(\d+)_bc(\d+)\.csv/)
                tuple(m[0][1].toInteger(), m[0][2].toInteger(), csv)
            }

        // Group all SIMPLE_DEMUX fastqs by (mm, bcd) across input BAMs
        ch_per_combo_fastqs = SIMPLE_DEMUX.out.fastqs
            .map { mm, bcd, fqs ->
                tuple("${mm}_${bcd}", mm, bcd, fqs instanceof List ? fqs : [fqs])
            }
            .groupTuple(by: 0)
            .map { key, mms, bcds, fq_lists ->
                tuple(mms[0], bcds[0], fq_lists.flatten())
            }

        ch_filter_in = ch_per_combo_csv
            .join(ch_per_combo_fastqs, by: [0, 1])
            .map { mm, bcd, csv, fastqs -> tuple(mm, bcd, csv, fastqs) }

        FILTER_SHARED_UMIS(ch_filter_in)
    }
}
