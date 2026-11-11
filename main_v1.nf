#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// nf-flairr-demux: Kinnex FLAIRR-seq parameter sweep
// ============================================================================
// Input: give ONE of --samplesheet, --input_bams or --input_fastqs
params.samplesheet  = null    // CSV/TSV, one row per sample: sample name, barcode, FASTQ or BAM path
params.input_bams   = null    // HiFi BAM(s) or a folder of them (+ --biosample)
params.input_fastqs = null    // FASTQ(s) or a folder of them (+ --biosample)
params.biosample    = null    // SampleName,BarcodeName CSV (not needed with --samplesheet)
params.kinnex       = true    // BAMs are Kinnex arrays: run skera first. false = demux BAMs as they are
params.mas_primers  = null    // Kinnex MAS adapters for skera (only for BAMs with --kinnex true)
params.barcodes     = null
params.tso_window   = 1       // positions from the read start searched for the TSO (1 = position 0 only)
params.bc_window    = 1       // positions from the read end searched for the sample barcode
params.bc_trim_len  = 0       // trim only this many barcode bases (0 = whole barcode entry);
                              // e.g. 16 keeps the C-region primer of "bc1004_TRA"-style entries
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
def n_inputs = [params.samplesheet, params.input_bams, params.input_fastqs].count { it }
if (n_inputs == 0) { error "ERROR: give one of --samplesheet, --input_bams or --input_fastqs" }
if (n_inputs > 1)  { error "ERROR: give only one of --samplesheet / --input_bams / --input_fastqs" }
if (!params.samplesheet && !params.biosample) { error "ERROR: --biosample is required with --input_bams / --input_fastqs" }
if (!params.barcodes) { error "ERROR: --barcodes is required" }
def KINNEX = params.kinnex.toString().toBoolean()

// ============================================================================
// Processes
// ============================================================================
process SKERA_SPLIT {
    tag "${bam.baseName}"
    label 'process_medium'
    container "${params.demux_container}"

    publishDir "${params.outdir}/skera", mode: 'copy', pattern: '*.skera.bam*'

    input:
    tuple val(input_id), path(bam), val(biosample_csv)
    path mas_primers

    output:
    tuple val("${bam.baseName}.skera"), path("${bam.baseName}.skera.bam"), val(biosample_csv), emit: skera_bam
    path "${bam.baseName}.skera.bam.pbi", emit: skera_pbi, optional: true

    script:
    """
    skera split \\
        ${bam} \\
        ${mas_primers} \\
        ${bam.baseName}.skera.bam
    """
}

process SIMPLE_DEMUX {
    tag "${input_id}_mm${mismatches}_bc${max_bc_dist}"
    label 'process_medium'
    container "${params.demux_container}"

    publishDir "${params.outdir}/demux/mm${mismatches}_bc${max_bc_dist}/${input_id}", mode: 'copy'

    input:
    tuple val(input_id), path(reads), val(biosample_csv), val(mismatches), val(max_bc_dist)
    path barcodes

    output:
    tuple val(mismatches), val(max_bc_dist), path("*.fastq.gz"),            emit: fastqs
    tuple val(mismatches), val(max_bc_dist), path("audit_*.csv"),           emit: audit_log
    tuple val(mismatches), val(max_bc_dist), path("read_length_*.csv"),     emit: rl_csv
    tuple val(mismatches), val(max_bc_dist), path("read_length_*.html"),    emit: rl_html
    path "demux.out"

    script:
    def trim_flag = params.no_trim ? '--no-trim' : ''
    def unique_suffix = "_${input_id}_mm${mismatches}_bc${max_bc_dist}"
    """
    # barcode -> sample map for this input only
    cat > biosample.csv <<'BIOSAMPLE_EOF'
${biosample_csv.trim()}
BIOSAMPLE_EOF

    python ${projectDir}/bin/simple_demux.py \\
        --input ${reads} \\
        --barcodes ${barcodes} \\
        --biosample biosample.csv \\
        --tso-window ${params.tso_window} \\
        --bc-window ${params.bc_window} \\
        --bc-trim-len ${params.bc_trim_len} \\
        --mismatches ${mismatches} \\
        --max-bc-dist ${max_bc_dist} \\
        --cpus ${task.cpus} \\
        --chunk-size ${params.chunk_size} \\
        --suffix '${unique_suffix}' \\
        ${trim_flag}

    cp .command.out demux.out

    # Rename audit log with unique suffix
    mv audit_log.csv audit${unique_suffix}.csv

    # Gzip all FASTQs written by the demux (not the staged input, when it is a FASTQ)
    for fq in *.fastq; do
        [ "\${fq}" = "${reads}" ] && continue
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

    // Step 1: demux inputs as (input_id, reads, biosample_csv_text).
    // The id ends up in every output name, and the shared-UMI filter splits
    // "<id>__<sample>" on '__', so ids must not contain a double underscore.
    def check_id = { String id ->
        if (id.contains('__')) error "ERROR: input name '${id}' contains '__'; rename the file (single underscores are fine)"
        id
    }
    def FQ_RE  = /\.(fastq|fq)(\.gz)?$/
    def id_of  = { f -> check_id(f.name.replaceAll(FQ_RE, '').replaceAll(/\.bam$/, '')) }
    def is_bam = { f -> f.name.endsWith('.bam') }

    def inputs = []     // list of [id, file, biosample_csv_text]
    if (params.samplesheet) {
        // ---- samplesheet: one row per sample; rows that share a file are demuxed together
        def bc_names = file(params.barcodes).readLines()
            .findAll { it.startsWith('>') }
            .collect { it.substring(1).trim().split(/\s+/)[0] } as Set
        def lines  = file(params.samplesheet, checkIfExists: true).readLines()
            .collect { it.replaceAll(/\r$/, '') }
            .findAll { it.trim() && !it.startsWith('#') }
        def sep    = lines[0].contains('\t') ? '\t' : ','
        def header = lines[0].split(sep, -1).collect { it.trim() }
        def pick   = { List names -> names.find { header.contains(it) } }
        def c_sample = pick(['SampleName', 'sample_name', 'sample'])
        def c_bc     = pick(['BarcodeName', 'primer_barcode', 'barcode'])
        def c_reads  = pick(['fastq', 'fastq_path', 'reads', 'input', 'bam', 'bam_path'])
        def c_locus  = pick(['locus'])
        def c_keep   = pick(['sequenced', 'include'])
        if (!c_sample || !c_bc || !c_reads)
            error "ERROR: --samplesheet needs columns for sample (SampleName / sample_name / sample), " +
                  "barcode (BarcodeName / primer_barcode / barcode) and reads (fastq / fastq_path / bam / bam_path). Found: ${header}"
        log.info "Samplesheet columns: sample=${c_sample}, barcode=${c_bc}, reads=${c_reads}" +
                 (c_locus ? ", locus=${c_locus}" : "") + (c_keep ? ", keep=${c_keep}" : "")

        def rows = lines.drop(1).collect { l ->
            def v = l.split(sep, -1).collect { it.trim() }
            def r = [:]; header.eachWithIndex { h, i -> r[h] = i < v.size() ? v[i] : '' }; r
        }.findAll { r -> !c_keep || !(r[c_keep].toLowerCase() in ['no', 'n', 'false', '0']) }
        if (!rows) error "ERROR: no rows to run in ${params.samplesheet}"

        // barcode name as it appears in --barcodes: "<bc>_<locus>" if present there, else "<bc>"
        def resolve_bc = { r ->
            def bc = r[c_bc]; def loc = c_locus ? r[c_locus] : ''
            def cands = (loc ? ["${bc}_${loc}".toString()] : []) + [bc]
            def hit = cands.find { bc_names.contains(it) }
            if (!hit) error "ERROR: barcode ${cands.join(' / ')} (sample ${r[c_sample]}) not found in ${params.barcodes}"
            hit
        }
        rows.groupBy { it[c_reads] }.each { path, rs ->
            def f = file(path, checkIfExists: true)
            def map = rs.collect { r -> [r[c_sample], resolve_bc(r)] }
            def dup = map.groupBy { it[1] }.findAll { bc, ms -> ms*.getAt(0).unique().size() > 1 }
            if (dup) error "ERROR: ${f.name}: one barcode is listed for several samples: ${dup.keySet()}"
            inputs << [id_of(f), f, "SampleName,BarcodeName\n" + map.collect { "${it[0]},${it[1]}" }.join('\n')]
        }
        log.info "Samplesheet: ${rows.size()} samples in ${inputs.size()} input file(s)"
    } else {
        // ---- folder / glob inputs, one barcode map for all of them (--biosample)
        def bios_text = file(params.biosample, checkIfExists: true).text
        def src  = params.input_fastqs ?: params.input_bams
        def glob = params.input_fastqs ? '**.{fastq,fq,fastq.gz,fq.gz}' : '**.bam'
        def files = file(src).isDirectory() ? files("${src}/${glob}") : files(src)
        files.unique().each { f -> inputs << [id_of(f), f, bios_text] }
    }
    if (!inputs) error "ERROR: no input files found"
    def ids = inputs.collect { it[0] }
    def clash = ids.findAll { id -> ids.count(id) > 1 }.unique()
    if (clash) error "ERROR: input files share a name (${clash.join(', ')}); rename them so each is unique"

    // Kinnex BAMs are deconcatenated with skera first; FASTQs (and BAMs with --kinnex false) go straight in
    def to_skera = KINNEX ? inputs.findAll { is_bam(it[1]) } : []
    def direct   = inputs - to_skera
    if (to_skera && !params.mas_primers) error "ERROR: --mas_primers is required for Kinnex BAMs (or pass --kinnex false)"

    ch_direct = Channel.fromList(direct.collect { tuple(it[0], it[1], it[2]) })
    if (to_skera) {
        ch_mas_primers = Channel.fromPath(params.mas_primers, checkIfExists: true).first()
        SKERA_SPLIT(Channel.fromList(to_skera.collect { tuple(it[0], it[1], it[2]) }), ch_mas_primers)
        ch_demux_in = ch_direct.mix(SKERA_SPLIT.out.skera_bam.map { id, bam, bios -> tuple(check_id(id), bam, bios) })
    } else {
        ch_demux_in = ch_direct
    }

    // One combined (input, barcode) -> sample map for the sweep summary / UMI filter
    ch_biosample = ch_demux_in
        .map { id, reads, bios -> bios.readLines().drop(1).findAll { it.trim() }.collect { "${it},${id}" }.join('\n') + '\n' }
        .collectFile(name: 'biosample_all.csv', seed: 'SampleName,BarcodeName,Input\n', sort: true)
        .first()

    // Step 2: Build parameter sweep grid.
    // Dedupe at the Groovy list level so a CLI like `--max_bc_dist_range 0,0`
    // doesn't fan the sweep out into two identical (mm, bcd) combos — when
    // cached, both branches resolve to the same SIMPLE_DEMUX work-dir, and
    // SIMPLE_DEMUX.out.audit_log then emits the same file path twice, blowing
    // up SWEEP_SUMMARY with an input-staging collision.
    ch_mismatches  = Channel.from(params.mismatches_range.toString().tokenize(',').collect { it.trim().toInteger() }.unique())
    ch_max_bc_dist = Channel.from(params.max_bc_dist_range.toString().tokenize(',').collect { it.trim().toInteger() }.unique())

    // Combine: every input × every (mismatches, max_bc_dist) pair
    ch_sweep = ch_demux_in
        .combine(ch_mismatches)
        .combine(ch_max_bc_dist)
        .map { id, reads, bios, mm, bcd -> tuple(id, reads, bios, mm, bcd) }

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
