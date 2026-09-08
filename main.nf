include { samplesheetToList }           from 'plugin/nf-schema'
include { FASTQC as FASTQC_RAW }        from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED }    from './modules/nf-core/fastqc/main'
include { FASTP }                       from './modules/nf-core/fastp/main'
include { PARABRICKS_RNAFQ2BAM }        from './modules/nf-core/parabricks/rnafq2bam/main'
include { SUBREAD_FEATURECOUNTS }       from './modules/nf-core/subread/featurecounts/main'
include { MULTIQC }                     from './modules/nf-core/multiqc/main'


workflow {
    if (!workflow.profile) { 
        exit 1, "ERROR: Please specify an execution profile (-profile standard or -profile portable)." 
    }
    if (!params.input) { exit 1, "ERROR: --input parameter is missing." }
    if (!params.fasta || !params.star_index || !params.gtf) { 
    exit 1, "ERROR: --fasta, --star_index, and --gtf are required." 
    }

    def clinical_samplesheet = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    
    ch_reads = Channel.fromList(clinical_samplesheet)
        .map { sample, timepoint, treatment, replicate, strandedness, fastq_1, fastq_2 ->
            def meta = [
                sample:       sample,
                timepoint:    timepoint,
                treatment:    treatment,
                replicate:    replicate,
                strandedness: strandedness,
                id:           "${sample}_${treatment}_d${timepoint}_r${replicate}",
                single_end:   false
            ]
            
            return tuple(meta, [ file(fastq_1), file(fastq_2) ])
        }

    ch_fasta = Channel.value(tuple([id: 'genome'], file(params.fasta, checkIfExists: true)))
    ch_star_index = Channel.value(tuple([id: 'star_index'], file(params.star_index, checkIfExists: true)))
    ch_gtf = Channel.value(file(params.gtf, checkIfExists: true))

    FASTQC_RAW(ch_reads)
    
    ch_fastp_input = ch_reads.map { meta, reads -> tuple(meta, reads, []) }
    FASTP(ch_fastp_input, false, false, false)

    FASTQC_TRIMMED(FASTP.out.reads)

    PARABRICKS_RNAFQ2BAM(
        FASTP.out.reads, 
        ch_fasta, 
        ch_star_index,
        true, 
        false  
    )

    ch_featurecounts_input = PARABRICKS_RNAFQ2BAM.out.bam
        .combine(ch_gtf)

    SUBREAD_FEATURECOUNTS(ch_featurecounts_input)

    ch_count_matrix_input = SUBREAD_FEATURECOUNTS.out.counts
        .ifEmpty { error "No successful featureCounts outputs were available to merge." }
        .collect(flat: false)
        .map { count_tuples ->
            def count_files = count_tuples.collect { item -> item[1] }
            tuple(count_tuples, count_files)
        }
    
    ch_expected_meta = ch_reads.map { meta, reads -> meta }.collect()

    GENERATE_COUNT_MATRIX(ch_count_matrix_input, ch_expected_meta)
    
    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.map { meta, logs -> logs })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.map { meta, logs -> logs })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { meta, logs -> logs })
    ch_multiqc_files = ch_multiqc_files.mix(PARABRICKS_RNAFQ2BAM.out.log_final.map { meta, log -> log })
    ch_multiqc_files = ch_multiqc_files.mix(PARABRICKS_RNAFQ2BAM.out.qc_metrics.map { meta, metrics -> metrics })
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.map { meta, logs -> logs })
    
    MULTIQC(
        ch_multiqc_files.collect().map { files -> [[id: 'clinical_study'], files, [], [], [], []] }
    )
}


process GENERATE_COUNT_MATRIX {
    tag 'merge featureCounts outputs'
    publishDir "${params.outdir}/counts", mode: 'copy'
    label 'process_single'
    
    input:
    tuple val(count_tuples), path(count_files)
    val expected_meta
    
    output:
    path 'gene_counts.tsv'
    path 'coldata.tsv'
    path 'count_manifest.tsv'
    path 'dropout_report.tsv'
    
    script:
    def manifest_rows = count_tuples.collect { item ->
        def meta = item[0]
        def count_file = item[1]
        "${meta.id}\t${count_file.name}\t${meta.timepoint}\t${meta.treatment}\t${meta.replicate}"
    }.join('\n')

    def expected_rows = expected_meta.collect { meta ->
        "${meta.id}\t${meta.sample}\t${meta.timepoint}\t${meta.treatment}\t${meta.replicate}"
    }.join('\n')
    
    """
cat > count_manifest.tsv <<'EOF'
sample_id\tcount_file\ttimepoint\ttreatment\treplicate
${manifest_rows}
EOF

cat > expected_manifest.tsv <<'EOF'
sample_id\tsample\ttimepoint\ttreatment\treplicate
${expected_rows}
EOF

python3 - <<'PYTHON'
import pandas as pd
from pathlib import Path

manifest = pd.read_csv('count_manifest.tsv', sep='\t')
expected_df = pd.read_csv('expected_manifest.tsv', sep='\t')

survived_ids = set(manifest['sample_id'])
expected_df['status'] = expected_df['sample_id'].apply(
    lambda sid: 'SURVIVED' if sid in survived_ids else 'DROPOUT'
)

expected_df.to_csv('dropout_report.tsv', sep='\t', index=False)

metadata_columns = {'Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'}
dataframes = []

for row in manifest.itertuples(index=False):
    df = pd.read_csv(row.count_file, sep='\t', comment='#')
    count_columns = [col for col in df.columns if col not in metadata_columns]
    sample_counts = df.set_index('Geneid')[count_columns[0]].rename(row.sample_id)
    dataframes.append(sample_counts)

matrix = pd.concat(dataframes, axis=1)

if matrix.isna().any().any():
    raise ValueError("Gene index mismatch detected. Samples contain differing gene lists.")

matrix = matrix.reset_index()
sample_columns = [col for col in matrix.columns if col != 'Geneid']
matrix[sample_columns] = matrix[sample_columns].astype('int64')
matrix.to_csv('gene_counts.tsv', sep='\t', index=False)

coldata = manifest[['sample_id', 'timepoint', 'treatment', 'replicate']].copy()
coldata.set_index('sample_id', inplace=True)
coldata.to_csv('coldata.tsv', sep='\t')
PYTHON
    """

stub:
    def manifest_rows = count_tuples.collect { item ->
        def meta = item[0]
        def count_file = item[1]
        "${meta.id}\t${count_file.name}\t${meta.timepoint}\t${meta.treatment}\t${meta.replicate}"
    }.join('\n')

    def expected_rows = expected_meta.collect { meta ->
        "${meta.id}\t${meta.sample}\t${meta.timepoint}\t${meta.treatment}\t${meta.replicate}"
    }.join('\n')
    
    def survived_ids = count_tuples.collect { it[0].id }
    def gene_counts_header = "Geneid\t" + survived_ids.join('\t')
    def gene_counts_row = "ENSG00000000003\t" + count_tuples.collect { "100" }.join('\t')
    
    def coldata_rows = count_tuples.collect { item ->
        "${item[0].id}\t${item[0].timepoint}\t${item[0].treatment}\t${item[0].replicate}"
    }.join('\n')

    """
cat > count_manifest.tsv <<'EOF'
sample_id\tcount_file\ttimepoint\ttreatment\treplicate
${manifest_rows}
EOF

cat > expected_manifest.tsv <<'EOF'
sample_id\tsample\ttimepoint\ttreatment\treplicate
${expected_rows}
EOF

python3 - <<'PYTHON'
import pandas as pd
import sys

manifest = pd.read_csv('count_manifest.tsv', sep='\t')

expected_df = pd.read_csv('expected_manifest.tsv', sep='\t')
survived_ids = set(manifest['sample_id'])

expected_df['status'] = expected_df['sample_id'].apply(lambda sid: 'SURVIVED' if sid in survived_ids else 'DROPOUT')
expected_df.to_csv('dropout_report.tsv', sep='\t', index=False)
PYTHON

cat > gene_counts.tsv <<'EOF'
${gene_counts_header}
${gene_counts_row}
EOF

cat > coldata.tsv <<'EOF'
sample_id\ttimepoint\ttreatment\treplicate
${coldata_rows}
EOF
    """
}