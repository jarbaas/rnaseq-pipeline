nextflow.enable.dsl=2

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
    
    Channel.fromList(clinical_samplesheet)
        .map { sample, timepoint, treatment, replicate, fastq_1, fastq_2 ->
            def sample_id = "${sample}_${treatment}_d${timepoint}_r${replicate}"
            
            def meta = [
                id:         sample_id,
                sample:     sample,
                single_end: false,
                timepoint:  timepoint,
                treatment:  treatment,
                replicate:  replicate
            ]
            
            return tuple(meta, [ file(fastq_1), file(fastq_2) ])
        }
        .set { ch_reads }

    ch_fasta = Channel.value(tuple([id: 'genome'], file(params.fasta)))
    ch_star_index = Channel.value(tuple([id: 'star_index'], file(params.star_index)))
    ch_gtf = Channel.value(tuple([id: 'genome'], file(params.gtf)))

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
        .map { meta, bam, meta_gtf, gtf -> 
            tuple(meta, bam, gtf) 
        }

    SUBREAD_FEATURECOUNTS(ch_featurecounts_input)

    ch_count_matrix_input = SUBREAD_FEATURECOUNTS.out.counts
        .map { meta, counts -> tuple(meta, counts) }
        .collect(flat: false)
        .map { count_tuples ->
            def count_files = count_tuples.collect { item -> item[1] }
            tuple(count_tuples, count_files)
        }
    

    GENERATE_COUNT_MATRIX(ch_count_matrix_input)
    
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
    
    output:
    path 'gene_counts.tsv'
    path 'count_manifest.tsv'
    
    script:
    def manifest_rows = count_tuples.collect { item ->
        def meta = item[0]
        def count_file = item[1]
        "${meta.id}\t${count_file.getName()}"
    }.join('\n')
    
    """
cat > count_manifest.tsv <<'EOF'
sample_id\tcount_file
${manifest_rows}
EOF

python3 - <<'PYTHON'
import pandas as pd
from pathlib import Path

manifest = pd.read_csv('count_manifest.tsv', sep='\\t')

if manifest.empty:
    raise SystemExit('No successful featureCounts outputs were available to merge.')

for filename in manifest['count_file']:
    if not Path(filename).is_file():
        raise FileNotFoundError(f'Expected staged featureCounts result was not found: {filename}')

matrix = None
metadata_columns = {'Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'}

for row in manifest.itertuples(index=False):
    df = pd.read_csv(row.count_file, sep='\\t', comment='#')
    
    if 'Geneid' not in df.columns:
        raise ValueError(f'{row.count_file} is not a featureCounts table: Geneid is absent.')
        
    count_columns = [col for col in df.columns if col not in metadata_columns]
    
    if len(count_columns) != 1:
        raise ValueError(f'{row.count_file} has {len(count_columns)} count columns; this expects exactly one.')

    sample_counts = df[['Geneid', count_columns[0]]].rename(columns={count_columns[0]: row.sample_id})
    
    matrix = sample_counts if matrix is None else matrix.merge(sample_counts, on='Geneid', how='outer')

matrix = matrix.fillna(0)
sample_columns = [col for col in matrix.columns if col != 'Geneid']
matrix[sample_columns] = matrix[sample_columns].astype('int64')
matrix.to_csv('gene_counts.tsv', sep='\\t', index=False)
PYTHON
    """
}