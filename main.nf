include { samplesheetToList } from 'plugin/nf-schema'
include { RNASEQ }          from './submodules/rnaseq/workflows/rnaseq.nf'

workflow {
    def clinical_samplesheet = params.input ? samplesheetToList(params.input, "assets/schema_input.json") : []
    
    Channel.fromList(clinical_samplesheet)
        .map { row ->
            def meta = [
                id:         row.sample,
                single_end: false,
                timepoint:  row.timepoint,
                treatment:  row.treatment
            ]
            def fastqs = [ file(row.fastq_1), file(row.fastq_2) ]
            
            return tuple(meta, fastqs)
        }
        .set { ch_reads }

    RNASEQ(ch_reads)
}