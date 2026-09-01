include { samplesheetToList } from 'plugin/nf-schema'

workflow {
    def clinical_samplesheet = params.input ? samplesheetToList(params.input, "assets/schema_input.json") : []
    
    Channel.fromList(clinical_samplesheet)
        .map { row ->
            def meta = [
                id:        row.sample,
                timepoint: row.timepoint,
                treatment: row.treatment,
                replicate: row.replicate
            ]
            def fastqs = [ file(row.fastq_1), file(row.fastq_2) ]
            
            return tuple(meta, fastqs)
        }
        .set { validated_reads_ch }

    validated_reads_ch
        .map { meta, fastqs ->
            def group_key = "Day${meta.timepoint}_${meta.treatment}"
            return tuple(group_key, meta, fastqs)
        }
        .groupTuple(by: 0)
        .map { group_key, meta_list, fastq_list ->
            log.info "Successfully grouped ${group_key} with ${meta_list.size()} replicates."
            return tuple(group_key, fastq_list)
        }
        .set { grouped_clinical_ch }
}