#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

String extractPathWithoutFilename(String filePath) {
    def file = new File(filePath)
    def parent = file.parent
    return parent ?: '.'
}

params.samplesFile = './data/samplesheet.csv'
params.referenceGenome = './ref_genome/dpyd.fa'
params.indexing = true

params.referenceGenomeRoot = extractPathWithoutFilename(params.referenceGenome)

// Read the CSV file
samples = file(params.samplesFile).splitCsv(header: true)

// Define the pipeline process
process processFastq {
    input:
    tuple val(patient_id), file(read1), file(read2) 

    //output:
    //file("${patient_id}_processed.fastq") 

    script:
    """
    echo  $read1 $read2 -o patient_id
    """
}


process createIndex {
    publishDir params.referenceGenomeRoot, mode: 'copy', overwrite: false
    input:
       file(genome) 

    output:
       file("*") 

    script:
    """
    bwa index ${genome}
    """
}


// Run the pipeline
workflow {
     Channel.fromPath(params.samplesFile) | splitCsv(header:true) | \
     map { row-> tuple(row.patient_id, file(row.R1), file(row.R2)) } | \
     processFastq

    if (params.indexing == true) { 
        Channel.fromPath(params.referenceGenome) | createIndex
        println "Creating index for reference genome!"
    }




     

}


