#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

String extractPathWithoutFilename(String filePath) {
    def file = new File(filePath)
    def parent = file.parent
    return parent ?: '.'
}

params.samplesFile = './data/samplesheet.csv'
params.referenceGenome = '/media/houcem/theDrum/BILIM/github/FHG_clinical_bioinformatics/assembly_ngs/ref_genome/dpyd.fa'
params.indexing = true
params.cpus=2

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
    conda 'bioconda::bwa-mem2=0.7.17'
    conda "bioconda::gatk4=4.4"
    tag "CREATING INDEX AND DICT FILES FOR REF GENOME"
    publishDir params.referenceGenomeRoot, mode: 'copy', overwrite: false
    input:
       file(genome) 

    output:
       val("indexing complete"), emit: indexing_control
       file("*.dict")

    script:
    """
    bwa-mem2 index ${genome}
    gatk CreateSequenceDictionary REFERENCE=${genome} \\
                                  OUTPUT=${genome}.dict
    """
}

process alignReadsToRef {
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    tag "ALIGNING GENOMES TO REFERENCE"
    label 'process_high'
    input: 
        tuple val(patient_id), file(read1), file(read2)

    output: 
        file("*_sorted.bam")
        val(patient_id)

    script: 
    """
    bwa-mem2 mem  ${params.referenceGenome} $read1 $read2 -t $params.cpus \\
                         | samtools view -Sb -@ 2 | samtools sort -o ${patient_id}_sorted.bam
    """
}

process assignReadGroup {
    conda "bioconda::picard=3.1.1"
    tag "ASSIGNING RG"
    input:
        file(aligned_bam)
        val(patient_id)

    script: 
    """
    picard AddOrReplaceReadGroups I=$aligned_bam \\
                 O=${patient_id}_sorted_labeled.bam \\ 
                 RGID=$patient_id \\
                 RGPL=ILLUMINA \\     #plateform
                 RGLB=unspec \\       #lib id
                 RGPU=unspec \\       #plateform unit 
                 RGSM=20 \\           #read group sample name
                 RGPM=unspec \\       #plateform model
                 RGDT=unspec  \\      #run date
                 RGCN=unspec          #sequecing center
    """

}

process markDuplicates {
    conda "bioconda::gatk4=4.4"
    tag "MARK DUPLICATES"

    input: 
        file(sorted_bam)
        val(sample_id)

    script:
    """
    gatk MarkDuplicates --INPUT $sorted_bam \\
                        --OUTPUT ${sample_id}_sorted_markduplicates.sam \\
                        --METRICS_FILE ${sample_id}.metrict \\
                        --TMP_DIR .
    """       

}

// Run the pipeline
workflow {
    // parsing sample datasheet
    sample = Channel.fromPath(params.samplesFile) | splitCsv(header:true) | \
        map { row-> tuple(row.patient_id, file(row.R1), file(row.R2)) } 

    // create index and align to reference 
    if (params.indexing == true) {
        ref_gen_channel = Channel.fromPath(params.referenceGenome)
        createIndex( ref_gen_channel)
        alignReadsToRef(sample) | markDuplicates

        
    }

    else {
       alignReadsToRef(sample)

    }

}


