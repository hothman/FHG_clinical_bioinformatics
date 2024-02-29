#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/samplesheet.csv"
params.Reference = "$projectDir/Ref_genome/*.fa"
params.outdir = "$projectDir/Results/"
params.cpus = 2

log.info """\
    =========================================
                  *** ALIGN *** 
    =========================================
    reads          : ${params.reads}  
    Reference      : ${params.Reference}
    Indexes        : ${params.outdir}/BWA/Indexes
    Alignment file : ${params.outdir}/BWA/Mem
    """
    .stripIndent()

// Check Indexing the Reference
 
process createIndex {
    conda 'bioconda::bwa-mem2=0.7.17'
    conda "bioconda::gatk4=4.4"
    tag "CREATING INDEX AND DICT FILES FOR REF GENOME"
    publishDir "${params.outdir}/BWA/Indexes", mode: 'copy'

    input:
        path ref

    output:
        path "*", emit: "indexes"
        path "*.dict" 
        
    script:
        """
        bwa-mem2 index ${ref} ${params.outdir}
        gatk CreateSequenceDictionary REFERENCE=${ref} \\
                                  OUTPUT=${ref}.dict
        """
}

// Alignment based-reference


process alignReadsToRef {
    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19"
    tag "ALIGNING GENOMES TO REFERENCE"
    publishDir "${params.outdir}/BWA/Mem", mode: 'copy'

    input:

        path refGenome
        path idx
        tuple val(patient_id), path(R1), path(R2)

   output:
        path "${patient_id}_sorted.bam" , emit : bam


    script:
        """
        bwa-mem2 mem -t ${params.cpus}  -M ${refGenome} ${R1} ${R2} \\
                         | samtools view -Sb -@ ${params.cpus} | samtools sort -@ ${params.cpus}  -o ${patient_id}_sorted.bam  
                         
                         
        """
        

}

process assignReadGroup {
    conda "bioconda::picard=3.1.1"
    tag "ASSIGNING RG"
    publishDir "${params.outdir}/BWA/Mem", mode: 'copy'
    input:
        path aligned_bam
        tuple val(patient_id), path(R1), path(R2)
         
    output:
        path "${patient_id}_sorted_labeled.bam"

    script: 
    """
picard AddOrReplaceReadGroups \\
                         -I ${aligned_bam} \\
                         -O ${patient_id}_sorted_labeled.bam \\
                         -RGID ${patient_id} \\
                         -RGLB unspec \\
                         -RGPL ILLUMINA \\
                         -RGPU unspec \\
                         -RGSM 20 \\
                         -RGPM unspec \\
                         -RGCN unspec
                            """
    

}


process markDuplicates {
    conda "bioconda::gatk4=4.4"
    tag "MARK DUPLICATES"
    publishDir "${params.outdir}/BWA/Mem", mode: 'copy'
    input: 
        path sorted_bam
        tuple val(patient_id), path(R1), path(R2)
    
    output:
        path "*"

    script:
    """
    picard MarkDuplicates --INPUT $sorted_bam \\
                        --OUTPUT ${patient_id}_sorted_markduplicates.sam \\
                        --METRICS_FILE ${patient_id}.metrict \\
                        --TMP_DIR .
    """       

}
// Channels 

// channels 
READS = Channel.fromPath(params.reads, checkIfExists: true)  
               .splitCsv(header: true)  
               .map { row -> tuple(row.patient_id, file(row.R1), file(row.R2)) } 
  
Channel.fromPath(params.Reference)
       .first()
       .set { ref_file }
 
    
// Main workflow

 workflow {
      createIndex(ref_file)
      alignReadsToRef(ref_file, createIndex.out.indexes.collect(), READS)
      assignReadGroup(alignReadsToRef.out.collectFile( sort: true), READS)
      markDuplicates(assignReadGroup.out.collectFile( sort: true), READS)
}


