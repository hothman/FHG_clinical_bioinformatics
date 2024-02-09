#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/sample.csv"
params.Reference = "$projectDir/Ref_genome/*.fa"
params.outdir = "$projectDir/Results/"

log.info """\
    =========================================
             *** ALIGN & BWA MEM *** 
    =========================================
    reads        : ${params.reads}  
    Reference    : ${params.Reference}
    Indexes      : ${params.outdir}/BWA/Indexes
    Sam file     : ${params.outdir}/BWA/Mem
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
        path "*", emit: 'indexes'
        path "*.dict" 
        
    script:
        """
        bwa index ${ref} ${params.outdir}/BWA/BWAIndex
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
        tuple val(sampleId), path(read1), path(read2)

   output:
        path "alignment_${sampleId}.sam"

    script:
        """
        bwa mem -t ${task.cpus} -M ${refGenome} ${read1} ${read2} \\
                         | samtools view -Sb -@ 2 | samtools sort -o alignment_${sampleId}.sam  
        """
}

// channels 
READS = Channel.fromPath(params.reads, checkIfExists: true)  
               .splitCsv(header: true)  
               .map { row -> tuple(row.sampleId, file(row.read1), file(row.read2)) } 
  
Channel.fromPath(params.Reference)
       .first()
       .set { ref_file }

// Main workflow

 workflow {
      createIndex(ref_file)
      alignReadsToRef(ref_file, createIndex.out.indexes , READS)
}


