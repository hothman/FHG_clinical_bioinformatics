#!/usr/bin/env nextflow



/*
 * NXF ver greater 19.08+ needed because of the use of tuple instead of set
 */

    if( !nextflow.version.matches('>=19.08') ) {
          println "This workflow requires Nextflow version 19.08 or greater and you are running version $nextflow.version"
          exit 1 }



/*
 * Parameters 
 */ 
 
  // Check if files with specified extensions exist in the project directory 
  
   if (file("$projectDir/Reads/*.{fq,fastq,fq.gz,fastq.gz}")) {
  
  // INPUT Parameters
        params.reads = "$projectDir/Reads/*_{1,2}.{fq,fastq,fq.gz,fastq.gz}"
     } 

 // ERROR Message     
       else { error "No input files with extensions .fq, .fastq, .fq.gz, or .fastq.gz found in $projectDir"}



params.outdir = "$projectDir/Results"
 
log.info """\
======================================================
      *** TRIMMING &  QUALITY CONTROL CHECK *** 
======================================================
reads             :${params.reads}  
Trimmed Reads     :${params.outdir}/trimmed
Trimmed Qual Ctrl :${params.outdir}/Reads_QC/trimmed/  
Multiqc Trimmed   :${params.outdir}/Reads_QC/trimmed/multiqc/
    """
    .stripIndent()    
 
 
/*
 * Trimming by Trimmomatic ;
 */
 
 // Default parameters , but could adjusted and specified on the command line " --something value" (eg. --LEADING 10) 

params.LEADING = 15
params.TRAILING = 15
params.SLIDINGWINDOW = '4:20'
params.MINLEN = 100


process Trimming {
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    
    input:
        tuple val(name), path(reads)
    output:
        path "*_paired.fq", emit: triMmed
        path "*_unpaired.fq"
        
    script:
        """        
        trimmomatic PE \
        ${reads} \
        ${reads[0].baseName}_paired.fq \
        ${reads[0].baseName}_unpaired.fq \
        ${reads[1].baseName}_paired.fq \
        ${reads[1].baseName}_unpaired.fq \
        LEADING:${params.LEADING} \
        TRAILING:${params.TRAILING} \
        SLIDINGWINDOW:${params.SLIDINGWINDOW} \
        MINLEN:${params.MINLEN}
        """
}


/*
 * Check Trimmed reads Quality ;
 */
 
process QC_Trim {
    publishDir "${params.outdir}/Reads_QC/trimmed", mode: 'copy'
        
        input:
             path(reads)

        output:
             path "*.{html,zip}"  

        script:
        
        """
        fastqc ${reads} 
        """
}


process Multiqc_Trim {
        publishDir "${params.outdir}/Reads_QC/trimmed/multiqc/" ,  mode:'copy'
        
        input:
            path (fastqc)  
            
        output:
            path "{multiqc_data,multiqc_report.html}" 
           
        script:
        """
        multiqc .
        """
}

/* 
 * Channels ;
 */

Channel.fromFilePairs(params.reads, checkIfExists: true)
       .set { reads_ch }   
  
    
workflow{
    Trimming(reads_ch)   
    QC_Trim(Trimming.out.triMmed)   
    Multiqc_Trim(QC_Trim.out.collect()) 
  
}
 
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}


