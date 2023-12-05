#!/usr/bin/env nextflow

/*
_________                  .__    .__  ___________           _________                     __                   .__   
\\_____  \\  __ ___ ____   |  |   |__| \\__    ___/__ .__.  \\_   ___ \\  ____    ____   _/  |________   _____  |  |  
//  / \\  \\|  |  \\  _ \\ |  |   |  |   |    |   /  \|  |  /    \\  \\/ /  _ \\ /    \\ \\   __\\___\\ /  _  \\|  |  
//  \\_//  \\  |  // /_\ \\|  |__ |  |   |    |   \\___  |  \\     \\___(  (_) )/   |  \\ |  |  |  |   (  (_)  )|  |__
\\_____\\ \\_/____/(____  /|____/ |__|   |____|   / ____ |   \\______  /\\____/ |__/|  /  |__|  |__|   \\_____/ |____/
        \__/           \/                         \/                 \/              \/                               

*/ 


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
        params.reads = "$projectDir/Reads/*.{fq,fastq,fq.gz,fastq.gz}"
     } 

 // ERROR Message     
       else { error "No input files with extensions .fq, .fastq, .fq.gz, or .fastq.gz found in $projectDir"}



params.outdir = "$projectDir/Results"
 
log.info """\
======================================================
           *** QUALITY CONTROL CHECK *** 
======================================================
reads           :${params.reads}  
Raw Qual Ctrl   :${params.outdir}/Reads_QC/RAW/  
Multiqc Raw     :${params.outdir}/Reads_QC/RAW/multiqc/
    """
    .stripIndent()    
 
 
/*
 * Check Raw reads Quality ;
 */
 
process QC_RAW {
        publishDir "${params.outdir}/Reads_QC/RAW" ,  mode:'copy'
        
        input:
             path(reads) 

        output:
             path "*.{html,zip}"  

        script:
        
        """
        fastqc ${reads} 
        """
}


process Multiqc_Raw {
        publishDir "${params.outdir}/Reads_QC/RAW/multiqc/" ,  mode:'copy'
        
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

Channel.fromPath(params.reads, checkIfExists: true)
       .set { reads_ch }   
  
    
workflow{
    QC_RAW(reads_ch)   
    Multiqc_Raw(QC_RAW.out.collect())
  
}
 
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}


