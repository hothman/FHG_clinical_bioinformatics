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
        params.reads = "$projectDir/Reads/*{1,2}.{fq,fastq,fq.gz,fastq.gz}"
     } 

 // ERROR Message     
       else { error "No input files with extensions .fq, .fastq, .fq.gz, or .fastq.gz found in $projectDir"}



params.outdir = "$projectDir/Results"
 
log.info """\
======================================================
           *** QUALITY CONTROL CHECK *** 
======================================================
reads           :${params.reads}  
Trimmed reads   :${params.outdir}/Trimmed_reads
Raw Qual Ctrl   :${params.outdir}/Reads_QC/RAW/  
Multiqc Raw     :${params.outdir}/Reads_QC/RAW/multiqc/
Trim Qual Ctrl  :${params.outdir}/Reads_QC/Trimmed/ 
Multiqc Trimmed :${params.outdir}/Reads_QC/Trimmed/multiqc/
    """
    .stripIndent()    
 
 
/*
 * Check Raw reads Quality ;
 */
 
process QC_RAW {
        publishDir "${params.outdir}/Reads_QC/RAW" ,  mode:'copy'
        
        input:
             tuple val(sample_id), path(reads) 

        output:
             file "*.{zip,html}"  

        script:
        
        """
        fastqc  ${reads} 
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
        multiqc ${fastqc} .
        """
}

/*
 * Trimming using fastq ;
 */
 
process Trimming {
        publishDir "${params.outdir}/Trimmed_reads", mode: 'copy'
        
        input:
             tuple val(name), path(reads)
        output:
             path (reads) 
        
        script:
        
        """        
        fastp -i ${reads[0]} -I ${reads[1]} -o ${reads[0].baseName}.gz -O ${reads[1].baseName}.gz
        """
}
 
/*
 * Trimmed reads quality control ;
 */
 
process QC_Trimmed {
        publishDir "${params.outdir}/Reads_QC/Trimmed" ,  mode:'copy'
      
        input:
             path (reads)  
        
        output:
             file "*.{zip,html}"  
             
        script:
        
        """
        fastqc ${reads} 
        """
}


process Multiqc_Trim {
        publishDir "${params.outdir}/Reads_QC/Trimmed/multiqc/" ,  mode:'copy'
        
        input:
            path (fastqc)  
            
        output:
           path "{multiqc_data,multiqc_report.html}" 
           
        script:
        """
        multiqc ${fastqc} .
        """
}


/* 
 * Channels ;
 */

Channel.fromFilePairs(params.reads, checkIfExists: true)
       .set { reads_pairs_channels }   
  
    
workflow {
   
    QC_RAW(reads_pairs_channels)   
    Multiqc_Raw(QC_RAW.out.collect())
    Trimming(reads_pairs_channels)
    QC_Trimmed(Trimming.out.collect())
    Multiqc_Trim(QC_Trimmed.out.collect())

}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
