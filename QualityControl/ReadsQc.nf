#!/usr/bin/env nextflow

/*
_________                  .__    .__  ___________           _________                     __                   .__   
\\_____  \\  __ ___ ____   |  |   |__| \\__    ___/__ .__.  \\_   ___ \\  ____    ____   _/  |________   _____  |  |  
//  / \\  \\|  |  \\  _ \\ |  |   |  |   |    |   /  \|  |  /    \\  \\/ /  _ \\ /    \\ \\   __\\___\\ /  _  \\|  |  
//  \\_//  \\  |  // /_\ \\|  |__ |  |   |    |   \\___  |  \\     \\___(  (_) )/   |  \\ |  |  |  |   (  (_)  )|  |__
\\_____\\ \\_/____/(____  /|____/ |__|   |____|   / ____ |   \\______  /\\____/ |__/|  /  |__|  |__|   \\_____/ |____/
        \__/           \/                         \/                 \/              \/                               

*/ 


// NXF ver greater 19.08+ needed because of the use of tuple instead of set

    if( !nextflow.version.matches('>=19.08') ) {
          println "This workflow requires Nextflow version 19.08 or greater and you are running version $nextflow.version"
          exit 1 }



// Parameters 
 
params.reads = "$projectDir/samplesheet.csv"
params.outdir = "$projectDir/outdir"
 
log.info """\
======================================================
           *** QUALITY CONTROL CHECK *** 
======================================================
reads           :${params.reads}  
Raw Qual Ctrl   :${params.outdir}/QualityControl/RAW/  
Multiqc Raw     :${params.outdir}/QualityControl/RAW/multiqc/
    """
    .stripIndent()    
 
 
// Check Raw reads Quality ;
                  
process RawQc {
	conda 'bioconda::fastqc=0.12.1'
    	tag "GENERATING QUALITY FOR RAW READS"
        publishDir "${params.outdir}/QualityControl/RAW/" ,  mode:'copy'
        
        input:
        	tuple val(patient_id), path(R1), path(R2)

        output:
		path "*.{html,zip}"  

        script:
        """
        fastqc ${R1} ${R2} 
        """
}                             

// Gathering Qc Reports ;

process ReadsMultiqc {
	conda 'bioconda::multiqc=1.17'

    	tag "Gathering Multiqc FOR RAW READS"
        publishDir "${params.outdir}/QualityControl/RAW/multiqc/" ,  mode:'copy'
        
        input:
            path (fastqc)  
            
        output:
            path "{multiqc_data,multiqc_report.html}" 
           
        script:
        """
        multiqc .
        """
}


// Channels ;

Channel.fromPath(params.reads, checkIfExists: true)  
       .splitCsv(header: true)  
       .map { row -> tuple(row.patient_id, file(row.R1), file(row.R2)) } 
       .set{ READS }

// Workflow 
   
workflow{
	RawQc(READS)   
    	ReadsMultiqc(RawQc.out.collect())
  
}
 

