#!/usr/bin/env nextflow

// Params

params.reads = "$projectDir/SamplesheetForTrimming.csv"
params.outdir = "$projectDir/outdir"
params.cpus = 2

log.info """\
======================================================
      *** TRIMMING &  QUALITY CONTROL CHECK *** 
======================================================
reads             :${params.reads}  
Trimmed Reads     :${params.outdir}/TrimmedREADS/
Trimmed Qual Ctrl :${params.outdir}/QualityControl/TRIMMED/  
Multiqc Trimmed   :${params.outdir}/QualityControl/TRIMMED/multiqc/
    """
    .stripIndent()    
 
 
//Trimming by Trimmomatic 

process Trimming {
	conda "bioconda::trimmomatic=0.39 "
    	tag "TRIMMIG READS"
	publishDir path: "${params.outdir}/TrimmedREADS/", mode: 'copy' 
   
    	input:
		tuple val(patient_id), path(R1), path(R2), val(MINLEN), val(LEADING), val(TRAILING), val(SLIDINGWINDOW)
    
    	output:
        	path "*.fastq"			, emit: paired  	// To be used in DOWNSTREAM Analysis
        	path "*_unpaired.fastq"		, emit: unpaired
        	
    	script:
        """        
        trimmomatic PE \\
       		-threads ${params.cpus} \\
       		${R1} ${R2} \\
       		${R1.baseName.takeWhile{ it != '.' }}.fastq \\
       		${R1.baseName.takeWhile{ it != '.' }}_unpaired.fastq \\
       	  	${R2.baseName.takeWhile{ it != '.' }}.fastq \\
       		${R2.baseName.takeWhile{ it != '.' }}_unpaired.fastq \\
       		MINLEN:${MINLEN} \\
       		LEADING:${LEADING} \\
       		TRAILING:${TRAILING} \\
        	SLIDINGWINDOW:${SLIDINGWINDOW}
        """
}

// Check Trimmed reads Quality ;
 
process TrimmedQC {
    	conda "bioconda::fastqc=0.12.1"
    	tag "CHECK TRIMMED READS QUALITY"
        publishDir "${params.outdir}/QualityControl/TRIMMED/", mode: 'copy'

        
        input:
        	path(reads)

        output:
             	path "*.{html,zip}"  

        script:
        """
        fastqc ${reads} 
        """
}

// Multiqc 

process MultiqcTrimmed {
	conda "bioconda::multiqc=1.17"
    	tag "GATHER TRIMMED QC REPORTS"
        publishDir "${params.outdir}/QualityControl/TRIMMED/multiqc/" ,  mode:'copy'

        
        input:
            	path (fastqc)  
            
        output:
            	path "{multiqc_data,multiqc_report.html}" 
           
        script:
        """
        multiqc .
        """
}

// Channels 
 
Channel.fromPath(params.reads, checkIfExists: true)  
       .splitCsv(header: true)  
       .map { row -> tuple(row.patient_id,
       		      file(row.R1), 
       		       file(row.R2), 
       		        row.MINLEN, 
       		         row.LEADING,
       		          row.TRAILING, 
       		           row.SLIDINGWINDOW ) 	}.set{ READS }
       
       
workflow {
    	Trimming	(	READS				) 
    	TrimmedQC	(	Trimming.out.paired 		)   
    	MultiqcTrimmed	(	TrimmedQC.out.collect() 	) 
}
 
 
 
 
