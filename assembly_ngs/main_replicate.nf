#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/samplesheet.csv"
params.Reference = "$projectDir/Ref_genome/Homo_sapiens_GRCh38_dna_chromosome20.fa"
params.outdir = "$projectDir/outdir/"
params.knwonSite1 = "$projectDir/knownsites/1000g_gold_standard.indels.filtered.vcf"
params.knwonSite2 = "$projectDir/knownsites/GCF.38.filtered.renamed.vcf"
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
 
process IndexRef {
	conda "bioconda::gatk4=4.4 bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19"
	tag "CREATING INDEX AND DICT FILES FOR REF GENOME"
        publishDir "${params.outdir}/Indexes/Reference", mode: 'copy'

     	input:
     		path ref
     		
   	output:
        	path "*.{0123,amb,ann,bwt.2bit.64,pac}" 	, emit: "idxREF"
        	path "*.dict"                           	, emit: "dicREF"
        	path "*.fai"                            	, emit: "samidxREF"
    	script:
        """
        bwa-mem2 index ${ref}  
        gatk CreateSequenceDictionary --REFERENCE ${ref}
        samtools faidx ${ref}  --output ${ref}.fai                         
        """
}

// Alignment based-reference


process alignReadsToRef {
	conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.19"
    	tag "ALIGNING GENOMES TO REFERENCE"
    	publishDir "${params.outdir}/Mapping", mode: 'copy'

    	input:
        	path refGenome
        	path idx
        	tuple val(patient_id), path(R1), path(R2)

   	output:
        	path "${patient_id}_sor.bam" , emit : sorted_bam //  Sorted Bam file
       
    	script:
        """
       	bwa-mem2 mem -t ${params.cpus}  -M ${refGenome} ${R1} ${R2} \
        		| samtools view -Sb -@ ${params.cpus} \
                        | samtools sort -@ ${params.cpus}  -o ${patient_id}_sor.bam           
        """
}

// Assigning ReadGroups

process assignReadGroup {
	conda "bioconda::picard=3.1.1"
    	tag "ASSIGNING READ GROUPS"
    	publishDir "${params.outdir}/Mapping", mode: 'copy'
    
    	input:
        	path aligned_bam
        	tuple val(patient_id), path(R1), path(R2)
         
    	output:
        	path "${patient_id}_sor@RG.bam" 	, emit : sorted_labeled_bam // Sorted_labeled bam file with RG

    	script: 
    	"""
	picard AddOrReplaceReadGroups \
        		-I ${aligned_bam} \
                        -O ${patient_id}_sor@RG.bam \
                        -RGID ${patient_id} \
                        -RGLB unspec \
                        -RGPL ILLUMINA \
                        -RGPU unspec \
                        -RGSM 20 \
                        -RGPM unspec \
                        -RGCN unspec
     	"""   
}

// Marking Duplicates

process markDuplicates {
	conda "bioconda::picard=3.1.1"
	tag "MARKING DUPLICATES"
    	publishDir "${params.outdir}/Mapping", mode: 'copy'
    
    	input: 
    	    path sorted_bam
    	    tuple val(patient_id), path(R1), path(R2)
    
    	output:
    	    path "${patient_id}_sor@RG@MD.bam" 	, emit : sorted_markduplicates_bam
    	    path "${patient_id}.metrict"

    	script:
    	"""
    	picard MarkDuplicates --INPUT $sorted_bam \
        		--OUTPUT ${patient_id}_sor@RG@MD.bam \
                        --METRICS_FILE ${patient_id}.metrict \
                        --TMP_DIR . 
    	"""       
}
// Generating Indexes of Bam files

process IndexBam{
    	conda "bioconda::samtools=1.19"
    	tag "CREATING INDEX FOR BAM FILES"
    	publishDir "${params.outdir}/Indexes/BamFiles", mode: 'copy'

    	input:
        	path BamFile 	  	// Bam file from sorted_markduplicates_bam 

    	output:
        	path "${BamFile}.bai"		, emit: "IDXBAM" 
    
    	script:
        """
        samtools index -@ ${params.cpus}  ${BamFile}  

       	"""
}

// Generate Statictics before & after Marking Duplicates

process GenerateStat {
    	conda "bioconda::samtools=1.19"
    	tag "Statistics for bam-sam files"
    	publishDir "${params.outdir}/Mapping/Metrics", mode: 'copy'
    	
    	input:
      		path sorted_labeled_bam
      		path sorted_markduplicates_bam
    
    	output:
     		path "*.flagstat" 

    	script: 
    	"""
      	samtools flagstat -@ {params.cpus} ${sorted_labeled_bam} > ${sorted_labeled_bam}.flagstat
      	samtools flagstat -@ {params.cpus}  ${sorted_markduplicates_bam} >  ${sorted_markduplicates_bam}.flagstat
    	"""
}

// Generating Indexes of known sites vcfs

process IndexKNownSites {
    	conda "bioconda::gatk4=4.4"
    	tag "CREATING INDEX for known sites vcf"
    	publishDir "${params.outdir}/Indexes/knownSites", mode: 'copy'

    	input:
        	path kn_site_File1 	  // known Sites file n°1 related to params.knwonSite1
        	path kn_site_File2 	  // known Sites file n°2 related to params.knwonSite2

    	output:
        	path "${kn_site_File1}.idx"	, emit: "IDXknS1" 
        	path "${kn_site_File2}.idx"     , emit: "IDXknS2"
    
    	script:
        """
        gatk IndexFeatureFile \
        		--input ${kn_site_File1} \
        		--output ${kn_site_File1}.idx\
        		--tmp-dir . 
        
       	gatk IndexFeatureFile \
        		--input ${kn_site_File2} \
        		--output ${kn_site_File2}.idx\
        		--tmp-dir . 
       	"""
}


// BaseRecalibration 

process BaseRecalibrator {
	conda "bioconda::gatk4=4.4"
    	tag "CREATING TABLE FOR BQSR"
    	publishDir "${params.outdir}/Mapping", mode: 'copy'

    	input:
       		path ref  
       		path dic
       		path fai
       		path sor_md_bam_file
       		path knownsiteFile1
       		path IDXknsF1 		  // index of known site file n° x
       		path knownsiteFile2
       		path IDXkGCF
      
    	output:
        	path "*bqsr.table" , emit: "BQSR_Table"

    	script:
        """
  	gatk BaseRecalibrator	\
  			--reference ${ref} \
  	  		--input ${sor_md_bam_file} \
  	  		--known-sites ${knownsiteFile1} \
  	  		--known-sites ${knownsiteFile2} \
  	  		--output ${sor_md_bam_file.baseName}.bqsr.table
       """
}

// Apply base Recalibration 

process ApplyBQSR {
	conda "bioconda::gatk4=4.4"
	tag "APPLYING BASE QUALITY SCORE RECALIBRATION"
    	publishDir "${params.outdir}/Mapping", mode: 'copy'

    	input:
       		path sor_md_bam_file 
       		path bqsrTABLE
     
    	output:
        	path "*.recal.bam" , emit: "recal_bam"

    	script:
        """
  	gatk ApplyBQSR \
  			--input ${sor_md_bam_file} \
  			--bqsr-recal-file ${bqsrTABLE} \
  			--output ${sor_md_bam_file.baseName}.recal.bam
       	"""
}  

// Generating Indexes of Recalibrated Bam files

process IndexRecalBam{
    	conda "bioconda::samtools=1.19"
    	tag "CREATING INDEX FOR Recalibrated BAM FILES"
    	publishDir "${params.outdir}/Indexes/BamFiles", mode: 'copy'

    	input:
        	path RecalBamFile	// Bam file from recal_bam

    	output:
        	path "${RecalBamFile}.bai"     	, emit: "IDXRECALBAM"
    
    	script:
        """
        samtools index -@ ${params.cpus}  ${RecalBamFile}  
       	"""
}

// GATK variant calling for raw mapped reads

process RawHaploCall {
	conda "bioconda::gatk4=4.4"
	tag "Raw Var Call with Gatk HaplotypeCaller"
    	publishDir "${params.outdir}/Mapping/Variants", mode: 'copy'

    	input:
       		path ref  
       		path dic
       		path fai
       		path BamFile
       		path BamBai
     
    	output:
        	path "*raw.HC.vcf" , emit: "vcf_HaplotypeCaller_Raw" // HC : HaplotypeCaller
    	script:
        """
  	gatk HaplotypeCaller \
			--reference ${ref} \
			--input ${BamFile} \
			--output ${BamFile.baseName}.raw.HC.vcf 
       	"""
}

// GATK variant calling for Recalibrated mapped reads  // SNP  // same processes will be replicated using include {NameProcess as NameProcess1} from './modules/NameProcess

process RecalHaploCall {
	conda "bioconda::gatk4=4.4"
	tag "Recal Var Call with Gatk HaplotypeCaller"
    	publishDir "${params.outdir}/Mapping/Variants", mode: 'copy'

    	input:
       		path ref  
       		path dic
       		path fai
       		path ReclBamFile
       		path ReclBamBai
     
    	output:
        	path "*.vcf" , emit: "vcf_HaplotypeCaller_Recal"

    	script:
        """
  	gatk HaplotypeCaller \
			--reference ${ref} \
			--input ${ReclBamFile} \
			--output ${ReclBamFile.baseName}.HC.vcf 
       	"""
}

// Variant to Table 

process VarToTable {
	conda "bioconda::gatk4=4.4"
	tag "Collect Variant in a Table using GATK4"
	publishDir "${params.outdir}/Mapping/Variants", mode: 'copy'
	
	input:
		path Rawvcf	// vcf files from RawHaploCall 
		path Recalvcf	// vcf files from RecalHaploCAll
		
	output:
		path "${Rawvcf}.table"
		path "${Recalvcf}.table"	
	script:
	"""
	gatk VariantsToTable \
			--fields CHROM -F POS -F TYPE -GF GT \
			--variant ${Rawvcf} \
			--output ${Rawvcf}.table
	gatk VariantsToTable \
			--fields CHROM -F POS -F TYPE -GF GT \
			--variant ${Recalvcf} \
			--output ${Recalvcf}.table
	"""
}

// SNP Filtering from Gatk vcf outputs. 

process SnpFilter {
conda "bioconda::gatk4=4.4"
	tag "Collect Variant in a Table using GATK4"
	publishDir "${params.outdir}/Mapping/Variants", mode: 'copy'
	
	input:
		path variants
		
	output:
		path "${variants.baseName}.SNP.vcf"
	script:
	"""
	gatk SelectVariants \
			--variant ${variants} \
			--select-type-to-include SNP \
			--output ${variants.baseName}.SNP.vcf
	"""
}
// channels 

// Reference file channel 

Channel.fromPath(params.Reference)
       .first()
       .set{ ref_file }
       
// Reads channel

Channel.fromPath(params.reads, checkIfExists: true)  
       .splitCsv(header: true)  
       .map { row -> tuple(row.patient_id, file(row.R1), file(row.R2)) } 
       .set{ READS }
       
// knwon file 1 channel for BQSR    
   
Channel.fromPath(params.knwonSite1, checkIfExists: true)  
       .first()
       .set{knwonSite1}    
       
// knwon file 2 channel for BQSR       
         
Channel.fromPath(params.knwonSite2, checkIfExists: true)  
       .first()
       .set{knwonSite2}        
    
// Main workflow

 workflow {
 	IndexRef	(	ref_file								)

	alignReadsToRef	(	ref_file,
      		      		IndexRef.out.idxREF.collect(),
      		      		READS									)	
      		      
      	assignReadGroup	(	alignReadsToRef.out.collectFile( sort: true),
                      		READS									)
                      
	markDuplicates	(	assignReadGroup.out.collectFile( sort: true),
           	     		READS									)
           	     		
     	IndexBam	(	markDuplicates.out.sorted_markduplicates_bam.collectFile(sort: true)	)   	
         	     	
      	GenerateStat	(	assignReadGroup.out.sorted_labeled_bam.collectFile(sort: true),
              	   		markDuplicates.out.sorted_markduplicates_bam.collectFile(sort: true)	)
              	   	
        IndexKNownSites	(	knwonSite1,knwonSite2							)     	
	
	BaseRecalibrator(	ref_file,
				IndexRef.out.dicREF,
				IndexRef.out.samidxREF,
				markDuplicates.out.sorted_markduplicates_bam.collectFile(sort: true),
				knwonSite1,
				IndexKNownSites.out.IDXknS1,
				knwonSite2,
				IndexKNownSites.out.IDXknS2					)
	
	ApplyBQSR	(	markDuplicates.out.sorted_markduplicates_bam.collectFile(sort: true),
				BaseRecalibrator.out.BQSR_Table.collectFile(sort:true)			)	
						

 	IndexRecalBam	(	ApplyBQSR.out.recal_bam.collectFile(sort: true)				)
 	
 	
	RawHaploCall	(	ref_file,
				IndexRef.out.dicREF,
				IndexRef.out.samidxREF,
				markDuplicates.out.sorted_markduplicates_bam.collectFile(sort: true),
				IndexBam.out.IDXBAM.collectFile(sort: true)				)
				
	RecalHaploCall	(	ref_file,
				IndexRef.out.dicREF,
				IndexRef.out.samidxREF,
				ApplyBQSR.out.recal_bam.collectFile(sort: true),
				IndexRecalBam.out.IDXRECALBAM.collectFile(sort: true)			)
   
   	VarToTable 	(	RawHaploCall.out.vcf_HaplotypeCaller_Raw.collectFile(sort: true),
   				RecalHaploCall.out.vcf_HaplotypeCaller_Recal.collectFile(sort: true)	)
   	SnpFilter 	(	RecalHaploCall.out.vcf_HaplotypeCaller_Recal.collectFile(sort: true)	)
}       






