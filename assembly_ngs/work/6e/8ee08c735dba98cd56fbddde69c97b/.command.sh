#!/bin/bash -ue
bwa-mem2 mem -R '@RG	ID:group1	SM:sample1	PL:illumina	LB:lib1	PU:unit1' \
                     /media/houcem/theDrum/BILIM/github/FHG_clinical_bioinformatics/assembly_ngs/ref_genome/dpyd.fa sample1_r1.fq.gz sample1_r2.fq.gz -t 2 \
                     | samtools view -Sb -@ 2 | samtools sort -o patient1_sorted.bam
