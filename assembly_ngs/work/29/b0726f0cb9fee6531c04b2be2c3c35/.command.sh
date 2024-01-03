#!/bin/bash -ue
bwa mem  /media/houcem/theDrum/BILIM/github/FHG_clinical_bioinformatics/assembly_ngs/ref_genome/dpyd.fa sample1_r1.fq.gz sample1_r2.fq.gz  | samtools view -Sb - > patient1_unsorted.bam
