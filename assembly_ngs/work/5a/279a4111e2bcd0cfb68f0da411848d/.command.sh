#!/bin/bash -ue
bwa-mem2 mem  /media/houcem/theDrum/BILIM/github/FHG_clinical_bioinformatics/assembly_ngs/ref_genome/dpyd.fa sample1_r1.fq.gz sample1_r2.fq.gz -t 2 \
                     | samtools addreplacerg -r ID:sampl1 -r PL:Illumina  \
                     | samtools view -Sb -@ 2 | samtools sort -o patient1_sorted.bam
