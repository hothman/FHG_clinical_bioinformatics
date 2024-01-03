#!/bin/bash -ue
mkdir patient1_fastqc
fastqc sample1_r1.fq.gz sample1_r2.fq.gz -o patient1_fastqc
