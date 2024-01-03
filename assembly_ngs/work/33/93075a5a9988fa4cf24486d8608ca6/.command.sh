#!/bin/bash -ue
# Check if the index file (.fai) exists in the output directory
if [ ! -e ./ref_genome/dpyd.fa ]; then
echo index file not found 
    # Generate the index file using samtools
    # Replace this with the appropriate command for your indexing tool
    samtools faidx dpyd -o ./ref_genome/dpyd.fa/dpyd.fai
fi
