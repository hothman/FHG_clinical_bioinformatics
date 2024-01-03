#!/bin/bash -ue
if [ ! -e dpyd.fa.bwt ]; then
    bwa index dpyd.fa
else
    echo "Index files for dpyd.fa already exist."
fi
