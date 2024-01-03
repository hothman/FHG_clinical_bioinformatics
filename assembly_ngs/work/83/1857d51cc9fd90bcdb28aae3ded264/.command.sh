#!/bin/bash -ue
bwa-mem2 index dpyd.fa
gatk CreateSequenceDictionary REFERENCE=dpyd.fa \
                              OUTPUT=dpyd.fa.dict
