#!/bin/bash -ue
gatk MarkDuplicates --INPUT patient1_sorted.bam \
                    --OUTPUT patient1_sorted_markduplicates.bam
                    --METRICS_FILE metrict.txt
                    --TMP_DIR .
