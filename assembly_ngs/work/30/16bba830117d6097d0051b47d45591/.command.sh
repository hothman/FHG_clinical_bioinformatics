#!/bin/bash -ue
gatk MarkDuplicates --INPUT patient1_sorted.bam \
                    --OUTPUT patient1_sorted_markduplicates.sam \
                    --METRICS_FILE patient1.metrict \
                    --TMP_DIR .
