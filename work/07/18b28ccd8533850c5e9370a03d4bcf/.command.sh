#!/bin/bash -ue
# echo -1 sample_R1_val_1.fq_trimmed.fq.gz -2 sample_R2_val_2.fq_trimmed.fq.gz
bismark --genome bismarkindex \
     --hisat2 --no-spliced-alignment \
     --multicore 8 \
     -1 sample_R1_val_1.fq_trimmed.fq.gz -2 sample_R2_val_2.fq_trimmed.fq.gz
