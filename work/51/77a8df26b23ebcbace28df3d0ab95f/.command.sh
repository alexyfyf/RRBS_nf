#!/bin/bash -ue
fastqc --quiet --threads 8 sample_R1_val_1.fq_trimmed.fq.gz sample_R2_val_2.fq_trimmed.fq.gz
