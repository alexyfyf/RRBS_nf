#!/bin/bash -ue
trim_galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC \
--paired sample_R1.fq.gz sample_R2.fq.gz --cores 8

python2 /home/ubuntu/software/NuMetRRBS/trimRRBSdiversityAdaptCustomers.py -1 sample_R1_val_1.fq.gz \
-2 sample_R2_val_2.fq.gz &> sample_trimpy.log
