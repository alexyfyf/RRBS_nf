#!/bin/bash -ue
bismark2report \
    --alignment_report sample_R1_val_1.fq_trimmed_bismark_hisat2_PE_report.txt \
    --splitting_report sample_R1_val_1.fq_trimmed_bismark_hisat2_pe_splitting_report.txt \
    --mbias_report sample_R1_val_1.fq_trimmed_bismark_hisat2_pe.M-bias.txt
