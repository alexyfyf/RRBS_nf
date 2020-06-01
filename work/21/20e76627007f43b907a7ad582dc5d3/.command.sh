#!/bin/bash -ue
bismark_methylation_extractor --comprehensive --merge_non_CpG \
--multicore 8 \
--cytosine_report --genome_folder /Myvolume/fyan0011/GitHub/RRBS_nf/data \
--ample_memory \
--no_overlap \
--bedGraph \
--gzip \
--report \
sample_R1_val_1.fq_trimmed_bismark_hisat2_pe.bam
