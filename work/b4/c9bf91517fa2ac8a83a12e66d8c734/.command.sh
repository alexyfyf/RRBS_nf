#!/bin/bash -ue
fastqc --quiet --threads 8 sample_R1.fq.gz sample_R2.fq.gz
