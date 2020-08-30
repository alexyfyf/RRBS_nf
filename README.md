# RRBS_nf
Nextflow pipeline for analyzing RRBS data from NuGen Ovation library

## Inputs
This pipeline works with FASTQ files
* If library is single-end, please use `--single_end --reads '*.fq.gz'`.
* If library is paired-end, please use `--reads '*{1,2}.fq.gz'`.
* Please use quotes for the input file names with wildcards.


## Bismark parameters
Support `hisat2` (default) and `bowtie2`

## Quick start
```nextflow main.nf --reads 'data/*R{1,2}.fq.gz'  --genomedir data/ref/ -profile slurm```

```nextflow main.nf --single_end --reads 'data/*.fq.gz'  --genomedir data/ref/ -profile slurm```

`-profile slurm` enbales the use of powerful resourse management software SLURM, and the relavant setting for specific task can be changed in `nextflow.config` file.

## Updates:
* Adapted from nf-core methyl-seq pipeline (https://github.com/nf-core/methylseq)

* Simplified and added NuGen specific trimming process (https://github.com/nugentechnologies/NuMetRRBS)

