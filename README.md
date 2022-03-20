# RRBS_nf
[![DOI](https://zenodo.org/badge/243173536.svg)](https://zenodo.org/badge/latestdoi/243173536)

Nextflow pipeline for analyzing RRBS data from NuGen Ovation library

## Analysis steps
1. FastQC before trimming (FastQC)
2. Trimming (trim_galore)
    * For NuGen Ovation kit, use adaptive trimming from manufacturer.
    * For Illumina Epic Truseq kit, use trim_galore auto trimming.
3. FastQC after trimming (FastQC)
4. Bismark alignment (Bismark)
5. Bismark methylation extraction (Bismark)
6. Generate bigwig files for visualization (ucsc utils)
7. Generate summary statistics based on group conditions 

## Command line options
- --genomedir     Directory for genome file, can be fa, fasta, fa.gz, fasta.gz.
- --library       Library kit used, can be `nugen` or `epic`. Default is `nugen`. This will affect how trimming is performed.
- --single_end    If library is single-end, use this flag. Default is not on.
- --reads         FASTQ files in fq, fastq, fq.gz or fastq.gz format. If library is single-end, please use `--single_end --reads '*.fq.gz'`. If library is paired-end, please use `--reads '*{1,2}.fq.gz'`.
- --outdir        Output directory. Default is `results`.
- --aligner       Aligner used by Bismark, can be `hisat2` or `bowtie2`. Default is `hisat2`.
- --species       Species to use in summary statistics, can be `mm10` or `hg38`. 
- --samplesheet   Sample sheet csv file contains sample grouping information. No heading, first column is fastq file name (for paired-end data, only use R1 file name), second column is group. Currently require > 2 samples in each group.

Optional augument related to nextflow
- -resume         Resume previous analysis.
- -profile        Profile to use. Profiles can be edited in `nextflow.config` file. Use `-profile slurm` to use slurm HPC scheduler.


## Inputs
This pipeline works with FASTQ files (fastq, fq, fastq.gz, fq.gz)
* If library is single-end, please use `--single_end --reads '*.fq.gz'`.
* If library is paired-end, please use `--reads '*{1,2}.fq.gz'`.
* Please use quotes for the input file names with wildcards.

## Bismark parameters
Support `hisat2` (default) and `bowtie2`

## Quick start
```nextflow main.nf --reads 'data/*R{1,2}.fq.gz'  --genomedir data/ref/ -profile slurm --samplesheet samplesheet.csv --species mm10```

```nextflow main.nf --single_end --reads 'data/*.fq.gz'  --genomedir data/ref/ -profile slurm --samplesheet samplesheet.csv --species mm10```

`-profile slurm` enbales the use of powerful resourse management software SLURM, and the relavant setting for specific task can be changed in `nextflow.config` file.

## Updates:
* Adapted from nf-core methyl-seq pipeline (https://github.com/nf-core/methylseq)

* Simplified and added NuGen specific trimming process (https://github.com/nugentechnologies/NuMetRRBS)

