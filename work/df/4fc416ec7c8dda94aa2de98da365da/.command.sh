#!/bin/bash -ue
mkdir bismarkindex
cp data/*.fa* bismarkindex/
bismark_genome_preparation --hisat2 --parallel 2--verbose bismarkindex
