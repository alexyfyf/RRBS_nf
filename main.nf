/* 
 * 'RRBS_nf' - A Nextflow pipeline for RRBS data analysis
 * 
 * This pipeline deals with RRBS data from NuGen Ovation library 
 * input is BCL files or Run folder, or reads
 * 
 * Feng Yan
 * feng.yan@monash.edu
 */


/*
 * Define the default parameters
 */ 
 
params.genomedir     = "$baseDir/data/"
params.numet      = "~/software/NuMetRRBS/trimRRBSdiversityAdaptCustomers.py"
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.results    = "results"

log.info """\
R R B S -  N F    v 1.0 
================================
genomedir: $params.genomedir
numet    : $params.numet
reads    : $params.reads
results  : $params.results
"""

/*
 *  Parse the input parameters
 */

reads_ch        = Channel.fromFilePairs(params.reads)


/**********
 * PART 1: Data preparation
 *
 * Process 1A: Create a bisulfite converted genome index (.fai) with bismark
 */

process '1A_prepare_bisulfite_genome' { 
  tag "$genome.baseName"
 
  input: 
      path genome from $params.genomedir
  
  script:
  """
  bismark_genome_preparation --hisat2 --verbose $genome
  """
}

 
 
