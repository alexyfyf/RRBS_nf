/*
 * 'RRBS_nf' - A Nextflow pipeline for RRBS data analysis
 *
 * This pipeline deals with RRBS data from NuGen Ovation library
 * input is reads in FASTQ format
 *
 * Feng Yan
 * feng.yan@monash.edu
 */


/*
 * Define the default parameters
 */

nextflow.enable.dsl=2

params.genomedir  = "$baseDir/data/ref/"
params.numet      = "$baseDir/bin/trimRRBSdiversityAdaptCustomers.py"
params.reads      = "$baseDir/data/sample_R{1,2}.fq.gz"
params.outdir     = "results"
params.aligner    = "bismark_hisat"
params.summary    = "$baseDir/bin/summaryV4.R"
params.species    = "mm10"
params.samplesheet= "$baseDir/data/samplesheet.csv"
params.library    = "nugen"
params.clip3prime = 5
params.trim       = "trim"
params.single_end = true

// Convert single_end to boolean
// params.single_end = params.single_end.toBoolean()

log.info """\
R R B S -  N F    v 2.0
================================
genomedir	: $params.genomedir
species  	: $params.species
reads    	: $params.reads
outdir   	: $params.outdir
samplesheet	: $params.samplesheet
library  	: $params.library
trim        : $params.trim
single_end  : $params.single_end
"""

/*
 * Create a channel for input files
 */
ch_reads = Channel
    .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!" }

ch_reads_meta = Channel
    .fromPath(params.reads)
    .map { id, reads ->
        def meta = [:]
        meta.id = id
        meta.single_end = params.single_end
        meta.library = params.library
        return [ meta, reads ]
    }
ch_genome = Channel.fromPath("${params.genomedir}/*.{fa,fasta,fna}")
ch_numet = Channel.fromPath(params.numet)
ch_samplesheet = Channel.fromPath(params.samplesheet)
ch_summary = Channel.fromPath(params.summary)
ch_species = Channel.value(params.species)

/*
 * PART 0: Preparation
 */
process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    executor 'local'

    output:
    file '*.txt'

    script:
    """
    echo "$workflow.manifest.version" &> v_ngi_methylseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    bismark --version &> v_bismark.txt
    bismark_genome_preparation --version &> v_bismark_genome_preparation.txt
    deduplicate_bismark --version &> v_deduplicate_bismark.txt
    bismark_methylation_extractor --version &> v_bismark_methylation_extractor.txt
    bismark2report --version &> v_bismark2report.txt
    bismark2summary --version &> v_bismark2summary.txt
    samtools --version &> v_samtools.txt
    hisat2 --version &> v_hisat2.txt
    bowtie2 --version &> v_bowtie2.txt
    multiqc --version &> v_multiqc.txt
    R --version &> v_R.txt
    """
}

/**********
 * PART 1: Preprocessing
 *
 * Process 1A: fastqc report for raw data
 */
process FASTQC {
    tag "$meta.id"
    label 'big'
    publishDir "${params.outdir}/pre_fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    tuple val(meta), path(reads)

    output:
    path "*_fastqc.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet --threads ${task.cpus} $reads
    """
}


/*
 * Process 1B: 2-step triming for NuGen RRBS data
 */
process TRIM {
    tag "$meta.id"
    label 'big'
    publishDir "${params.outdir}/trim"

    input:
    tuple val(meta), path(reads)
    path numet
    // val(clip3) from params.clip3prime

    output:
    tuple val(meta), path('*.fq_trimmed.fq.gz'), emit: trimmed_reads
    path '*report.txt', emit: trim_reports
    path '*.log', optional: true
     
    script:
    // def single_end = params.single_end
    // def library = params.library
    // def trim = params.trim

    if( params.library == "nugen" ) {
    if( params.single_end ) {
            """
            trim_galore -a AGATCGGAAGAGC $reads --cores ${task.cpus}
            echo $trimpy
            python2 $trimpy -1 ${reads.simpleName}_trimmed.fq.gz &> ${reads.simpleName}_trimpy.log
            """
        } else {
            """
            trim_galore -a AGATCGGAAGAGC -a2 AAATCAAAAAAAC \\
            --paired $reads --cores ${task.cpus}
            echo $trimpy
            python2 $trimpy -1 ${reads[0].simpleName}_val_1.fq.gz \\
            -2 ${reads[1].simpleName}_val_2.fq.gz &> ${name}_trimpy.log
            """
        }
    } else if (params.library == "epic" && params.trim == "trim") { // adde trim for Epic becasue we noticed bias toward 3 prime end
    if( params.single_end ) {
            // clip3 = params.clip3prime == 0 ? "" : "--three_prime_clip_R1 $params.clip3prime"
            """
            ## leave to auto detection
            trim_galore $reads --cores ${task.cpus}
            mv ${reads.simpleName}_trimmed.fq.gz ${reads.simpleName}.fq_trimmed.fq.gz
            """
        } else {
            // clip3 = params.clip3prime == 0 ? "" : "--three_prime_clip_R1 $params.clip3prime --three_prime_clip_R2 $params.clip3prime"
            """
            trim_galore --paired $reads --cores ${task.cpus}
            mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}.fq_trimmed.fq.gz 
            mv ${reads[1].simpleName}_val_2.fq.gz ${reads[1].simpleName}.fq_trimmed.fq.gz
            """
        }
    } else if (params.library == "epic" && params.trim == "skip") {
    if( params.single_end ) {
            """
            mv ${reads} ${reads.simpleName}.fq_trimmed.fq.gz
            """
        } else {
            """
            mv ${reads[0]} ${reads[0].simpleName}.fq_trimmed.fq.gz
            mv ${reads[1]} ${reads[1].simpleName}.fq_trimmed.fq.gz
            """
        }
    } else if (params.library == 'emseq' && params.trim == 'trim') {
       if( params.single_end ) {
            // --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10
            """
            ## leave to auto detection
            trim_galore --clip_R1 10 --three_prime_clip_R1 10 $reads --cores ${task.cpus}
            mv ${reads.simpleName}_trimmed.fq.gz ${reads.simpleName}.fq_trimmed.fq.gz
            """
        } else {
            // --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10
            """
            trim_galore --paired --clip_R1 10 --clip_R2 10 --three_prime_clip_R1 10 --three_prime_clip_R2 10 $reads --cores ${task.cpus}
            mv ${reads[0].simpleName}_val_1.fq.gz ${reads[0].simpleName}.fq_trimmed.fq.gz 
            mv ${reads[1].simpleName}_val_2.fq.gz ${reads[1].simpleName}.fq_trimmed.fq.gz
            """
        } 
    }
}

/**********
 * Process 1C: fastqc report for trimmed data
 */
process '1C_post_fastqc' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/post_fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from clean_reads_fastqc_ch

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc2_results_for_multiqc
     
    when:
    params.trim == "trim"
 
    script:
    """
    fastqc --quiet --threads ${task.cpus} $reads
    """
}


/**********
 * PART 2: Bismark
 *
 * Process 2A: Create a bisulfite converted genome index (.fai) with bismark
 */

process '2A_prepare_bisulfite_genome' {
  tag "$genome.baseName"
  label 'bismark'

  input:
      file genome from genomedir
  output:
      file 'bismarkindex' into genome_dir_ch

  script:
  aligner = params.aligner == 'bismark_hisat' ? '--hisat2' : '--bowtie2'
  """
  mkdir bismarkindex
  cp ${genome}/*.fa* bismarkindex/
  bismark_genome_preparation $aligner --parallel ${task.cpus / 4} --verbose bismarkindex
  """
}



/**********
 * Process 2B: Align RRBS reads to the genome
 */

process '2B_mapping_bismark' {
  tag "$name"
  label 'bismark'
  publishDir "${params.outdir}/bismark_align", mode: 'copy'

  input:
      file genomeDir from genome_dir_ch
      set val(name), file(reads) from clean_reads_bismark_ch

  output:
      set val(name), file('*.bam') into aligned_bam_ch, ch_bam_for_bismark_summary
      set val(name), file('*report.txt') into ch_bismark_align_log_for_multiqc, ch_bismark_align_log_for_bismark_report, ch_bismark_align_log_for_bismark_summary, ch_bismark_align_log_for_Rsummary

  script:
  // Paired-end or single end input files
  input = params.single_end ? reads : "-1 ${reads[0]} -2 ${reads[1]}"

  // Choice of read aligner
  aligner = params.aligner == "bismark_hisat" ? "--hisat2" : "--bowtie2"
  hisat2 = params.aligner == "bismark_hisat" ? "--no-spliced-alignment" : ""
  """
  # echo $input
  bismark --genome $genomeDir \\
       $aligner $hisat2 \\
       --multicore ${task.cpus / 4} \\
       $input
  ## can use sorted bam to reduce size, but for paired end data, need to use read name sorting, which doesn't really reduce size
  # samtools sort -n -o ${name}_sorted.bam -@ ${task.cpus} ${name}*bismark*.bam
  # rm ${name}*bismark*.bam 
  """
}



/***********
 * Process 2C: - Bismark methylation extraction
 */
process '2C_bismark_methXtract' {
    tag "$name"
    label 'bismark'
    publishDir "${params.outdir}/bismark_methylation", mode: 'copy'

    input:
    set val(name), file(bam) from aligned_bam_ch

    output:
    set val(name), file("*splitting_report.txt") into ch_bismark_splitting_report_for_bismark_report, ch_bismark_splitting_report_for_multiqc, ch_bismark_splitting_report_for_bismark_summary
    set val(name), file("*.M-bias.txt") into ch_bismark_mbias_for_bismark_report, ch_bismark_mbias_for_multiqc, ch_bismark_mbias_for_bismark_summary
    set val(name), file("*.cov.gz") into coverage_bismark_ch, covgz_for_Rsummary
    set val(name), file("*.bedGraph.gz") into bedgraph_bismark_ch
    // file '*.{png,gz}'
    file '*.bedGraph.gz'
    file '*.cov.gz'
    file '*.CpG_report.txt.gz'
    file '*.txt'

    script:
    // cytosine_report = params.cytosine_report ? "--cytosine_report --genome_folder ${index} " : ''

    """
    bismark_methylation_extractor --comprehensive --merge_non_CpG \\
    --multicore ${task.cpus / 4} \\
    --cytosine_report --genome_folder $genomedir \\
    --ample_memory \\
    --no_overlap \\
    --bedGraph \\
    --gzip \\
    --report \\
    $bam
    """
}

ch_bismark_align_log_for_bismark_report
  .join(ch_bismark_splitting_report_for_bismark_report)
  .join(ch_bismark_mbias_for_bismark_report)
  .set{ ch_bismark_logs_for_bismark_report }

// bedgraph_bismark_ch.view()

/*
 * Process 2D: Bismark Sample Report
 */
process '2D_bismark_report' {
    tag "$name"
    publishDir "${params.outdir}/bismark_reports", mode: 'copy'

    input:
    set val(name), file(align_log), file(splitting_report), file(mbias) from ch_bismark_logs_for_bismark_report

    output:
    file '*{html,txt}'

    script:
    """
    bismark2report \\
        --alignment_report $align_log \\
        --splitting_report $splitting_report \\
        --mbias_report $mbias
    """
}


/*
 * Process 2E - Bismark Summary Report
 */
process '2E_bismark_summary' {
    publishDir "${params.outdir}/bismark_summary", mode: 'copy'

    input:
    file ('*') from ch_bam_for_bismark_summary.collect()
    file ('*') from ch_bismark_align_log_for_bismark_summary.collect()
    file ('*') from ch_bismark_splitting_report_for_bismark_summary.collect()
    file ('*') from ch_bismark_mbias_for_bismark_summary.collect()

    output:
    file '*{html,txt}'

    script:
    """
    bismark2summary
    """
}

/**********
 * PART 3: Summary
 *
 * Process 3A: MultiQC
 */
process '3A_multiqc' {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    // only use fastqc from trimmed reads
    // file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
    file ('post_fastqc/*') from ch_fastqc2_results_for_multiqc.collect().ifEmpty([])
    file ('trim/*') from ch_trimgalore_results_for_multiqc.collect().ifEmpty([])
    file ('bismark_align/*') from ch_bismark_align_log_for_multiqc.collect().ifEmpty([])
    file ('bismark_methylation/*') from ch_bismark_splitting_report_for_multiqc.collect().ifEmpty([])
    file ('bismark_methylation/*') from ch_bismark_mbias_for_multiqc.collect().ifEmpty([])

    output:
    file "*multiqc_report.html"
    file "*_data"

    script:
    """
    multiqc -f .
    """
}

/**********
 * PART 4: Visualization
 *
 * Process 4A: Generate fasta index
 */
process '4A_faidx' {
    tag "$fasta.baseName"
    label 'big'

    input:
    file fasta from fasta_ch

    output:
    file "chrom.sizes" into chr_size_ch

    script:
    if( fasta.extension ==~ /fa|fasta/ ) {
            """
            samtools faidx ${fasta}
            cut -f1,2 ${fasta}.fai | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' > chrom.sizes
            """
       } else if( fasta.extension == 'gz' ) {
	    """
            zcat ${fasta} | bgzip -c > ${fasta.simpleName}.fa.bgz
            samtools faidx ${fasta.simpleName}.fa.bgz
            cut -f1,2 ${fasta.simpleName}.fa.bgz.fai | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' > chrom.sizes
            """
       }
}

// chr_size_ch.view()
// bedgraph_bismark_ch.view()
// bedgraph_bismark_ch.combine(chr_size_ch).view()

/**********
 * Process 4B: Generate bigwig files
 */
process '4B_toBigWig' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    //set val(name), file(bedgraph) from bedgraph_bismark_ch
    //file chrsize from chr_size_ch
    set val(name), file(bedgraph), file(chrsize) from bedgraph_bismark_ch.combine(chr_size_ch)

    output:
    file "*.bw"

    script:
    """
    zcat $bedgraph | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' | sort -k1,1 -k2,2n > ${bedgraph.simpleName}.bedGraph
    bedGraphToBigWig ${bedgraph.simpleName}.bedGraph ${chrsize} ${bedgraph.simpleName}.bw
    """
}

//covgz_for_Rsummary
//  .join(ch_bismark_align_log_for_Rsummary).collect().view()

/**********
 * Process 4C: Generate summary statistics
 */
process '4c_toRSummary' {
    tag "summaryplot"
    label 'bismark'
    publishDir "${params.outdir}/summaryplot", mode: 'copy'

    input:
    file("*") from covgz_for_Rsummary.join(ch_bismark_align_log_for_Rsummary).collect()
    file(samplesheet) from samplesheet
    val(species) from species
    file(summary) from summary
    
    output:
    file "*.png"
    // file "*.RData"

    script:
    """
    # module load R
    Rscript --vanilla $summary $samplesheet $species
    """
}


workflow {
    // Main workflow
    GET_SOFTWARE_VERSIONS()
    
    // Create metadata structure
    ch_reads
        .map { id, reads ->
            def meta = [:]
            meta.id = id
            meta.single_end = params.single_end
            return [ meta, reads ]
        }
        .set { ch_reads_meta }

    // Run FastQC on raw reads
    FASTQC(ch_reads_meta)

    // Run trimming
    TRIM(ch_reads_meta, ch_numet)

    // Run post-trim FastQC
    FASTQC(TRIM.out.trimmed_reads)

    // Continue with other processes...
}