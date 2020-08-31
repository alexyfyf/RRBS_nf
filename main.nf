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

params.genomedir  = "$baseDir/data/ref/"
params.numet      = "$baseDir/bin/trimRRBSdiversityAdaptCustomers.py"
params.reads      = "$baseDir/data/sample_R{1,2}.fq.gz"
params.outdir     = "results"
params.aligner    = "bismark_hisat"
params.summary    = "$baseDir/bin/summaryV4.R"
params.species    = "mm10"
params.samplesheet= "$baseDir/data/samplesheet.csv"


log.info """\
R R B S -  N F    v 1.0
================================
genomedir: $params.genomedir
species  : $params.species
reads    : $params.reads
outdir   : $params.outdir
sample   : $params.samplesheet
"""

/*
 *  Parse the input parameters
 */

genomedir       = file(params.genomedir)
numet           = file(params.numet)
Channel
        .fromPath( "${genomedir}/*.fa" )
        .set{ fasta_ch }
samplesheet     = file(params.samplesheet)
species         = Channel.from(params.species)
summary         = file(params.summary)

/*
 * PART 0: Preparation
 */
process '0A_get_software_versions' {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    executor 'local'

    output:
    file '*.txt'

    script:
    """
    echo "$workflow.manifest.version" &> v_ngi_methylseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    bismark_genome_preparation --version &> v_bismark_genome_preparation.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    bismark --version &> v_bismark.txt
    deduplicate_bismark --version &> v_deduplicate_bismark.txt
    bismark_methylation_extractor --version &> v_bismark_methylation_extractor.txt
    bismark2report --version &> v_bismark2report.txt
    bismark2summary --version &> v_bismark2summary.txt
    samtools --version &> v_samtools.txt
    hisat2 --version &> v_hisat2.txt
    # bwa &> v_bwa.txt 2>&1 || true
    # bwameth.py --version &> v_bwameth.txt
    # picard MarkDuplicates --version &> v_picard_markdups.txt 2>&1 || true
    # MethylDackel --version &> v_methyldackel.txt
    # qualimap --version &> v_qualimap.txt || true
    # preseq &> v_preseq.txt
    multiqc --version &> v_multiqc.txt
    # scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Create a channel for input read files
 */
Channel
        .fromFilePairs( params.reads, size: params.single_end ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { raw_reads_fastqc_ch; raw_reads_trim_ch }


/**********
 * PART 1: Preprocessing
 *
 * Process 1A: fastqc report for raw data
 */
process '1A_pre_fastqc' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from raw_reads_fastqc_ch

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc_results_for_multiqc

    script:
    """
    fastqc --quiet --threads ${task.cpus} $reads
    """
}


/*
 * Process 1B: 2-step triming for NuGen RRBS data
 */
process '1B_trim' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/trim"

    input:
    set val(name), file(reads) from raw_reads_trim_ch
    file(trimpy) from numet

    output:
    set val(name), file('*.fq_trimmed.fq.gz') into clean_reads_bismark_ch, clean_reads_fastqc_ch
    file('*report.txt') into ch_trimgalore_results_for_multiqc
    file('*.log')

    script:

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
}

/**********
 * Process 1C: fastqc report for trimmed data
 */
process '1C_post_fastqc' {
    tag "$name"
    label 'big'
    publishDir "${params.outdir}/fastqc2", mode: 'copy',
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }

    input:
    set val(name), file(reads) from clean_reads_fastqc_ch

    output:
    file '*_fastqc.{zip,html}' into ch_fastqc2_results_for_multiqc

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
  publishDir "${params.outdir}/bismark_align"

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
    file '*.{png,gz}'

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
    file ('fastqc/*') from ch_fastqc_results_for_multiqc.collect().ifEmpty([])
    file ('fastqc2/*') from ch_fastqc2_results_for_multiqc.collect().ifEmpty([])
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
    """
    samtools faidx ${fasta}
    cut -f1,2 ${fasta}.fai | sed -e 's/\\(^[0-9XY]\\)/chr\\1/' -e 's/^MT/chrM/' | grep '^chr' > chrom.sizes
    """
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
    label 'big'
    publishDir "${params.outdir}/summaryplot", mode: 'copy'

    input:
    file("*") from covgz_for_Rsummary.join(ch_bismark_align_log_for_Rsummary).collect()
    file(samplesheet) from samplesheet
    val(species) from species
    file(summary) from summary
    
    output:
    file "*.png"
    file "*.RData"

    script:
    """
    module load R
    Rscript --vanilla $summary $samplesheet $species
    """
}
