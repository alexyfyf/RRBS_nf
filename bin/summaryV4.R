#!/usr/bin/env Rscript

## args[1] is the samplesheet
## args[2] is the genome version

args <- commandArgs(trailingOnly=TRUE)

# args <- c("samplesheet.csv","mm10")
## lns <- length(args)

print(args)
## print(lns)

if (!require(tidyverse)) {
  install.packages(tidyverse)
}
if (!require(ggridges)) {
  install.packages(ggridges)
}
if (!require(methylKit)) {
  BiocManager::install(methylKit)
}
if (!require(ChIPseeker)) {
  BiocManager::install(ChIPseeker)
}
if (!require(ChIPseeker)) {
  BiocManager::install(ChIPseeker)
}
if (!require(annotatr)) {
  BiocManager::install(annotatr)
}
if (!require(R.utils)) {
  install.packages(R.utils)
}


# suppressPackageStartupMessages({
#   library(tidyverse)
#   library(methylKit)
#   library(R.utils)
#   library(ChIPseeker)
#   library(annotatr)
#   library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#   library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#   library(org.Mm.eg.db)
#   library(org.Hs.eg.db)
#   library(ggridges)
# })

samplesheet <- read.csv(args[1], header = FALSE)
group <- samplesheet$V2
cond <- levels(group)

# samplesheet

for (i in cond) {
  dir.create(i)
  filebasename <- paste0(samplesheet$V1[samplesheet$V2==i]) %>% as.character() %>%
    gsub(pattern = "\\.(fq|fastq)\\.gz", replacement = "") ## support both fq.gz and fastq.gz files

  lapply(list.files(pattern = paste0(filebasename, collapse = "|")), function(x){
    file.symlink(file.path("..",x), file.path(i))
  })

}

dirs <- file.path(cond) %>% paste0("/")
lns <- length(cond)

## source("readBismarkFiles.R")

## QC1: Read numbers
## change1: use bismark report instead of multiqc report, for better portablibity
## change2: change group to basename of the folder, for generalizability

qclist <- lapply(1:lns, function(x){
  path <- file.path(dirs[x])
  files <- list.files(path, pattern = "[SP]E_report.txt")
  readnum <- sapply(files, function(f) grep("unique best hit", readLines(paste0(path,f)), value = T))
  readnum <- as.numeric((readnum %>% str_split("\t",simplify = T))[,2])
  qc <- data.frame(readnum=readnum, group= basename(dirs[x]))
  qc })

# saveRDS(qclist,"qclist.rds")


Reduce(rbind, qclist) %>% mutate(group=factor(group)) %>%
  ggplot(aes(x = group, y = readnum / 1e6,
             fill = group)) +
  geom_boxplot() + ylim(0, NA) +
  stat_summary(fun=mean, colour="darkred", geom="text", show.legend = FALSE,
               vjust=-0.7, aes( label=round(..y.., digits=1)))

print("Save plot for aligned read numbers")
ggsave("qc_aligned_readnum.png")


## CpGs covered

cpgslist <- lapply(1:lns, function(x){
  path <- list.files(dirs[x], pattern = "*.cov.gz")
  cpgs <- sapply(paste0(dirs[x], path), countLines)
  cpgs <- data.frame(cpg_1x=cpgs,
                     group=basename(dirs[x]))
})
# saveRDS(cpgslist, "cpgslist.rds")

Reduce(rbind, cpgslist) %>% mutate(group=factor(group)) %>%
  ggplot(aes(x=group, y=cpg_1x/1e6, fill=group))+
  geom_boxplot() + ylim(0, NA) +
  stat_summary(fun=mean, colour="darkred", geom="text", show.legend = FALSE,
               vjust=-0.7, aes( label=round(..y.., digits=1)))

print("Save plot for number of CpGs covered at 1x")
ggsave("CpG_1x.png")


## CpGs at 10x
covlist <- lapply(1:lns, function(x){
  path <- list.files(dirs[x], pattern = "*.cov.gz")

  cov <- methRead(paste0(dirs[x], path) %>% as.list(),
                  assembly="mm10",mincov = 10,
                  pipeline = "bismarkCoverage",
                  sample.id=as.character(1:length(path)) %>% as.list(),
                  treatment=rep(1, length(path)))
})
# saveRDS(covlist,"covlist.rds")

cpg10xlist <- lapply(1:lns, function(x){
  cpg10x <- data.frame(cpg_10x=sapply(covlist[[x]], nrow),
                       group=basename(dirs[x]))
})

Reduce(rbind, cpg10xlist) %>% mutate(group=factor(group)) %>%
  ggplot(aes(x=group, y=cpg_10x/1e6, fill=group))+
  geom_boxplot() + ylim(0, NA) +
  stat_summary(fun=mean, colour="darkred", geom="text", show.legend = FALSE,
               vjust=-0.7, aes( label=round(..y.., digits=1)))

print("Save plot for number of CpGs covered at 10x")
ggsave("CpG_10x.png")

## ridges plot
histlist <- lapply(1:lns, function(x){
  cov <- lapply(covlist[[x]],function(m) {
    getData(m) %>% mutate(sampleid=m@sample.id)}) %>%
    Reduce(rbind,.) %>%
    mutate(group=basename(dirs[x]))
})

Reduce(rbind, histlist) %>% dplyr::select(coverage, sampleid,group) %>%
  ggplot(aes(x=log10(coverage),y=as.factor(group),
             fill=factor(stat(quantile))))+
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE,
    bandwidth = 0.025, from = 1
  ) +
  scale_fill_viridis_d(name = "Quartiles")

print("Save plot: CpG coverage ridge plot")
ggsave("ridge_plot.png")

## regions of CpGs

if (args[2]=="mm10") {
  if (!require(TxDb.Mmusculus.UCSC.mm10.knownGene)) {
    BiocManager::install(TxDb.Mmusculus.UCSC.mm10.knownGene)
  }
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  # orgdb <- "org.Mm.eg.db"
} else if (args2=="hg38") {
  if (!require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
    BiocManager::install(TxDb.Hsapiens.UCSC.hg38.knownGene)
  }
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  # orgdb <- "org.Hs.eg.db"
} else print("not supported and stop")

## CpGs had more than 10x coverage in all samples
## currently need 2 samples in each group at least
methlist <- lapply(1:lns, function(x){
  meth <- methylKit::unite(covlist[[x]],
                             # min.per.group = 2L,
                             destrand=FALSE)
  # meth <- meth[data.frame(meth)$chr %in% c(1:19,"X","Y","MT"),]
})

grlist <- lapply(1:lns, function(x) {
  gr <- as(methlist[[x]],"GRanges")
  seqlevelsStyle(gr) <- "UCSC"
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  mcols(gr) <- NULL
  gr
})

annolist <- lapply(grlist, function(x){
  annotatePeak(x, tssRegion=c(-3000, 3000),
               # annoDb=orgdb,
               TxDb=txdb)
})

annostat <- lapply(1:lns, function(x) {
  num <- annolist[[x]]@annoStat$Frequency/100 * (annolist[[x]]@anno %>% data.frame() %>% nrow())
  region <- annolist[[x]]@annoStat$Feature %>% as.character()
  df <- data.frame(region=region, num=num,
                   group=basename(dirs[x]))
})

annostat %>% Reduce(rbind,.) %>% mutate(group=factor(group)) %>%
  ggplot(aes(x = group, y = num/1e6, fill = region)) +
  geom_bar(stat = "identity", position = "stack")


print("Save plot annotation distribution using ChIPseeker")
ggsave("AnnotationStat_stack.png")

annostat %>% Reduce(rbind,.) %>% mutate(group=factor(group)) %>%
  ggplot(aes(x = group, y = num/1e6, fill = region)) + geom_bar(stat = "identity", position = "fill")

ggsave("AnnotationStat_pct.png")

## number of genes covered (count genes with any promoter coverage)

gene1k <- sapply(1:lns, function(x){
  annolist[[x]]@anno %>% data.frame() %>%
    filter(annotation %in% "Promoter (<=1kb)") %>%
    pull(geneId) %>% unique() %>% length()
})

gene2k <- sapply(1:lns, function(x){
  annolist[[x]]@anno %>% data.frame() %>%
    filter(annotation %in% c("Promoter (<=1kb)","Promoter (1-2kb)")) %>%
    pull(geneId) %>% unique() %>% length()
})

gene3k <- sapply(1:lns, function(x){
  annolist[[x]]@anno %>% data.frame() %>%
    filter(annotation %in% c("Promoter (<=1kb)","Promoter (1-2kb)","Promoter (2-3kb)")) %>%
    pull(geneId) %>% unique() %>% length()
})

cbind(p1k=gene1k, p2k=gene2k, p3k=gene3k) %>% data.frame() %>%
  mutate(group=cond) %>%
  gather(distance, num_gene, 1:3) %>%
  ggplot(aes(x=group,y=num_gene, group=distance, col=distance)) +
  geom_point() +geom_line() +
  stat_summary(fun=mean, colour="darkred", geom="text", show.legend = FALSE,
               vjust=-0.7, aes( label=round(..y.., digits=1)))

print("Save plot numbers of gene promoters covered using different TSS")
ggsave("gene_covered_promoter.png")

## another way to annotate

cpgs_info <- build_annotations(genome = args[2], annotations = paste(args[2],"cpgs",sep = "_"))

cpgannolist <- lapply(1:lns, function(x){
  anno <- annotate_regions(regions = grlist[[x]],
                           annotations = cpgs_info,
                           ignore.strand = TRUE,
                           quiet = FALSE)
})

cpgannostat <- lapply(1:lns, function(x){
  num <- cpgannolist[[x]] %>% summarize_annotations() %>% pull(n)
  region <- c("Inter","Islands","Shelves","Shores")
  df <- data.frame(region=region, num=num,
                   group=basename(dirs[x]))
})

cpgannostat %>% Reduce(rbind,.) %>% mutate(group=factor(group)) %>%
  ggplot(aes(x = group, y = num/1e6, fill = region)) + geom_bar(stat = "identity", position = "stack")

print("Save plot CpG island annotation distribution using")
ggsave("CpGStat_stack.png")

cpgannostat %>% Reduce(rbind,.) %>% mutate(group=factor(group)) %>%
  ggplot(aes(x = group, y = num/1e6, fill = region)) + geom_bar(stat = "identity", position = "fill")

ggsave("CpGStat_pct.png")

save.image(file = "alldata.RData")
