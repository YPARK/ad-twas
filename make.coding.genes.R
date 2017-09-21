#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) {
    q()
}

options(stringsAsFactors = FALSE)
library(dplyr)
library(rtracklayer)

## select coding genes in GTF
gtf.file <- argv[1]       # e.g., gtf.file <- 'gencode.v19.genes.patched_contigs.gtf.gz'
data.gene.file <- argv[2] # e.g., data.gene.file <- 'data/rnaseq.rows.txt.gz'
out.file <- argv[3]

gtf.tab <- readGFF(gtf.file, tags = c('gene_id', 'gene_name', 'transcript_name', 'gene_type'))

strip.ensg <- function(x) strsplit(x, split = '[.]')[[1]][1]

coding.genes <-
    gtf.tab %>% mutate(chr = seqid, ensg = gene_id) %>%
        filter(gene_type == 'protein_coding', type == 'transcript') %>%
            select(chr, start, end, strand, ensg, gene_name) %>%
                arrange(chr, start) %>%
                    mutate(ensg = sapply(ensg, strip.ensg))

data.genes <- read.table(data.gene.file, col.names = c('ensg', 'data.loc')) %>%
    mutate(ensg = sapply(ensg, strip.ensg))

out.tab <- coding.genes %>%
    left_join(data.genes, by = 'ensg') %>%
        na.omit() %>%
            filter(chr %in% 1:22)

write.table(out.tab, file = gzfile(out.file), sep = '\t',
            quote = FALSE, row.names = FALSE, col.names = FALSE)
