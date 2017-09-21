#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)
gene.file <- argv[1] # e.g., 'data/rnaseq/chr21-genes.txt.gz'

options(stringsAsFactors = FALSE)
library(dplyr)

block.size0 <- 1e6
max.size <- 20

gene.cols <- c('chr', 'tss', 'tes', 'strand', 'ensg', 'hgnc', 'remove')

## 1. assign genes into blocks
genes <- read.table(gene.file, col.names = gene.cols) %>%
    mutate(loc = ifelse(strand == '+', tss, tes), idx = 1:n()) %>%
        mutate(block.size = block.size0, block = round(loc/block.size)) %>%
            dplyr::select(-remove) %>%
                arrange(loc)

## 2. break large blocks into smaller chunks
while(TRUE) {
    big.blocks <- genes %>% group_by(block) %>%
        summarize(block.size = n()) %>%
            filter(block.size > max.size) %>% select(block)

    if(nrow(big.blocks) == 0) break

    ret1 <- genes %>% filter(! block %in% big.blocks$block)

    ret2 <-
        genes %>% right_join(big.blocks, by = 'block') %>%
            mutate(block.size = block.size / 2) %>%
                mutate(block = round(loc / block.size * 4) / 4)

    genes <- rbind(ret1, ret2) %>% arrange(loc)
}

out <- genes %>% group_by(block) %>%
    summarize(idx.vec = paste('c(',paste(idx, collapse = ','),')', sep=''))

cat(paste(out$idx.vec, collapse = '\n'), file = stdout())
cat('\n', file = stdout())
