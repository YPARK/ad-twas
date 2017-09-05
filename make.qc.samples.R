#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

data.sample.file <- argv[1] # e.g., 'data/rnaseq.samples.txt.gz'
sample.info.file <- argv[2] # e.g., 'sample.info.txt'
out.file <- argv[3]

options(stringsAsFactors = FALSE, check.names = FALSE)
library(dplyr)

rm.3rd.dash <- function(x) paste(strsplit(x, split = '_')[[1]][1:2], collapse = '_')

sample.info.tab <- read.table(sample.info.file, header = TRUE, sep = '\t') %>%
    mutate(sample.id = sapply(collaborator.sample.id, rm.3rd.dash)) %>%
        mutate(RIN = RIN.by.BSP) %>% 
            select(sample.id, Projid, RIN)

data.samples <- read.table(data.sample.file, col.names = c('sample.id', 'data.cols')) %>%
    mutate(sample.id = sapply(sample.id, rm.3rd.dash))

matched <- data.samples %>% left_join(sample.info.tab, by = 'sample.id') %>%
    na.omit() %>% filter(RIN >= 6.0) # e.g., GTEx pipeline

write.table(matched, file = gzfile(out.file),
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
