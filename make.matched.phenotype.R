#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

pheno.file <- argv[1] # e.g., 'phenotypes/pheno_cov_n3033_032315.csv'
sample.file <- argv[2] # e.g., 'data/qc.samples.txt.gz'
out.file <- argv[3]

options(stringsAsFactors = FALSE)
library(dplyr)

pheno.tab <- read.table(pheno.file, header = TRUE, sep = ',')
pheno.tab[pheno.tab == -9] <- NA

sample.tab <- read.table(sample.file, sep = '\t',
                         col.names = c('sample.id', 'data.col', 'projid', 'RIN'))

pheno.out <- sample.tab %>% left_join(pheno.tab, by = 'projid') %>%
    select(data.col, projid, amyloid_sqrt, tangles_sqrt, cogn_ep_random_slope, age_death, studyn, msex, educ) %>%
        filter(!is.na(amyloid_sqrt)) %>%
            group_by(projid) %>%
                filter(row_number() == 1)

write.table(pheno.out, file = gzfile(out.file), sep = '\t', row.names = FALSE,
            col.names = TRUE, quote = FALSE)
