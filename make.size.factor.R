#!/usr/bin/env Rscript
## Estimate size factor of samples following Anders & Huber (2010)

source('util.R')

argv <- commandArgs(trailingOnly = TRUE)

data.hdr <- argv[1]      # e.g., data.hdr <- 'data/rnaseq/chr'
data.suffix <- argv[2]   # e.g., data.suffix <- '-count.txt.gz'
pheno.file <- argv[3]    # e.g., pheno.file <- 'data/pheno.txt.gz'
out.file <- argv[4]

data.files <- data.hdr %&&% 1:22 %&&% data.suffix
pheno.tab <- read.table(pheno.file, header = TRUE)

Y <- do.call(cbind, lapply(data.files, read.table))

## 1. Treat 0 as missing values
remove.zero <- function(Y) {
    Y[ Y < 1 ] <- NA
    return(Y)    
}

Y <- remove.zero(Y)

## 2. take gene-wise geometric mean
g.mean <- exp(apply(log(Y), 2, mean, na.rm = TRUE))

## 3. size.factor[i] = median_g { K[i, g] / geometric mean_j K[j, g] }
size.factor <- apply(sweep(Y, 2, g.mean, `/`), 1, median, na.rm = TRUE)

out <- cbind(data.col = pheno.tab$data.col, size.factor)
write.mat(data.frame(out), file = gzfile(out.file))
