#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = TRUE)

pheno.file <- argv[1] # e.g., 'data/pheno.txt.gz'
size.factor.file <- argv[2] # e.g., 'data/rnaseq/size_factors.txt.gz'
out.file <- argv[3]

library(dplyr)
library(vsn)
library(peer)
source('util.R')
options(stringsAsFactors = FALSE)

pheno.tab <- read.table(pheno.file, header = TRUE)
size.factor <- read.table(size.factor.file)[, 2]
data.files <- 'data/rnaseq/chr' %&&% 1:22 %&&% '-count.txt.gz'

## Take most variable 1000 genes and adjust size factor
Y <- do.call(cbind, lapply(data.files, read.table))
Y <- sweep(Y, 1, size.factor, `/`)

y.sd <- apply(Y, 2, sd, na.rm = TRUE)
idx <- order(y.sd, decreasing = TRUE)[1:1000]
Y <- Y %c% idx

## log2 transformation followed by vsn
Y.log2 <- log2(0.5 + t(Y))
vsn.fit <- vsn2(Y.log2)

Y.vsn <- predict(vsn.fit, Y.log2) %>% t() %>% scale()
C <- pheno.tab %>% select(-data.col, -projid) %>%
    scale() %>% as.matrix()
C[is.na(C)] <- 0

peer.model <- PEER()
PEER_setPhenoMean(peer.model, Y.vsn)
PEER_setCovariates(peer.model, C)
PEER_setNk(peer.model, 30)
PEER_update(peer.model)

## confounder matrix
W <- PEER_getX(peer.model) %c% -(1:ncol(C))
W <- W %c% order(apply(W, 2, sd), decreasing = TRUE)

write.mat(W, file = gzfile(out.file))
