#!/usr/bin/env Rscript
## Estimate QTL effects with polygenic model

argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 3) q()

y.file <- argv[1] # e.g., y.file = 'result/qtl/3/b12-hs-lm.resid.gz'
data.hdr <- argv[2] # e.g.,  data.hdr = '/broad/hptmp/ypp/AD/twas/qtl//3/data-12'
out.hdr <- argv[3] # e.g., out.hdr <- 'temp'

do.perm <- FALSE
if(length(argv) > 3){
    rseed <- as.integer(argv[4])
    do.perm <- TRUE
    set.seed(rseed)
}

lodds.cutoff <- -2

################################################################
library(dplyr)
library(feather)
library(fqtl)
library(Matrix)
library(methods)
source('util.R')

.read.ft.tab <- function(...) read_feather(...) %>% as.data.frame()
.read.ft.mat <- function(...) read_feather(...) %>% as.matrix()
.cor <- function(...) cor(..., use = 'pairwise.complete.obs')

x.data.file <- data.hdr %&&% '.x.ft'
x.bim.data.file <- data.hdr %&&% '.x.bim.ft'
gene.data.file <- data.hdr %&&% '.y.gene.ft'
cov.data.file <- data.hdr %&&% '.cov.ft'
peer.data.file <- data.hdr %&&% '.peer.ft'

out.max.file <- out.hdr %&&% '.max-effect.gz'
out.effect.file <- out.hdr %&&% '.effect.gz'

X <- .read.ft.mat(x.data.file) %>% scale() %>% as.matrix()
x.bim <- .read.ft.tab(x.bim.data.file)[, -3]
colnames(x.bim) <- c('chr', 'snp', 'snp.loc', 'a1', 'a2')
genes <- .read.ft.tab(gene.data.file) %>% select(-data.pos)

if(do.perm) {
    X <- X %r% sample(nrow(X))
    log.msg('Permuted X with Rseed = %d\n\n', rseed)
}

Y <- read.table(y.file) %>% scale() %>% as.matrix()
valid.genes <- which(apply(!is.na(Y), 2, mean) > .8)

if(length(valid.genes) < 1) {
    log.msg('No valid genes\n\n')
    system("printf \"\" | gzip > " %&&% out.max.file)
    system("printf \"\" | gzip > " %&&% out.effect.file)
    q()
}

genes <- genes %r% valid.genes
Y <- Y %c% valid.genes

snp.loc <- x.bim$snp.loc
ensg <- genes$ensg

sparse2triple <- function(.sp, val.name = '', .r = snp.loc, .c = ensg) {
    require(Matrix)
    require(dplyr)
    ret <- Matrix::summary(.sp)        
    ret <- data.frame(ret, snp.loc = .r[ret$i], ensg = .c[ret$j]) %>%
        select(-i, -j)
    colnames(ret)[1] <- val.name
    return(ret)
}

melt.effect <- function(mean.obj) {
    require(dplyr)
    ret <- sparse2triple(mean.obj$theta, 'theta') %>%
        select(ensg, snp.loc, theta) %>% 
            left_join(sparse2triple(sqrt(mean.obj$theta.var), 'se'), by = c('snp.loc', 'ensg')) %>%
                left_join(sparse2triple(mean.obj$lodds, 'lodds'), by = c('snp.loc', 'ensg'))
}

################################################################
vb.opt <- list(vbiter = 2500,
               gammax = 1e4,
               tol = 1e-8,
               rate = 1e-2,
               decay = -1e-2,
               pi.ub = -0,
               pi.lb = -2,
               out.residual = FALSE,
               model = 'gaussian')

out.1 <- fqtl.regress(y = Y, x.mean = X,
                      y.loc = genes$tss,
                      y.loc2 = genes$tes,
                      x.mean.loc = snp.loc,
                      options = vb.opt)

X.safe <- X
X.safe[is.na(X)] <- 0
v.g <- apply(X.safe %*% out.1$mean$theta, 2, var, na.rm = TRUE)
v.tot <- apply(Y, 2, var, na.rm = TRUE)
h2.g <- v.g / v.tot

effect.tot <- melt.effect(out.1$mean) %>%
        left_join(x.bim, by = 'snp.loc') %>%
            left_join(genes %>% select(-chr), by = 'ensg')

## gene-wise maximum stat
effect.max <- effect.tot %>%
    group_by(ensg) %>% slice(which.max(lodds))

effect.out <- effect.tot %>%
    filter(lodds > lodds.cutoff)

write.tab(effect.max, file = gzfile(out.max.file))
write.tab(effect.out, file = gzfile(out.effect.file))
