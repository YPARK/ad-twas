#!/usr/bin/env Rscript
## Generate QTL data
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 6) {
    q()
}

expr.gene.file <- argv[1]                # e.g., expr.gene.file = 'data/rnaseq/chr2-genes.txt.gz'
expr.file <- argv[2]                     # e.g., expr.file = 'data/rnaseq/chr2-count.txt.gz'
expr.cols <- eval(parse(text = argv[3])) # e.g., expr.cols = 1:33
n.top.y0 <- as.integer(argv[4])          # e.g., n.top.y0 = 3
plink.hdr <- argv[5]                     # e.g., plink.hdr = 'rosmap-geno/gen/impute/rosmap1709-chr2'
out.hdr <- argv[6]                       # e.g., out.hdr = 'temp'

cis.dist <- 1e6

library(fqtl)
library(feather)
library(dplyr)
source('util.R')
options(stringsAsFactors = FALSE)

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)
dir.create(out.hdr, recursive = TRUE)
temp.dir <- system('mktemp -d ' %&&% out.hdr %&&% '/temp.XXXX', intern = TRUE, ignore.stderr = TRUE)

pheno.file <- 'data/pheno.txt.gz'
peer.file <- 'data/rnaseq/peer_factors.txt.gz'
size.factor.file <- 'data/rnaseq/size_factors.txt.gz'

################################################################
y1.out.file <- out.hdr %&&% '.y1.ft'
y0.out.file <- out.hdr %&&% '.y0.ft'
x.out.file <- out.hdr %&&% '.x.ft'
x.bim.out.file <- out.hdr %&&% '.x.bim.ft'
gene.out.file <- out.hdr %&&% '.y.gene.ft'
cov.out.file <- out.hdr %&&% '.cov.ft'
peer.out.file <- out.hdr %&&% '.peer.ft'

size.factor <- read.table(size.factor.file)[, 2]
peer.factor <- read.table(peer.file) %>% as.matrix()

################################################################
## Find correlated genes in other chromosomes

gene.cols <- c('chr', 'tss', 'tes', 'dir', 'ensg', 'hgnc', 'data.pos')
genes <- read.table(expr.gene.file, sep='\t', col.names=gene.cols)
genes <- genes %r% expr.cols
chr <- unique(genes$chr)

log.msg('genes:\n%s\n\n', paste(genes$hgnc, collapse=', '))

pheno.tab <- read.table(pheno.file, header = TRUE)

plink.lb <- max(min(genes$tss) - cis.dist, 0)
plink.ub <- max(genes$tes) + cis.dist

plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, chr, plink.lb, plink.ub, temp.dir %&&% '/plink')
system(plink.cmd)
plink <- read.plink(temp.dir %&&% '/plink')

colnames(plink$FAM) <- c('fid', 'iid', '.', '..', 'msex', '...')
fam.tab <- data.frame(plink$FAM %>% select(fid, iid, msex), geno.pos = 1:nrow(plink$FAM)) %>%
    mutate(projid = sapply(iid, gsub, pattern = 'ROS', replacement = '')) %>%
        mutate(projid = sapply(projid, gsub, pattern = 'MAP', replacement = '')) %>%
            mutate(projid = as.integer(projid))

if(nrow(plink$BIM) < 1) {
    log.msg('No variants in cis\n\n')
    q()
}

log.msg('Read genotypes\n\n')

Y1 <- read.table(expr.file, sep = '\t') %c% expr.cols %>% as.matrix()
colnames(Y1) <- genes$ensg
Y1 <- sweep(Y1, 1, size.factor, `/`)

take.y0 <- function(.chr, y1, n.top) {

    .expr.file <- gsub(expr.file, patter = 'chr' %&&% chr, replacement = 'chr' %&&% .chr)

    y0 <- read.table(.expr.file, sep = '\t') %>% as.matrix()
    y0 <- sweep(y0, 1, size.factor, `/`)

    valid.y0 <- which(apply(y0 == 0, 2, sum) == 0)
    y0 <- y0 %c% valid.y0

    log.msg('correlation between y1 and y0 in chr%d\n', .chr)

    y1.std <- scale(log2(0.5 + y1))
    y0.std <- scale(log2(0.5 + y0))

    abs.cor.mat <- t(abs(fast.cor(y1.std, y0.std)))

    y0.idx <- apply(abs.cor.mat, 2, function(x) order(x, decreasing = TRUE)[1:n.top])
    y0.idx <- unique(as.vector(y0.idx))

    ret <- y0[, y0.idx, drop = FALSE]
    rm(y0); rm(abs.cor.mat); rm(y0)
    log.msg('constructed y0 in chr%d\n', .chr)

    return(ret)
}

other.chr <- setdiff(1:22, chr)
y0.list <- lapply(other.chr, take.y0, y1 = Y1, n.top = n.top.y0)
Y0 <- do.call(cbind, y0.list)
gc()

## further refine Y0 -> n x (n.top * n.cpg) at most
abs.cor.y10 <- t(abs(fast.cor(scale(log2(0.5 + Y1)), scale(log2(0.5 + Y0)))))

y0.idx <- apply(abs.cor.y10, 2, function(x) order(x, decreasing=TRUE)[1:n.top.y0])
Y0.ref <- do.call(cbind, lapply(1:ncol(y0.idx), function(j) Y0 %c% y0.idx[, j]))

################################################################
## Take individuals with genotypes
samples <- pheno.tab %>% select(-data.col) %>%
    mutate(expr.pos = 1:nrow(pheno.tab))

sample.info <- samples %>%
    left_join(fam.tab, by = 'projid') %>%
        na.omit()

Y1 <- Y1 %r% sample.info$expr.pos %>% as.data.frame()
Y0 <- Y0.ref %r% sample.info$expr.pos %>% as.data.frame()
peer.out <- peer.factor %r% sample.info$expr.pos %>% as.data.frame()

x.bim <- plink$BIM
X <- plink$BED %r% sample.info$geno.pos %>% scale() %>% as.data.frame()
colnames(X) <- x.bim[, 2]

################################################################
## Write them down
write_feather(sample.info, path = cov.out.file)
write_feather(peer.out, path = peer.out.file)
write_feather(Y1, path = y1.out.file)
write_feather(Y0, path = y0.out.file)
write_feather(X, path = x.out.file)
write_feather(x.bim, path = x.bim.out.file)
write_feather(genes, path = gene.out.file)

system('rm -r ' %&&% temp.dir)
log.msg('Successfully generated data\n')

