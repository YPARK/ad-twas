#!/usr/bin/env Rscript
## Estimate QTL effects
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

data.hdr <- argv[1] # e.g., data.hdr = '/broad/hptmp/ypp/AD/twas/qtl/19/data-2'
out.hdr <- argv[2]  # e.g., out.hdr = 'temp2'

dir.create(dirname(out.hdr), recursive = TRUE, showWarnings = FALSE)

################################################################
source('util.R')

y1.data.file <- data.hdr %&&% '.y1.ft'
y0.data.file <- data.hdr %&&% '.y0.ft'
x.data.file <- data.hdr %&&% '.x.ft'
x.bim.data.file <- data.hdr %&&% '.x.bim.ft'
gene.data.file <- data.hdr %&&% '.y.gene.ft'
cov.data.file <- data.hdr %&&% '.cov.ft'
peer.data.file <- data.hdr %&&% '.peer.ft'

in.files <- c(y1.data.file,
              y0.data.file,
              x.data.file,
              x.bim.data.file,
              gene.data.file,
              cov.data.file,
              peer.data.file)

if(!all(sapply(in.files, file.exists))) {
    log.msg('Insufficient input files:\n%s\n',
            paste(in.files, collapse='\n'))
    q()
}

################################################################
library(dplyr)
library(feather)
library(fqtl)
.read.tab <- function(...) read_feather(...) %>% as.data.frame()
.read.mat <- function(...) read_feather(...) %>% as.matrix()
.cor <- function(...) cor(..., use = 'pairwise.complete.obs')

## standadize scale
stdize.count <- function(xx) {
    xx[xx < .5] <- NA
    xx.med <- apply(xx, 2, median, na.rm = TRUE)
    xx.scaled <- sweep(xx, 2, xx.med, `/`)
    ret <- xx.scaled * 10
}

Y1 <- .read.mat(y1.data.file)
Y0 <- .read.mat(y0.data.file)
X <- .read.mat(x.data.file) %>% scale()
x.bim <- .read.tab(x.bim.data.file)
colnames(X) <- x.bim[, 4]

qq <- qnorm(1:nrow(Y1) / (nrow(Y1) + 1))
make.quant <- function(y) {
    ret <- y
    ret[order(y)] <- qq
    return(ret)
}
Y1.quant <- apply(Y1, 2, make.quant)

Y1 <- stdize.count(Y1)
Y0 <- stdize.count(Y0)

genes <- .read.tab(gene.data.file) %>% select(-data.pos)
colnames(Y1) <- gene.names <- genes$ensg
colnames(Y1.quant) <- gene.names

V <- .read.tab(cov.data.file)
peer <- .read.mat(peer.data.file) %>% scale()

pheno <- V %>%
    select(amyloid_sqrt, tangles_sqrt, cogn_ep_random_slope, msex.y, studyn, age_death, educ) %>%
        scale() %>% as.matrix()

pheno <- cbind(pheno, intercept = 1)

## 1. W = Y0 - X * theta
## 2. Y1 = f(W + R)
estimate.confounder <- function(yy, yy.ctrl,
                                xx = NULL,
                                clean.confounder = FALSE,
                                out.tag = '',
                                .pheno = pheno,
                                .iterations = 5000,
                                .factored = TRUE,
                                .gene.names = gene.names) {

    k <- pmax(ncol(yy) - 1, 1)

    vb.opt <- list(vbiter = ceiling(.iterations/2),
                   gammax = 1e4,
                   out.residual = TRUE,
                   tol = 1e-8,
                   rate = 1e-2,
                   decay = -1e-2,
                   pi.ub = -0,
                   pi.lb = -2,
                   model = 'nb',
                   svd.init = TRUE,
                   jitter = 0.01,
                   k = k)
    
    if(clean.confounder) {
        stopifnot(!is.null(xx))
        out.W <- fqtl.regress(y = yy.ctrl,
                              x.mean = xx,
                              c.mean = .pheno,
                              factored = FALSE,
                              options = vb.opt)
        
        W.proxy <- out.W$resid$theta %>% scale()
    } else {
        W.proxy <- yy.ctrl %>% scale()
    }

    ## Remove Y1 ~ W + .Pheno + 1    
    out.1 <- fqtl.regress(y = yy,
                          x.mean = W.proxy, 
                          c.mean = .pheno,
                          factored = .factored,
                          options = vb.opt)

    residual <- out.1$resid$theta
    residual[is.na(yy)] <- NA
    colnames(residual) <- .gene.names
    pve <- cbind(.gene.names, take.pve(out.1, W.proxy, .pheno))

    list(R = residual, W = W.proxy, PVE = pve,
         mean.left = out.1$mean.left,
         mean.right = out.1$mean.right,
         mean = out.1$mean,
         out.tag = out.tag)
}

estimate.lm <- function(yy, yy.ctrl,
                        out.tag = '',
                        .pheno = pheno,
                        .gene.names = gene.names) {

    .pheno[is.na(.pheno)] <- 0
    lm.out <- lm(yy ~ yy.ctrl + .pheno -1)

    out <- list()
    out$resid$theta <- lm.out$residuals
    out$mean$theta <- lm.out$coefficients %r% 1:ncol(yy.ctrl) %>% as.matrix()
    out$mean.cov$theta <- lm.out$coefficients %r% (ncol(yy.ctrl) + 1:ncol(.pheno)) %>% as.matrix()

    pve <- cbind(.gene.names, take.pve(out, as.matrix(yy.ctrl), as.matrix(.pheno)))
    
    list(R = lm.out$residuals, PVE = pve, out.tag = out.tag, lm = lm.out)
}

write.confounder <- function(.obj, .out.hdr) {
    resid.out.file <- .out.hdr %&&% '.resid.gz'
    pve.out.file <- .out.hdr %&&% '.pve.gz'
    write.tab(.obj$R, file = gzfile(resid.out.file))
    write.tab(.obj$PVE, file = gzfile(pve.out.file))
}

write.tab.gz <- function(.tab, .out.file) {
    write.tab(.tab, file = gzfile(.out.file))
}

filter.qtl <- function(qtl.tab, cis.dist = 1e6) {

    qtl.tab %>% mutate(ensg = as.character(gene)) %>%
        left_join(genes, by = 'ensg') %>%
            filter(snp > (tss - cis.dist)) %>%
                filter(snp < (tes + cis.dist)) %>%
                    select(snp, gene, beta, beta.z)
}

Y0.std <- log2(0.5 + Y0) %>% scale()
Y0.std[is.na(Y0.std)] <- 0

conf.1 <- estimate.confounder(Y1, Y0.std, .factored = TRUE, out.tag = 'hs-fqtl')
conf.2 <- estimate.confounder(Y1, Y0.std, .factored = FALSE, out.tag = 'hs-sqtl')
conf.3 <- estimate.confounder(Y1, peer, .factored = FALSE, out.tag = 'peer-sqtl')
conf.4 <- estimate.lm(Y1.quant, Y0.std, out.tag = 'hs-lm')
conf.5 <- estimate.lm(Y1.quant, peer, out.tag = 'peer-lm')
conf.list <- list(conf.1, conf.2, conf.3, conf.4, conf.4, conf.5)

sapply(conf.list, function(cc) write.confounder(cc, out.hdr %&&% '-' %&&% cc$out.tag))

qtl.0 <- get.marginal.qtl(X, Y1.quant) %>% filter.qtl()
write.tab.gz(qtl.0, out.hdr %&&% '.qtl-raw.gz')

sapply(conf.list, function(cc) {
    write.tab.gz(get.marginal.qtl(X, cc$R) %>% filter.qtl(),
                 out.hdr %&&% '.qtl-' %&&% cc$out.tag %&&% '.gz')
})

x.perm <- X %r% sample(nrow(X)) %>% as.matrix()
sapply(conf.list, function(cc) {
    write.tab.gz(get.marginal.qtl(x.perm, cc$R) %>% filter.qtl(),
                 out.hdr %&&% '.perm-' %&&% cc$out.tag %&&% '.gz')
})

write.tab.gz(x.bim, out.hdr %&&% '.snps.gz')
write.tab.gz(genes, out.hdr %&&% '.genes.gz')

log.msg('Finished\n')
