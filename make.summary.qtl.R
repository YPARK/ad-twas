#!/usr/bin/env Rscript
## Estimate QTL effects
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) < 2) q()

data.hdr <- argv[1] # e.g., data.hdr = '/broad/hptmp/ypp/AD/twas/qtl/1/data-11'
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
library(zqtl)
library(fqtl)
library(glmnet)
library(methods)

.read.tab <- function(...) read_feather(...) %>% as.data.frame()
.read.mat <- function(...) read_feather(...) %>% as.matrix()
.cor <- function(...) cor(..., use = 'pairwise.complete.obs')

stdize.count <- function(xx) {
    x.min <- -4
    x.max <- 4
    xx <- log2(1/2 + xx) %>% scale()
    xx[xx > x.max] <- x.max
    xx[xx < x.min] <- x.min
    xx[is.na(xx)] <- x.min
    xx <- 2^(xx * 2)

    return(xx)
}

nb2gauss <- function(Y) {
    require(fqtl)
    require(dplyr)

    vb.opt <- list(vbiter = 3000,
                   gammax = 1e4,
                   out.residual = TRUE,
                   tol = 1e-8,
                   rate = 0.01,
                   decay = -0.01,
                   model = 'nb',
                   pi = -0,
                   tau = -10,
                   do.hyper = FALSE,
                   nsample = 5,
                   print.interv = 100,
                   adam.rate.m = 0.9,
                   adam.rate.v = 0.99)

    y <- stdize.count(Y)
    y.bar <- matrix(apply(y, 1, mean), ncol = 1)
    out <- fqtl.regress(y, x.mean = y.bar, options = vb.opt)

    ret <- out$resid$theta %>% scale() %>% as.matrix()
    ret[is.na(Y)] <- NA
    return(ret)
}

take.theta <- function(qtl.out) {
    if('mean.left' %in% names(qtl.out)){
        theta <- qtl.out$mean.left$theta %*% t(qtl.out$mean.right$theta)
    } else {
        theta <- qtl.out$mean$theta
    }
    return(theta)
}

take.pve <- function(qtl.out, xx, cc) {

    theta <- take.theta(qtl.out)
    theta.cov <- qtl.out$mean.cov$theta
    resid <- qtl.out$resid$theta

    xx[is.na(xx)] <- 0
    theta[is.na(theta)] <- 0
    theta.cov[is.na(theta.cov)] <- 0
    cc[is.na(cc)] <- 0
    eta.1 <- xx %*% theta
    eta.2 <- cc %*% theta.cov

    data.frame(v1 = apply(eta.1, 2, var, na.rm = TRUE) %>% as.vector(),
               v2 = apply(eta.2, 2, var, na.rm = TRUE) %>% as.vector(),
               vr = apply(resid, 2, var, na.rm = TRUE) %>% as.vector())
}

run.glmnet <- function(y, x, alpha = 1){
    valid <- !is.na(y)
    xx <- x[valid,,drop=FALSE]
    yy <- as.matrix(y[valid])

    cv.out <- cv.glmnet(x=xx, y=yy, alpha=alpha, nfolds=5)
    bet <- glmnet(x=xx, y=yy, alpha=alpha, lambda=cv.out$lambda.min)$beta
    bet <- matrix(bet)

    cat(mean(abs(bet) > 0), '\n')

    if(sum(abs(bet) > 0) == 0) {
        return(list(beta = bet, resid = y))
    }

    ## re-estimate lm on non-zero coefficients
    resid <- matrix(NA, nrow = length(y), ncol = 1)
    xx.sub <- xx %c% which(abs(bet)>0)
    lm.out <- lm(yy ~ xx.sub - 1)

    resid[valid,] <- as.matrix(lm.out$residuals)

    return(list(beta = bet, resid = resid))
}

estimate.lm <- function(yy, yy.ctrl,
                        out.tag = '',
                        .pheno = pheno,
                        .gene.names = gene.names) {

    .pheno[is.na(.pheno)] <- 0
    .nc <- ncol(yy.ctrl)
    .xx <- cbind(yy.ctrl, .pheno)
    .xx[is.na(.xx)] <- 0

    glmnet.list <- apply(yy, 2, run.glmnet, x = .xx, alpha = 0.5)

    out <- list()
    out$resid$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$resid)) %>% as.matrix()
    out$mean$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$beta %r% 1:.nc)) %>% as.matrix()
    out$mean.cov$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$beta %r% (-(1:.nc)))) %>% as.matrix()

    pve <- cbind(.gene.names, take.pve(out, as.matrix(yy.ctrl), as.matrix(.pheno)))
    colnames(out$resid$theta) <- .gene.names

    list(R = out$resid$theta, PVE = pve, out.tag = out.tag, lm = glmnet.list)
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

################################################################
Y1 <- .read.mat(y1.data.file)
Y0 <- .read.mat(y0.data.file)
genes <- .read.tab(gene.data.file) %>% select(-data.pos)

## basic Q/C : remove genes with too many zeros
valid.y1 <- which(apply(Y1 == 0, 2, mean) < .5)
Y1 <- Y1 %c% valid.y1
genes <- genes %r% valid.y1

if(ncol(Y1) == 0) {
    log.msg('No valid gene\n\n')
    system('printf "" | gzip > ' %&&% (out.hdr %&&% '.genes.gz'))
    q()
}
if(ncol(Y0) == 0) Y0 <- matrix(0, nrow = nrow(Y1), ncol = 1)

################################################################
X <- .read.mat(x.data.file) %>% scale()
x.bim <- .read.tab(x.bim.data.file)
colnames(X) <- x.bim[, 4]

################################################################
## Transform NB data to Gaussian fitting null model
Y.std <- nb2gauss(cbind(Y1, Y0))

Y1.std <- Y.std %c% (1:ncol(Y1))
Y0.std <- Y.std %c% (-(1:ncol(Y1)))

colnames(Y1.std) <- gene.names <- genes$ensg

peer <- .read.mat(peer.data.file) %>% scale()

pheno <- .read.tab(cov.data.file) %>%
    dplyr::select(amyloid_sqrt, tangles_sqrt, cogn_ep_random_slope, msex.y, studyn, age_death, educ) %>%
        scale() %>% as.data.frame() %>%
            dplyr::mutate(intercept = 1) %>%
                as.matrix()

################################################################
conf.1 <- estimate.lm(Y1.std, Y0.std, out.tag = 'hs-lm', pheno, gene.names)
conf.2 <- estimate.lm(Y1.std, peer, out.tag = 'peer-lm', pheno, gene.names)

conf.list <- list(conf.1, conf.2)

sapply(conf.list, function(cc) write.confounder(cc, out.hdr %&&% '-' %&&% cc$out.tag))


qtl.y1 <- get.marginal.qtl(X, Y1.std) %>% filter.qtl()
write.tab.gz(qtl.y1, out.hdr %&&% '.qtl-raw-y1.gz')

sapply(conf.list, function(cc) {
    write.tab.gz(get.marginal.qtl(X, cc$R) %>% filter.qtl(),
                 out.hdr %&&% '.qtl-' %&&% cc$out.tag %&&% '.gz')
})

write.tab.gz(x.bim, out.hdr %&&% '.snps.gz')
write.tab.gz(genes, out.hdr %&&% '.genes.gz')

log.msg('Finished\n')
