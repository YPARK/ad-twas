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
library(fqtl)
library(glmnet)
library(methods)

.read.tab <- function(...) read_feather(...) %>% as.data.frame()
.read.mat <- function(...) read_feather(...) %>% as.matrix()
.cor <- function(...) cor(..., use = 'pairwise.complete.obs')

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

    data.frame(v1 = apply(eta.1, 2, var, na.rm = TRUE),
               v2 = apply(eta.2, 2, var, na.rm = TRUE),
               vr = apply(resid, 2, var, na.rm = TRUE))
}

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

    vb.opt <- list(vbiter = .iterations,
                   gammax = 1e4,
                   out.residual = TRUE,
                   tol = 1e-8,
                   rate = 0.01,
                   decay = -0.01,
                   model = 'nb',
                   pi = -0,
                   tau = -5,
                   do.hyper = FALSE,
                   nsample = 10,
                   print.interv = 100,
                   adam.rate.m = 0.5,
                   adam.rate.v = 0.9)

    if(clean.confounder) {
        stopifnot(!is.null(xx))

        out.W <- fqtl.regress(y = yy.ctrl, x.mean = xx, factored = FALSE, options = vb.opt)
        W.proxy <- out.W$resid$theta
        W.proxy[is.na(yy.ctrl)] <- NA
        W.proxy <- scale(W.proxy)
    } else {
        W.proxy <- yy.ctrl %>% scale()
    }

    k <- max(min(ncol(yy) - 1, 5), 1)

    vb.opt$k <- k
    vb.opt$svd.init <- TRUE

    ## Remove Y1 ~ W + .Pheno + 1
    out.1 <- fqtl.regress(y = yy, x.mean = W.proxy, c.mean = .pheno, factored = .factored,
                          options = vb.opt)

    residual <- out.1$resid$theta
    residual[is.na(yy)] <- NA
    colnames(residual) <- .gene.names
    pve <- cbind(.gene.names, take.pve(out.1, W.proxy, .pheno))

    list(R = residual,
         W = W.proxy,
         PVE = pve,
         mean.left = out.1$mean.left,
         mean.right = out.1$mean.right,
         mean = out.1$mean,
         var = out.1$var,
         out.tag = out.tag)
}

run.glmnet <- function(y, x, alpha = 1){
    valid <- !is.na(y)
    xx <- x[valid,,drop=FALSE]
    yy <- as.matrix(y[valid])
    cv.out <- cv.glmnet(x=xx, y=yy, alpha=alpha, nfolds=5)
    bet <- glmnet(x=xx, y=yy, alpha=alpha, lambda=cv.out$lambda.min)$beta
    cat(mean(abs(bet) > 0), '\n')

    resid <- matrix(NA, nrow = length(y), ncol = 1)
    resid[valid,] <- as.matrix(yy - xx %*% bet)
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
    out$resid$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$resid))
    out$mean$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$beta %r% 1:.nc))
    out$mean.cov$theta <- do.call(cbind, lapply(glmnet.list, function(x) x$beta %r% (-(1:.nc))))

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

################################################################
Y1 <- .read.mat(y1.data.file)
Y0 <- .read.mat(y0.data.file)
genes <- .read.tab(gene.data.file) %>% select(-data.pos)

## basic Q/C : remove genes with too many zeros
valid.y1 <- which(apply(Y1 == 0, 2, mean) < .5)
Y1 <- Y1 %c% valid.y1
genes <- genes %r% valid.y1

if(ncol(Y1) == 0) q()

Y1 <- stdize.count(Y1)
Y0 <- stdize.count(Y0)

if(ncol(Y0) == 0) Y0 <- matrix(0, nrow = nrow(Y1), ncol = 1)

################################################################
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
Y1.quant[is.na(Y1)] <- NA

Y0.quant <- apply(Y0, 2, make.quant)
Y0.quant[is.na(Y0)] <- NA

colnames(Y1) <- gene.names <- genes$ensg
colnames(Y1.quant) <- gene.names

peer <- .read.mat(peer.data.file) %>% scale()

pheno <- .read.tab(cov.data.file) %>%
    dplyr::select(amyloid_sqrt, tangles_sqrt, cogn_ep_random_slope, msex.y, studyn, age_death, educ) %>%
        scale() %>% as.data.frame() %>%
            dplyr::mutate(intercept = 1) %>%
                as.matrix()

################################################################
conf.1 <- estimate.confounder(Y1, Y0, xx = X, clean.confounder = TRUE,
                              .factored = TRUE, out.tag = 'hs-fqtl')

conf.2 <- estimate.lm(Y1.quant, Y0.quant, out.tag = 'hs-lm')

conf.3 <- estimate.lm(Y1.quant, peer, out.tag = 'peer-lm')

conf.list <- list(conf.1, conf.2, conf.3)

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
