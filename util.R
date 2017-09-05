options(stringsAsFactors = FALSE)

log.msg <- function(...) {
    cat('[', date() ,']', sprintf(...), file = stderr())
}

`%&&%` <- function(a, b) paste(a, b, sep = '')

`%c%` <- function(a, b) a[, b, drop = FALSE]

`%r%` <- function(a, b) a[b, , drop = FALSE]

.NA <- function(nrow, ncol) {
    matrix(NA, nrow, ncol)
}

.rnorm <- function(nrow, ncol) {
    matrix(rnorm(nrow * ncol), nrow, ncol)
}

fast.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / n.obs
}

fast.z.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / sqrt(n.obs)
}

get.marginal.qtl <- function(xx, yy) {
    require(reshape2)
    require(dplyr)

    .xx <- scale(xx)
    .yy <- scale(yy)

    beta <- fast.cov(.xx, .yy) %>% as.matrix() %>% melt()
    beta.z <- fast.z.cov(.xx, .yy) %>% as.matrix() %>% melt()

    colnames(beta) <- c('rs', 'cg', 'beta')
    colnames(beta.z) <- c('rs', 'cg', 'beta.z')
    out <- beta %>% left_join(beta.z, by = c('rs', 'cg')) %>%
        as.data.frame()
    return(out)
}

fast.cor <- function(x, y) {
    x.sd <- apply(x, 2, sd, na.rm = TRUE)
    y.sd <- apply(y, 2, sd, na.rm = TRUE)
    ret <- fast.cov(x, y)
    ret <- sweep(sweep(ret, 1, x.sd, `/`), 2, y.sd, `/`)    
}

write.tab <- function(x, ...) {
    write.table(x, sep = '\t', quote = FALSE,
                row.names = FALSE, col.names = FALSE, ...)
}

write.mat <- function(x, digits = 4, ...) {
    write.tab(round(x, digits), ...)
}
