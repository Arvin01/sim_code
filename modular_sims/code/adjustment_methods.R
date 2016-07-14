#' @param calibrate a logical for whether to use MAD calibrated t-stats.
cate_nc <- function(Y, X, num_sv, control_genes, calibrate = FALSE) {
    calibrate <- as.logical(calibrate)
    cate_nc <- cate::cate(~Treatment, Y = Y,
                          X.data = data.frame(Treatment = X[, 2]),
                          r = num_sv, fa.method = "ml", adj.method = "nc",
                          nc = as.logical(control_genes), calibrate = calibrate)

    betahat   <- c(cate_nc$beta)
    sebetahat <- c(sqrt(cate_nc$beta.cov.row * cate_nc$beta.cov.col) /
                   sqrt(nrow(X)))
    df        <- rep(Inf, length = ncol(Y))
    pvalues   <- c(cate_nc$beta.p.value)
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

#' @param calibrate a logical for whether to use MAD calibrated t-stats.
cate_rr <- function(Y, X, num_sv, calibrate = FALSE) {
    calibrate <- as.logical(calibrate)
    cate_rr <- cate::cate(~Treatment, Y = Y,
                          X.data = data.frame(Treatment = X[, 2]),
                          r = num_sv, fa.method = "ml", adj.method = "rr",
                          calibrate = calibrate)
    cate_rr_out           <- list()
    betahat   <- c(cate_rr$beta)
    sebetahat <- c(sqrt(cate_rr$beta.cov.row * cate_rr$beta.cov.col) / sqrt(nrow(X)))
    pvalues   <- c(cate_rr$beta.p.value)
    df        <- rep(Inf, length = ncol(Y))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

ols <- function(Y, X) {
    limma_out <- limma::lmFit(object = t(Y), design = X)
    betahat   <- limma_out$coefficients[, 2]
    sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
    df        <- limma_out$df.residual
    tstats    <- betahat / sebetahat
    pvalues   <- 2 * pt(-abs(tstats), df = df)
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

ruv2 <- function(Y, X, num_sv, control_genes) {
    ruv_ruv2 <- ruv::RUV2(Y = Y, X = as.matrix(X[, 2]),
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = as.matrix(X[, -2]))
    sebetahat <- sqrt(ruv_ruv2$sigma2 * ruv_ruv2$multiplier)
    betahat   <- c(ruv_ruv2$betahat)
    pvalues   <- c(ruv_ruv2$p)
    df        <- rep(ruv_ruv2$df, length = ncol(Y))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

ruv4 <- function(Y, X, num_sv, control_genes) {
    ruv_ruv4 <- ruv::RUV4(Y = Y, X = as.matrix(X[, 2]),
                          ctl = as.logical(control_genes),
                          k = num_sv, Z = as.matrix(X[, -2]))
    sebetahat <- sqrt(ruv_ruv4$sigma2 * ruv_ruv4$multiplier)
    betahat   <- c(ruv_ruv4$betahat)
    pvalues   <- c(ruv_ruv4$p)
    df        <- rep(ruv_ruv4$df, length = ncol(Y))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

ruvinv <- function(Y, X, control_genes) {
    ruv_ruvinv <- ruv::RUVinv(Y = Y, X = as.matrix(X[, 2]),
                              ctl = as.logical(control_genes),
                              Z = as.matrix(X[, -2]))
    sebetahat <- sqrt(ruv_ruvinv$sigma2 * ruv_ruvinv$multiplier)
    betahat   <- c(ruv_ruvinv$betahat)
    pvalues   <- c(ruv_ruvinv$p)
    df        <- rep(ruv_ruvinv$df, length = ncol(Y))
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

succotash <- function(Y, X, num_sv) {
    succ_out <- succotashr::succotash(Y = Y, X = X, k = num_sv,
                                      fa_method = "pca", num_em_runs = 3,
                                      optmethod = "em",
                                      two_step = TRUE,
                                      var_scale_init_type = "null_mle",
                                      z_init_type = "null_mle")
    betahat <- succ_out$betahat
    lfdr    <- succ_out$lfdr
    pi0hat  <- succ_out$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}

sva <- function(Y, X, num_sv) {
    trash     <- capture.output(sva_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv))
    X.sv      <- cbind(X, sva_out$sv)
    limma_out <- limma::lmFit(object = t(Y), design = X.sv)
    betahat   <- limma_out$coefficients[, 2]
    sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
    df        <- limma_out$df.residual
    tstats    <- betahat / sebetahat
    pvalues   <- 2 * pt(-abs(tstats), df = df)
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}

vruv4 <- function(Y, X, num_sv, control_genes, adjust_bias = FALSE) {
    vout <- vicar::vruv4(Y = Y, X = X, k = num_sv, ctl = as.logical(control_genes),
                         limmashrink = TRUE, cov_of_interest = 2, adjust_bias = adjust_bias)
    betahat   <- c(vout$betahat)
    sebetahat <- c(vout$sebetahat)
    pvalues   <- c(vout$pvalues)
    df        <- vout$degrees_freedom
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}
