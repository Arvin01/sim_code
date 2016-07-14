ols <- function(Y, X, num_sv) {
    limma_out <- limma::lmFit(object = t(Y), design = X)
    betahat   <- limma_out$coefficients[, 2]
    sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
    df        <- limma_out$df.residual
    tstats    <- betahat / sebetahat
    pvalues   <- 2 * pt(-abs(tstats), df = df)
    return(list(betahat = betahat, sebetahat = sebetahat, df = df,
                pvalues = pvalues))
}
