ruvinv <- function(Y, X) {
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
