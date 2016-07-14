ruv4 <- function(Y, X, num_sv) {
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
