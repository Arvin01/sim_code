ruv_ruv2 <- ruv::RUV2(Y = Y, X = as.matrix(X[, 2]),
                      ctl = as.logical(control_genes),
                      k = num_sv, Z = as.matrix(X[, -2]))
sebetahat <- sqrt(ruv_ruv2$sigma2 * ruv_ruv2$multiplier)
betahat   <- c(ruv_ruv2$betahat)
pvalues   <- c(ruv_ruv2$p)
df        <- rep(ruv_ruv2$df, length = ncol(Y))
