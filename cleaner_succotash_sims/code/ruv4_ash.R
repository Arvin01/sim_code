## RUV4 then ASH -------------------------------------------------------
ruv_ruv4 <- ruv::RUV4(Y = Y, X = as.matrix(X[, ofinterest]),
                      ctl = as.logical(control_genes),
                      k = num_sv, Z = as.matrix(X[, -ofinterest]))
ruv4sebetahat <- sqrt(ruv_ruv4$sigma2 * ruv_ruv4$multiplier)
ashr_ruv4_out <- ashr::ash(betahat = ruv_ruv4$betahat, ruv4sebetahat)

betahat <- ashr_ruv4_out$PosteriorMean
lfdr    <- ashr_ruv4_out$lfdr
pi0hat  <- ashr_ruv4_out$fitted.g$pi[1]
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
