## RUV2 then ASH -------------------------------------------------------
ruv_ruv2 <- ruv::RUV2(Y = Y, X = as.matrix(X[, ofinterest]),
                      ctl = as.logical(control_genes),
                      k = num_sv, Z = as.matrix(X[, -ofinterest]))
ruv2sebetahat <- sqrt(ruv_ruv2$sigma2 * ruv_ruv2$multiplier)
ashr_ruv2_out <- ashr::ash(betahat = ruv_ruv2$betahat, ruv2sebetahat)

betahat <- ashr_ruv2_out$PosteriorMean
lfdr    <- ashr_ruv2_out$lfdr
pi0hat  <- ashr_ruv2_out$fitted.g$pi[1]
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
