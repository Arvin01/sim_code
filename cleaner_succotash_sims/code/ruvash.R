## RUVASH ---------------------------------------------------------
ruvash_out <- vicar::ash_ruv4(Y = Y, X = X, k = num_sv, ctl = as.logical(control_genes),
                              limmashrink = TRUE, likelihood = "normal")
betahat <- ruvash_out$PosteriorMean
lfdr    <- ruvash_out$lfdr
pi0hat  <- ruvash_out$fitted.g$pi[1]
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
