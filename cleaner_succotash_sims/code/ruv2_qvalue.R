## RUV2-------------------------------------------------------
ruv_ruv2 <- ruv::RUV2(Y = Y, X = as.matrix(X[, ofinterest]),
                      ctl = as.logical(control_genes),
                      k = num_sv, Z = as.matrix(X[, -ofinterest]))
ruv_ruv2_list         <- list()
ruv_ruv2_list$betahat <- c(ruv_ruv2$betahat)
ruv_ruv2_list$pvalue  <- c(ruv_ruv2$p)
ruv_ruv2_list$df      <- rep(ruv_ruv2$df, length = ncol(Y))
ruv2_qv               <- fit_freq_methods(out_obj = ruv_ruv2_list)

betahat <- ruv_ruv2_list$betahat
lfdr    <- ruv2_qv$q_storey$lfdr
pi0hat  <- c(pi0hat_vec, ruv2_qv$q_storey$pi0)
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
