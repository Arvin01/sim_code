## RUV4 -------------------------------------------------------
ruv_ruv4 <- ruv::RUV4(Y = Y, X = as.matrix(X[, ofinterest]),
                      ctl = as.logical(control_genes),
                      k = num_sv, Z = as.matrix(X[, -ofinterest]))
ruv_ruv4_list         <- list()
ruv_ruv4_list$betahat <- c(ruv_ruv4$betahat)
ruv_ruv4_list$pvalue  <- c(ruv_ruv4$p)
ruv_ruv4_list$df      <- rep(ruv_ruv4$df, length = ncol(Y))
ruv4_qv               <- fit_freq_methods(out_obj = ruv_ruv4_list)

betahat <- ruv_ruv4_list$betahat
lfdr    <- ruv4_qv$q_storey$lfdr
pi0hat  <- c(pi0hat_vec, ruv4_qv$q_storey$pi0)
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
