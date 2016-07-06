trash      <- capture.output(sva_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv))
X.sv       <- cbind(X, sva_out$sv)
sva_ols    <- get_ols(log_counts = t(Y), condition = X.sv[, -1])
sva_ols_qv <- fit_freq_methods(out_obj = sva_ols)

betahat <- sva_ols$betahat
lfdr    <- sva_ols_qv$q_storey$lfdr
pi0hat  <- sva_ols_qv$q_storey$pi0
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
