cate_rr <- cate::cate(~Treatment, Y = Y,
                      X.data = data.frame(Treatment = X[, 2]),
                      r = num_sv, fa.method = "ml", adj.method = "rr",
                      calibrate = FALSE)
cate_rr_out           <- list()
cate_rr_out$betahat   <- cate_rr$beta
cate_rr_out$sebetahat <- sqrt(cate_rr$beta.cov.row * cate_rr$beta.cov.col) / sqrt(nrow(X))
cate_rr_out$pvalue    <- cate_rr$beta.p.value
cate_qv               <- fit_freq_methods(out_obj = cate_rr_out)

betahat <- cate_rr_out$betahat
lfdr    <- cate_qv$q_storey$lfdr
pi0hat  <- cate_qv$q_storey$pi0
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
