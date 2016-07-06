cate_nc <- cate::cate(~Treatment, Y = Y,
                      X.data = data.frame(Treatment = X[, 2]),
                      r = num_sv_cate, fa.method = "ml", adj.method = "nc",
                      nc = as.logical(control_genes), calibrate = TRUE)
cate_nc_out           <- list()
cate_nc_out$betahat   <- cate_nc$beta
cate_nc_out$sebetahat <- sqrt(cate_nc$beta.cov.row * cate_nc$beta.cov.col) /
    sqrt(nrow(X))
cate_nc_out$pvalue    <- cate_nc$beta.p.value
cate_nc_qv            <- fit_freq_methods(out_obj = cate_nc_out)

betahat_df$cate_nc <- cate_nc_out$betahat
lfdr_df$cate_nc    <- cate_nc_qv$q_storey$lfdr
pi0hat_vec         <- c(pi0hat_vec, cate_nc_qv$q_storey$pi0)
