## function(Y, X, num_sv, control_genes = NULL) {
ols_fit <- get_ols(log_counts = t(Y), condition = X[, 2])
ols_out <- fit_freq_methods(out_obj = ols_fit)
betahat <- ols_fit$betahat
lfdr    <- ols_out$q_storey$lfdr
pi0hat  <- ols_out$q_storey$pi0
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
