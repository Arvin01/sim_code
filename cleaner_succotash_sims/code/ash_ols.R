## function(Y, X, num_sv, control_genes = NULL) {
ols_fit <- get_ols(log_counts = t(Y), condition = X[, 2])
ols_out <- fit_freq_methods(out_obj = ols_fit)
ash_ols <- ashr::ash(betahat = ols_fit$betahat, sebetahat = ols_fit$sebetahat)
betahat <- ash_ols$PosteriorMean
lfdr    <- ash_ols$lfdr
pi0hat  <- ash_ols$fitted.g$pi[1]
return_list <- list(betahat = betahat,
                    lfdr = lfdr,
                    pi0hat = pi0hat)
return_list
