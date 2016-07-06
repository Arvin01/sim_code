vout <- vicar::vruv4(Y = Y, X = X, k = num_sv, ctl = as.logical(control_genes),
                     limmashrink = TRUE, cov_of_interest = 2)
betahat   <- c(vout$betahat)
sebetahat <- c(vout$sebetahat)
pvalues   <- c(vout$pvalues)
df        <- vout$degrees_freedom
