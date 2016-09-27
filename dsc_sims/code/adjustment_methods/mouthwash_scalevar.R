vout <- vicar::mouthwash(Y = Y, X = X, k = num_sv, include_intercept = FALSE,
                         limmashrink = TRUE, likelihood = "t", mixing_dist = "uniform",
                         scale_var = TRUE)
betahat <- vout$result$PosteriorMean
lfdr    <- vout$result$lfdr
pi0hat  <- vout$pi0
