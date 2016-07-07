calibrate <- as.logical(calibrate)
cate_rr <- cate::cate(~Treatment, Y = Y,
                      X.data = data.frame(Treatment = X[, 2]),
                      r = num_sv, fa.method = "ml", adj.method = "rr",
                      calibrate = calibrate)
cate_rr_out           <- list()
betahat   <- c(cate_rr$beta)
sebetahat <- c(sqrt(cate_rr$beta.cov.row * cate_rr$beta.cov.col) / sqrt(nrow(X)))
pvalues   <- c(cate_rr$beta.p.value)
df        <- rep(Inf, length = ncol(Y))
