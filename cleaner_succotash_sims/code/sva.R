trash      <- capture.output(sva_out <- sva::sva(dat = t(Y), mod = X, n.sv = num_sv))
X.sv       <- cbind(X, sva_out$sv)
limma_out <- limma::lmFit(object = t(Y), design = X.sv)
betahat   <- limma_out$coefficients[, 2]
sebetahat <- limma_out$stdev.unscaled[, 2] * limma_out$sigma
df        <- limma_out$df.residual
tstats    <- betahat / sebetahat
pvalues   <- 2 * pt(-abs(tstats), df = df)
