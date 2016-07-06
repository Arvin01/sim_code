suppressWarnings(ashout <- ashr::ash(betahat = betahat, sebetahat = sebetahat, df = df, model = model, outputlevel = 3))
betahat <- c(ashout$PosteriorMean)
lfdr    <- c(ashout$lfdr)
pi0hat  <- c(ashr::get_pi0(ashout))
