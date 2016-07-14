
fit_ash <- function(betahat, sebetahat, df, model = "EE") {
    outputlevel <- 2
    suppressWarnings(ashout <- ashr::ash(betahat = betahat, sebetahat = sebetahat, df = df, model = model, outputlevel = outputlevel))
    if (outputlevel > 2) {
        betahat <- c(ashout$PosteriorMean)
    }
    lfdr    <- c(ashout$lfdr)
    pi0hat  <- c(ashr::get_pi0(ashout))
    return(list(lfdr = lfdr, pi0hat = pi0hat))
}
