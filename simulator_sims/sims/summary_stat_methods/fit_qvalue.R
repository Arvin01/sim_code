fit_qvalue <- function(pvalues) {
    qout    <- qvalue::qvalue(p = pvalues)
    lfdr    <- c(qout$lfdr)
    pi0hat  <- c(qout$pi0)
    betahat <- c(betahat)
    return(list(lfdr = lfdr, pi0hat = pi0hat))
}
