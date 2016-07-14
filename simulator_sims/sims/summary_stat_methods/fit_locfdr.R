
fit_locfdr <- function(pvalues) {
    pvalues[pvalues == 0] <- min(pvalues[pvalues != 0])
    pvalues[pvalues == 1] <- max(pvalues[pvalues != 1])
    zvalues <- qnorm(p = pvalues)
    locout <- locfdr::locfdr(zvalues)
    lfdr <- locout$fdr
    return(lfdr)
}
