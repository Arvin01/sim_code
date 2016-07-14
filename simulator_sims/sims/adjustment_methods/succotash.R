succotash <- function(Y, X, num_sv) {
    succ_out <- succotashr::succotash(Y = Y, X = X, k = num_sv,
                                      fa_method = "pca", num_em_runs = 3,
                                      optmethod = "em",
                                      two_step = TRUE,
                                      var_scale_init_type = "null_mle",
                                      z_init_type = "null_mle")
    betahat <- succ_out$betahat
    lfdr    <- succ_out$lfdr
    pi0hat  <- succ_out$pi0
    return(list(betahat = betahat, lfdr = lfdr, pi0hat = pi0hat))
}
