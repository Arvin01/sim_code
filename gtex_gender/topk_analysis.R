library(ggplot2)
tissue_vec <- c("adiposetissue", "bladder", "bloodvessel", "breast",
                "colon", "kidney", "lung", "nerve", "pancreas",
                "skin", "spleen", "adrenalgland", "blood", "brain",
                "esophagus", "heart", "liver", "muscle", "pituitary",
                "salivarygland", "smallintestine", "stomach", "thyroid")

num_sv_vec <- rep(NA, length = length(tissue_vec))
for (tissue_index in 1:length(tissue_vec)) {
    current_tissue <- tissue_vec[tissue_index]
    ## large dat ------------------------------------------------------------------
    dat <- readRDS(paste0("./output/cleaned_gtex_data/", current_tissue, ".Rds"))
    onsex <- dat$chrom == "X" | dat$chrom == "Y"

    dat$ctl[onsex] <- FALSE

    num_sv <- sva::num.sv(dat = dat$Y, mod = dat$X)
    num_sv_vec[tissue_index] <- num_sv
    cat(num_sv, "\n")
    ruv2out <- ruv::RUV2(Y = t(dat$Y), X = dat$X[, 2, drop = FALSE], ctl = dat$ctl,
                         k = num_sv, Z = dat$X[, -2, drop = FALSE])
    ruv3out <- vicar::ruv3(Y = t(dat$Y), X = dat$X, ctl = dat$ctl, cov_of_interest = ncol(dat$X),
                           k = num_sv, include_intercept = FALSE)
    ruv4out <- ruv::RUV4(Y = t(dat$Y), X = dat$X[, 2, drop = FALSE], ctl = dat$ctl,
                         k = num_sv, Z = dat$X[, -2, drop = FALSE])
    cateout <- cate::cate.fit(X.primary = dat$X[, 2, drop = FALSE],
                              X.nuis = dat$X[, -2, drop = FALSE],
                              Y = t(dat$Y), r = num_sv,
                              adj.method = "nc",
                              nc = dat$ctl, calibrate = FALSE)
    ruvbout <- vicar::ruvb(Y = t(dat$Y), X = dat$X, ctl = dat$ctl, k = num_sv,
                           cov_of_interest = ncol(dat$X), include_intercept = FALSE,
                           fa_func = vicar::bfa_gs_linked, return_mcmc = TRUE,
                           fa_args = list(use_code = "r", nsamp = 20000, thin = 20))
    saveRDS(object = ruvbout, file = paste0("./output/ruvbout/ruvbout_", current_tissue, ".Rds"))
}

ruvbt <- c(ruvbout$posterior_means / ruvbout$posterior_sd)
ruvbp <- pnorm(-abs(ruvbt)) * 2

pless <- pnorm(q = 0, mean = ruvbout$posterior_means, sd = ruvbout$posterior_sd)
lfsr2   <- pmin(pless, 1 - pless)

nplook <- 100
which_mostb <- order(ruvbp)[1:nplook]
which_mostcate <- order(cateout$beta.p)[1:nplook]


pdf(file = "cate_v_ruvb.pdf", height = 3, width = 3, family = "Times")
qplot(ruvbp, cateout$beta.p[!dat$ctl], size = I(0.2)) + xlab("RUVB") + ylab("CATE") +
    theme_bw() + geom_abline(slope = 1, intercept = 0, lty = 2)
dev.off()

## small dat ------------------------------------------------------------------
nsamps         <- 10
which_samps    <- sample(1:nrow(dat$X), size = nsamps)
smalldat       <- list()
smalldat$Y     <- dat$Y[, which_samps, drop = FALSE]
smalldat$X     <- dat$X[which_samps,, drop = FALSE]
smalldat$ctl   <- dat$ctl
smalldat$chrom <- dat$chrom
num_sv <- sva::num.sv(dat = smalldat$Y, mod = smalldat$X)
num_sv
ruv2out_small <- ruv::RUV2(Y = t(smalldat$Y),
                           X = smalldat$X[, 2, drop = FALSE],
                           ctl = smalldat$ctl, k = num_sv,
                           Z = smalldat$X[, -2, drop = FALSE])
ruv3out_small <- vicar::ruv3(Y = t(smalldat$Y), X = smalldat$X,
                             ctl = smalldat$ctl,
                             cov_of_interest = ncol(smalldat$X),
                             k = num_sv, include_intercept = FALSE)
ruv4out_small <- ruv::RUV4(Y = t(smalldat$Y),
                           X = smalldat$X[, 2, drop = FALSE],
                           ctl = smalldat$ctl, k = num_sv,
                           Z = smalldat$X[, -2, drop = FALSE])
cateout_small <- cate::cate.fit(X.primary = smalldat$X[, 2, drop = FALSE],
                                X.nuis = smalldat$X[, -2, drop = FALSE],
                                Y = t(smalldat$Y), r = num_sv,
                                adj.method = "nc", nc = smalldat$ctl)
ruvbout_small <- vicar::ruvb(Y = t(smalldat$Y), X = smalldat$X,
                             ctl = smalldat$ctl, k = num_sv,
                             cov_of_interest = ncol(smalldat$X),
                             include_intercept = FALSE,
                             fa_func = vicar::bfa_gs_linked,
                             return_mcmc = TRUE,
                             fa_args = list(use_code = "r", nsamp = 100, thin = 1))


sum((ruv2out_small$betahat - ruv2out$betahat) ^ 2)
sum((ruv3out_small$betahat - ruv3out$betahat) ^ 2)
sum((ruv4out_small$betahat - ruv4out$betahat) ^ 2)
sum((cateout_small$beta - cateout$beta) ^ 2)
sum((ruvbout_small$posterior_means - ruvbout$posterior_means) ^ 2)


ruvbt_small <- c(ruvbout_small$posterior_means / ruvbout_small$posterior_sd)
ruvbp_small <- 1 - 2 * pnorm(abs(ruvbt_small))

cor(ruvbp_small[which_mostb], ruvbp[which_mostb], method = "kendall")
cor(cateout_small$beta.p[which_mostcate], cateout$beta.p[which_mostcate], method = "kendall")



plot(ruvbt, cateout$beta.t[!dat$ctl])
abline(0, 1)



mean(c(ruvbout$posterior_lower < 0 & ruvbout$posterior_upper > 0))

r2lower <- ruv2out$betahat - 2 * sqrt(ruv2out$sigma2 * ruv2out$multiplier)
r2upper <- ruv2out$betahat + 2 * sqrt(ruv2out$sigma2 * ruv2out$multiplier)
mean(c(r2lower < 0 & r2upper > 0))


hist(c(ruvbout$lfsr))
hist(c(ruv4out$p))

ruv2order <- order(c(ruv2out$p)[!dat$ctl])
ruv3order <- order(c(ruv3out$pvalues_adjusted)[!dat$ctl])
ruv4order <- order(c(ruv4out$p)[!dat$ctl])
cateorder <- order(c(cateout$beta.p.value)[!dat$ctl])
ruvborder <- order(c(ruvbp))

r2 <- cumsum(onsex[!dat$ctl][ruv2order])
r3 <- cumsum(onsex[!dat$ctl][ruv3order])
r4 <- cumsum(onsex[!dat$ctl][ruv4order])
ca <- cumsum(onsex[!dat$ctl][cateorder])
rb <- cumsum(onsex[!dat$ctl][ruvborder])

numlook <- 700
plot(r2[1:numlook], type = "l")
lines(r3[1:numlook], col = 2)
lines(r4[1:numlook], col = 3)
lines(ca[1:numlook], col = 4)
lines(rb[1:numlook], col = 5)


c(33, 1, 32, 20, 14, 6, 25, 24, 19, 21, 15, 17, 20, 16, 20)
