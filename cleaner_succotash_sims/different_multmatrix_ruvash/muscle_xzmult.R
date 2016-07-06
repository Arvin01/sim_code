library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
    source("../code/datamaker_only_counts.R")
    args_val <- append(current_params, new_params)
    set.seed(new_params$current_seed)

    d_out <- datamaker_counts_only(args_val)

    which_null <- d_out$meta$null

    half_null <- which_null
    half_null[half_null == 1][sample(1:sum(which_null), size = sum(which_null) / 2)] <- 0

    beta_true <- rep(0, length = args_val$Ngene)
    beta_true[!which_null] <- d_out$meta$true_log2foldchange

    X <- as.matrix(model.matrix(~d_out$input$condition))
    colnames(X) <- c("Intercept", "Treatment")
    Y <- t(log2(as.matrix(d_out$input$counts + 1)))

    num_sv <- sva::num.sv(t(Y), mod = X, method = "be")

    start.time <- proc.time()
    ash_args <- list()
    ruvash_new_u <- vicar::ash_ruv4(Y = Y, X = X, ctl = as.logical(half_null),
                               k = num_sv,
                               likelihood = "normal",
                               limmashrink = FALSE)
    ruvash_new_s <- vicar::ash_ruv4(Y = Y, X = X, ctl = as.logical(half_null),
                               k = num_sv,
                               likelihood = "normal",
                               limmashrink = TRUE)
    ruvash_old_un <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                    k = num_sv,
                                    likelihood = "normal",
                                    limmashrink = FALSE,
                                    posthoc_inflate = FALSE)
    ruvash_old_sn <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                   k = num_sv,
                                   likelihood = "normal",
                                   limmashrink = TRUE,
                                   posthoc_inflate = FALSE)
    ruvash_old_ui <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                   k = num_sv,
                                   likelihood = "normal",
                                   limmashrink = FALSE,
                                   posthoc_inflate = TRUE)
    ruvash_old_si <- ashr::ash_ruv(Y = Y, X = X, ctl = as.logical(half_null),
                                   k = num_sv,
                                   likelihood = "normal",
                                   limmashrink = TRUE,
                                   posthoc_inflate = TRUE)

    ruvinv <- ruv::RUVinv(Y = Y, X = X[, 2, drop = FALSE], ctl = as.logical(half_null),
                          Z = X[, 1, drop = FALSE])
    ash_ruvinv <- ashr::ash(c(ruvinv$betahat), sqrt(ruvinv$sigma2 * ruvinv$multiplier))

    ruvrinv <- ruv::RUVrinv(Y = Y, X = X[, 2, drop = FALSE], ctl = as.logical(half_null),
                          Z = X[, 1, drop = FALSE])
    ash_ruvrinv <- ashr::ash(c(ruvrinv$betahat), sqrt(ruvrinv$sigma2 * ruvrinv$multiplier))

    ruv4out <- ruv::RUV4(Y = Y, X = X[, 2, drop = FALSE], ctl = as.logical(half_null),
                         Z = X[, 1, drop = FALSE], k = num_sv)
    ash_ruv4 <- ashr::ash(c(ruv4out$betahat), sqrt(ruv4out$sigma2 * ruv4out$multiplier))

    tot.time <- proc.time() - start.time

    postmean_df <- data.frame(ruvash_new_u  = ruvash_new_u$PosteriorMean,
                              ruvash_new_s  = ruvash_new_s$PosteriorMean,
                              ruvash_old_un = ruvash_old_un$PosteriorMean,
                              ruvash_old_sn = ruvash_old_sn$PosteriorMean,
                              ruvash_old_ui = ruvash_old_ui$PosteriorMean,
                              ruvash_old_si = ruvash_old_si$PosteriorMean,
                              ash_ruvinv    = ash_ruvinv$PosteriorMean,
                              ash_ruvrinv   = ash_ruvrinv$PosteriorMean,
                              ash_ruv4      = ash_ruv4$PosteriorMean)
    pi0vec <- c(ruvash_new_u$fitted.g$pi[1],
                ruvash_new_s$fitted.g$pi[1],
                ruvash_old_un$fitted.g$pi[1],
                ruvash_old_sn$fitted.g$pi[1],
                ruvash_old_ui$fitted.g$pi[1],
                ruvash_old_si$fitted.g$pi[1],
                ash_ruvinv$fitted.g$pi[1],
                ash_ruvrinv$fitted.g$pi[1],
                ash_ruv4$fitted.g$pi[1])
    lfdr_df <- data.frame(ruvash_new_u  = ruvash_new_u$lfdr,
                          ruvash_new_s  = ruvash_new_s$lfdr,
                          ruvash_old_un = ruvash_old_un$lfdr,
                          ruvash_old_sn = ruvash_old_sn$lfdr,
                          ruvash_old_ui = ruvash_old_ui$lfdr,
                          ruvash_old_si = ruvash_old_si$lfdr,
                          ash_ruvinv    = ash_ruvinv$lfdr,
                          ash_ruvrinv   = ash_ruvrinv$lfdr,
                          ash_ruv4      = ash_ruv4$lfdr)

    nmeth <- length(pi0vec)

    mse <- colMeans((postmean_df - beta_true) ^ 2)
    auc <- rep(NA, nmeth)
    if(args_val$nullpi < 1) {
        for (index in 1:nmeth) {
            auc[index] <- pROC::roc(predictor = lfdr_df[, index], response = which_null)$auc
        }
    }

    scale_val <- c(ruvash_new_u$ruv4$multiplier,
                   ruvash_new_s$ruv4$multiplier,
                   ruvash_old_un$ruv$multiplier,
                   ruvash_old_sn$ruv$multiplier,
                   ruvash_old_ui$ruv$multiplier,
                   ruvash_old_si$ruv$multiplier)

    return_vec <- c(pi0vec, mse, auc, scale_val)

    names(return_vec) <- c("pi0_ruvash_new_u", "pi0_ruvash_new_s",
                           "pi0_ruvash_old_un", "pi0_ruvash_old_sn",
                           "pi0_ruvash_old_ui", "pi0_ruvash_old_si",
                           "pi0_ash_ruvinv", "pi0_ash_ruvrinv",
                           "pi0_ash_ruv4",
                           "mse_ruvash_new_u", "mse_ruvash_new_s",
                           "mse_ruvash_old_un", "mse_ruvash_old_sn",
                           "mse_ruvash_old_ui", "mse_ruvash_old_si",
                           "mse_ash_ruvinv", "mse_ash_ruvrinv",
                           "mse_ash_ruv4",
                           "auc_ruvash_new_u", "auc_ruvash_new_s",
                           "auc_ruvash_old_un", "auc_ruvash_old_sn",
                           "auc_ruvash_old_ui", "auc_ruvash_old_si",
                           "auc_ash_ruvinv", "auc_ash_ruvrinv",
                           "auc_ash_ruv4",
                           "scale_ruvash_new_u", "scale_ruvash_new_s",
                           "scale_ruvash_old_un", "scale_ruvash_old_sn",
                           "scale_ruvash_old_ui", "scale_ruvash_old_si")


    return(return_vec)
}

itermax <- 100

## these change
Nsamp_seq  <- c(5, 10, 20)
nullpi_seq <- c(0.5, 0.9, 1)

par_vals <- expand.grid(list(1:itermax, Nsamp_seq, nullpi_seq))
colnames(par_vals) <- c("current_seed", "Nsamp", "nullpi")
par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

par_list <- list()
for (list_index in 1:nrow(par_vals)) {
    par_list[[list_index]] <- list()
    for (inner_list_index in 1:ncol(par_vals)) {
        par_list[[list_index]][[inner_list_index]] <- par_vals[list_index, inner_list_index]
        names(par_list[[list_index]])[inner_list_index] <- colnames(par_vals)[inner_list_index]
    }
}

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 1
args_val$tissue       <- "muscle"
args_val$path         <- "../../../data/gtex_tissue_gene_reads/"
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores()-1)
sout <- t(parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val))
stopCluster(cl)

## ## on RCC, use this
## np <- mpi.universe.size() - 1
## cluster <- makeMPIcluster(np)
## sout <- t(snow::parSapply(cl = cluster, X = par_list, FUN = one_rep, current_params = args_val))
## stopCluster(cluster)
## mpi.exit()


save(sout, file = "ruvash_new.Rd")




colnames(sout) <- rep(c("ruvash_adinf_TF", "ruvash_adinf_TT",
                        "ruvash_adinf_FF", "ruvash_adinf_FT"), times = 5)

pi0_mat  <- cbind(par_vals, sout[, 1:9])
mse_mat  <- cbind(par_vals, sout[, 10:18])
auc_mat  <- cbind(par_vals, sout[, 19:27])
scale_mat <- cbind(par_vals, sout[, 28:33])

write.csv(pi0_mat, file = "pi0_ruvash_adinf.csv", row.names = FALSE)
write.csv(mse_mat, file = "mse_ruvash_adinf.csv", row.names = FALSE)
write.csv(auc_mat, file = "auc_ruvash_adinf.csv", row.names = FALSE)
write.csv(scale_mat, file = "scale_ruvash_adinf.csv", row.names = FALSE)
