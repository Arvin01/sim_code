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
    succ_nullmle <- succotashr::succotash(Y = Y, X = X, k = num_sv,
                                          two_step = TRUE,
                                          z_init_type = "null_mle",
                                          var_scale_init_type = "null_mle",
                                          optmethod = "em")
    succ_random <- succotashr::succotash(Y = Y, X = X, k = num_sv,
                                         two_step = TRUE,
                                         z_init_type = "random",
                                         var_scale_init_type = "random",
                                         optmethod = "em")
    tot.time <- proc.time() - start.time

    postmean_df <- data.frame(succ_nullmle = succ_nullmle$betahat,
                              succ_random = succ_random$betahat)
    pi0vec      <- c(succ_nullmle = succ_nullmle$pi0,
                     succ_random = succ_random$pi0)
    scale_val_vec <-  c(succ_nullmle = succ_nullmle$scale_val,
                        succ_random = succ_random$scale_val)
    llike_vec <-  c(succ_nullmle = succ_nullmle$llike,
                        succ_random = succ_random$llike)
    null_llike_vec <-  c(succ_nullmle = succ_nullmle$null_llike,
                        succ_random = succ_random$null_llike)
    lfdr_df     <- data.frame(succ_nullmle = succ_nullmle$lfdr,
                              succ_random = succ_random$lfdr)
    nmeth <- length(pi0vec)

    mse <- colMeans((postmean_df - beta_true) ^ 2)
    auc <- rep(NA, nmeth)
    if(args_val$nullpi < 1) {
        for (index in 1:nmeth) {
            auc[index] <- pROC::roc(predictor = lfdr_df[, index], response = which_null)$auc
        }
    }

    return_vec <- c(pi0vec, mse, auc, scale_val_vec, llike_vec, null_llike_vec)

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


save(sout, file = "sout_succ_diff_init.Rd")

colnames(sout) <- rep(c("succ_null_mle", "succ_random"), times = 6)

pi0_mat <- cbind(par_vals, sout[, 1:2])
mse_mat <- cbind(par_vals, sout[, 3:4])
auc_mat <- cbind(par_vals, sout[, 5:6])
scale_val_mat <- cbind(par_vals, sout[, 7:8])
llike_mat <- cbind(par_vals, sout[, 9:10])
null_llike_mat <- cbind(par_vals, sout[, 11:12])

write.csv(pi0_mat, file = "pi0_ruvash_alpha1.csv", row.names = FALSE)
write.csv(mse_mat, file = "mse_ruvash_alpha1.csv", row.names = FALSE)
write.csv(auc_mat, file = "auc_ruvash_alpha1.csv", row.names = FALSE)
write.csv(scale_val_mat, file = "scale_val_ruvash_alpha1.csv", row.names = FALSE)
write.csv(llike_mat, file = "llike_ruvash_alpha1.csv", row.names = FALSE)
write.csv(null_llike_mat, file = "null_llike_ruvash_alpha1.csv", row.names = FALSE)
