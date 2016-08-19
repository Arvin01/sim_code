library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
    source("../code/data_generators.R")
    source("../code/adjustment_methods.R")
    source("../code/summary_stat_methods.R")
    source("../code/evaluations.R")
    args_val <- append(current_params, new_params)
    set.seed(new_params$current_seed)
    d_out <- datamaker_change_ncovs(Ngene = args_val$Ngene,
                                    Nsamp = args_val$Nsamp, ncovs = 5,
                                    nullpi = args_val$nullpi,
                                    ncontrols = args_val$ncontrols,
                                    alpha_sd = 2, sigsd = 1/2)

    which_null <- d_out$which_null
    control_genes <- as.logical(which_null)
    nnull         <- sum(control_genes)
    control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrol)] <- FALSE
    nnull         <- sum(control_genes)
    beta_true <- d_out$beta[2, ]

    X <- d_out$X
    Y <- log2(as.matrix(d_out$Y + 1))
    num_sv <- sva::num.sv(t(Y), mod = X, method = "be")

    start.time <- proc.time()
    method_list               <- list()
    method_list$ols           <- lm(Y ~ X)$coefficients
    method_list$ruv2          <- ruv::RUV2(Y = Y, X = X, ctl = control_genes, k = num_sv,
                                           Z = NULL)$betahat
    method_list$vruv4         <- vicar::vruv4(Y = Y, X = X, k = num_sv,
                                              ctl = control_genes, include_intercept = FALSE,
                                              cov_of_interest = 1:ncol(X))$betahat
    method_list$ruv4          <- ruv::RUV4(Y = Y, X = X, ctl = control_genes, k = num_sv,
                                           Z = NULL)$betahat
    method_list$ruv3          <- vicar::ruv3(Y = Y, X = X, k = num_sv, ctl = control_genes,
                                             cov_of_interest = 1:ncol(X))$betahat
    method_list$ruvimpute     <- vicar::ruvimpute(Y = Y, X = X, ctl = control_genes, k = num_sv,
                                                  cov_of_interest = 1:ncol(X))$betahat

    mse_wrap <- function(args, beta_true, ctl) {
        mean((args[!ctl] - beta_true[!ctl]) ^ 2)
    }

    mse_vec <- sapply(method_list, mse_wrap, beta_true = beta_true, ctl = control_genes)

    tot.time <- proc.time() - start.time
    return(mse_vec)
}



itermax <- 100

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(10, 20, 40)
ncontrol_seq <- c(50, 100)

par_vals <- expand.grid(list(1:itermax, nullpi_seq, Nsamp_seq, ncontrol_seq))
colnames(par_vals) <- c("current_seed", "nullpi", "Nsamp", "ncontrols")
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
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0

## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores() - 1)
sout <- t(snow::parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val))
stopCluster(cl)

## ## on RCC, use this
## np <- mpi.universe.size() - 1
## cluster <- makeMPIcluster(np)
## sout <- t(snow::parSapply(cl = cluster, X = par_list, FUN = one_rep, current_params = args_val))
## stopCluster(cluster)
## mpi.exit()


save(sout, file = "general_sims.Rd")

qvalue_auc_mat <- cbind(par_vals, sout[, 1:6])
ash_auc_mat <- cbind(par_vals, sout[, 7:12])
ash_pi0hat_mat <- cbind(par_vals, sout[, 13:18])
qvalue_pi0hat_mat <- cbind(par_vals, sout[, 19:24])
mse_mat <- cbind(par_vals, sout[, 25:31])

write.csv(qvalue_auc_mat, file = "qvalue_auc_mat.csv", row.names = FALSE)
write.csv(ash_auc_mat, file = "ash_auc_mat.csv", row.names = FALSE)
write.csv(ash_pi0hat_mat, file = "ash_pi0hat_mat.csv", row.names = FALSE)
write.csv(qvalue_pi0hat_mat, file = "qvalue_pi0hat_mat.csv", row.names = FALSE)
write.csv(mse_mat, file = "mse_mat.csv", row.names = FALSE)
