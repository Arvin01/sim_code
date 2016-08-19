library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
    source("../code/data_generators.R")
    source("../code/adjustment_methods.R")
    source("../code/summary_stat_methods.R")
    source("../code/evaluations.R")
    args_val <- append(current_params, new_params)
    set.seed(new_params$current_seed)
    d_out <- datamaker_counts_only(args_val)
    which_null <- d_out$meta$null
    control_genes <- as.logical(which_null)
    nnull         <- sum(control_genes)
    control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrol)] <- FALSE

    beta_true <- rep(0, length = args_val$Ngene)
    beta_true[!which_null] <- d_out$meta$true_log2foldchange

    X <- as.matrix(model.matrix(~d_out$input$condition))
    colnames(X) <- c("Intercept", "Treatment")
    Y <- t(log2(as.matrix(d_out$input$counts + 1)))
    num_sv <- sva::num.sv(t(Y), mod = X, method = "be")

    start.time <- proc.time()
    method_list               <- list()
    method_list$ols           <- ols(Y = Y, X = X)
    method_list$ruv2          <- ruv2(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
    method_list$ruv4          <- ruv4(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
    method_list$ruvinv        <- ruvinv(Y = Y, X = X, control_genes = control_genes)
    method_list$sva           <- sva(Y = Y, X = X, num_sv = num_sv)
    method_list$cate_nc_nocal <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                         control_genes = control_genes, calibrate = FALSE)
    method_list$cate_nc_cal   <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes, calibrate = TRUE)
    method_list$cate_rr_nocal <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = FALSE)
    method_list$cate_rr_cal   <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = TRUE)
    method_list$vruv4         <- vruv4(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes, adjust_bias = FALSE)
    method_list$vruv2         <- vruv2(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
    method_list$vruvinv       <- vruvinv(Y = Y, X = X, control_genes = control_genes)

    ash_wrap <- function(args) {
        new_args <- list()
        new_args$betahat <- args$betahat
        new_args$sebetahat <- args$sebetahat
        new_args$df <- args$df
        do.call(what = fit_ash, args = new_args)
    }
    ash_list <- lapply(method_list[-c(7, 9)], ash_wrap)

    qvalue_wrap <- function(args) {
        fit_qvalue(args$pvalues)
    }
    qvalue_list <- lapply(method_list, qvalue_wrap)


    auc_wrap <- function(args, which_null) {
        get_auc(lfdr = args$lfdr, which_null = which_null)
    }

    fdp_wrap <- function(args, which_null, fdr_level) {
        get_fdr(args$pvalues, which_null, fdr_level)
    }

    qvalue_auc    <- sapply(qvalue_list, auc_wrap, which_null = which_null)
    ash_auc       <- sapply(ash_list, auc_wrap, which_null = which_null)
    ash_pi0hat    <- sapply(ash_list, FUN = function(x){x$pi0hat})
    qvalue_pi0hat <- sapply(qvalue_list, FUN = function(x){x$pi0hat})
    fdp_vec       <- sapply(method_list, fdp_wrap, which_null = as.logical(which_null),
                            fdr_level = 0.1)


    return_vec <- c(qvalue_auc, ash_auc, ash_pi0hat, qvalue_pi0hat, fdp_vec)

    tot.time <- proc.time() - start.time
    return(return_vec)
}

itermax <- 100

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(5, 10, 20)
ncontrol_seq <- c(100, 1000)

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
args_val$tissue       <- "muscle"
args_val$path         <- "../../../data/gtex_tissue_gene_reads/"
args_val$Ngene        <- 10000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

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


auc_pvalues   <- cbind(par_vals, sout[, 1:12])
auc_ash       <- cbind(par_vals, sout[, 13:22])
pi0hat_ash    <- cbind(par_vals, sout[, 23:32])
pi0hat_qvalue <- cbind(par_vals, sout[, 33:44])
fdp_mat       <- cbind(par_vals, sout[, 45:56])

write.csv(auc_pvalues, file = "auc_pvalues.csv", row.names = FALSE)
write.csv(auc_ash, file = "auc_ash.csv", row.names = FALSE)
write.csv(pi0hat_ash, file = "pi0hat_ash.csv", row.names = FALSE)
write.csv(pi0hat_qvalue, file = "pi0hat_qvalue.csv", row.names = FALSE)
write.csv(fdp_mat, file = "fdp_mat.csv", row.names = FALSE)
