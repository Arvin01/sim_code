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
    method_list$ruv2          <- ruv2_rsvar_ebayes(Y = Y, X = X, num_sv = num_sv,
                                                   control_genes = control_genes)
    method_list$ruv4          <- ruv4_rsvar_ebayes(Y = Y, X = X, num_sv = num_sv,
                                                   control_genes = control_genes)
    method_list$cate_nc_cal   <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes, calibrate = TRUE)
    method_list$vruv4         <- vruv4(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes, adjust_bias = TRUE)
    method_list$ruv3          <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                      control_genes = control_genes, multiplier = FALSE)
    method_list$ruvimpute     <- ruvimpute(Y = Y, X = X, control_genes = control_genes,
                                           num_sv = num_sv)

    vout <- vicar::vruv4(Y = Y, X = X, ctl = control_genes, k = num_sv)
    hist(vout$sigma2)
    lowrankest <- vout$Zhat %*% t(vout$alphahat)
    sd(c(lowrankest))
    sd(beta_true[!as.logical(which_null)])
    summary(vout$sigma2)

    library(ggplot2)
    qplot(beta_true, method_list$ruv3$betahat) +
        geom_abline(intercept = 0, slope = 1, color = 2, lty = 2)
    qplot(beta_true, c(method_list$ruvimpute$betahat)) +
        geom_abline(intercept = 0, slope = 1, color = 2, lty = 2)
    qplot(beta_true, method_list$ruv4$betahat) +
        geom_abline(intercept = 0, slope = 1, color = 2, lty = 2)
    qplot(beta_true, method_list$ruv2$betahat) +
        geom_abline(intercept = 0, slope = 1, color = 2, lty = 2)
    qplot(beta_true, method_list$ols$betahat) +
        geom_abline(intercept = 0, slope = 1, color = 2, lty = 2)

    ash_wrap <- function(args, ctl) {
        new_args <- list()
        new_args$betahat <- args$betahat[!ctl]
        new_args$sebetahat <- args$sebetahat[!ctl]
        new_args$df <- args$df
        do.call(what = fit_ash, args = new_args)
    }
    ash_list <- lapply(method_list[-length(method_list)], ash_wrap, ctl = control_genes)

    qvalue_wrap <- function(args, ctl) {
        fit_qvalue(args$pvalues[!ctl])
    }
    qvalue_list <- lapply(method_list[-length(method_list)], qvalue_wrap, ctl = control_genes)


    auc_wrap <- function(args, which_null, ctl) {
        get_auc(lfdr = args$lfdr, which_null = which_null[!ctl])
    }

    mse_wrap <- function(args, beta_true, ctl) {
        mean((args$betahat[!ctl] - beta_true[!ctl]) ^ 2)
    }

    qvalue_auc    <- sapply(qvalue_list, auc_wrap, which_null = which_null, ctl = control_genes)
    ash_auc       <- sapply(ash_list, auc_wrap, which_null = which_null, ctl = control_genes)
    ash_pi0hat    <- sapply(ash_list, FUN = function(x){x$pi0hat})
    qvalue_pi0hat <- sapply(qvalue_list, FUN = function(x){x$pi0hat})
    mse_vec       <- sapply(method_list, mse_wrap, beta_true = beta_true, ctl = control_genes)
    mse_vec

    true_pi <- (round(args_val$nullpi * args_val$Ngene) - args_val$ncontrol) /
        (args_val$Ngene - args_val$ncontrol)
    true_pi

    return_vec <- c(qvalue_auc, ash_auc, ash_pi0hat, qvalue_pi0hat, mse_vec)

    tot.time <- proc.time() - start.time
    return(return_vec)
}

itermax <- 100

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(10, 20)
ncontrol_seq <- c(100, 400)

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
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 5

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
