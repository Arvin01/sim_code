library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
    return_vec <- tryCatch(expr = {
        source("../code/data_generators.R")
        source("../code/adjustment_methods.R")
        source("../code/summary_stat_methods.R")
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

        counts_adjusted <- d_out$input$counts
        Y <- t(log2(as.matrix(counts_adjusted + 1)))


        num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

        method_list            <- list()
        method_list$ols        <- ols(Y = Y, X = X)
        method_list$ruv2       <- ruv2(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
        method_list$ruv2_rsvar <- ruv2(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
        method_list$ruv3       <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes,
                                       multiplier = FALSE)
        method_list$ruv4       <- ruv4(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
        method_list$ruv4_rsvar <- ruv4_rsvar_ebayes(Y = Y, X = X, num_sv = num_sv,
                                                    control_genes = control_genes)
        method_list$catenc     <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                          control_genes = control_genes,
                                          calibrate = TRUE)
        method_list$caterr     <- cate_rr(Y = Y, X = X, num_sv = num_sv)
        method_list$ruv4v      <- vruv4(Y = Y, X = X, num_sv = num_sv,
                                        control_genes = control_genes)
        method_list$sva        <- sva(Y = Y, X = X, num_sv = num_sv)

        method_names <- names(method_list)


        method_list_ash <- list()
        for (list_index in 1:length(method_list)) {
            method_list_ash[[list_index]] <- fit_ash(c(method_list[[list_index]]$betahat),
                                                     c(method_list[[list_index]]$sebetahat),
                                                     method_list[[list_index]]$df[1])
        }

        method_list_qvalue <- list()
        for (list_index in 1:length(method_list)) {
            method_list_qvalue[[list_index]] <- fit_qvalue(method_list[[list_index]]$pvalues)
        }




        get_mse <- function(args, beta_true) {
            mean((args$betahat - beta_true) ^ 2)
        }

        get_auc_pvalues <- function(args, which_null) {
            if (sum(which_null) == length(which_null)) {
                return(NA)
            }
            pROC::roc(predictor = c(args$pvalues), response = which_null)$auc
        }

        get_auc_lfdr <- function(args, which_null) {
            if (sum(which_null) == length(which_null)) {
                return(NA)
            }
            pROC::roc(predictor = c(args$lfdr), response = which_null)$auc
        }

        get_pi0hat <- function(args) {
            args$pi0hat
        }

        auc_ash <- sapply(method_list_ash, FUN = get_auc_lfdr, which_null = which_null)
        auc_pvalues <- sapply(method_list, FUN = get_auc_pvalues, which_null = which_null)
        names(auc_ash)     <- paste0("auc_ash_", method_names)
        names(auc_pvalues) <- paste0("auc_pvalues_", method_names)

        pi0hat_ash <- sapply(method_list_ash, get_pi0hat)
        pi0hat_qvalue <- sapply(method_list_qvalue, get_pi0hat)
        names(pi0hat_ash)     <- paste0("pi0_ash_", method_names)
        names(pi0hat_qvalue) <- paste0("pi0_qvalue_", method_names)

        mse_ash <- sapply(method_list_ash, get_mse, beta_true = beta_true)
        mse_default <- sapply(method_list, get_mse, beta_true = beta_true)
        names(mse_ash)     <- paste0("mse_ash_", method_names)
        names(mse_default) <- paste0("mse_default_", method_names)

        ## SUCCOTASH results ----------------------------------------------------
        mouth_noscale  <- mouthwash_noscale(Y = Y, X = X, num_sv = num_sv)
        mouth_scale    <- mouthwash_scale(Y = Y, X = X, num_sv = num_sv)
        mouthn_noscale <- mouthwash_normal_noscale(Y = Y, X = X, num_sv = num_sv)
        mouthn_scale   <- mouthwash_normal_scale(Y = Y, X = X, num_sv = num_sv)
        mouth_list     <- list(mouth_noscale, mouth_scale, mouthn_noscale, mouthn_scale)

        mse_mouth <- sapply(mouth_list, get_mse, beta_true = beta_true)
        names(mse_mouth) <- c("mse_mouthTns", "mse_mouthTs", "mse_mouthNns", "mse_mouthNs")

        auc_mouth <- sapply(mouth_list, get_auc_lfdr, which_null = which_null)
        names(auc_mouth) <- c("auc_mouthTns", "auc_mouthTs", "auc_mouthNns", "auc_mouthNs")

        pi0hat_mouth <- sapply(mouth_list, get_pi0hat)
        names(pi0hat_mouth) <- c("pi0_mouthTns", "pi0_mouthTs", "pi0_mouthNns", "pi0_mouthNs")

        ## set return values -----------------------------------------------------

        return_vec <- c(mse_default, mse_ash, mse_mouth,
                        auc_pvalues, auc_ash, auc_mouth,
                        pi0hat_qvalue, pi0hat_ash, pi0hat_mouth)

        return(return_vec)

    }, error = function(e){rep(NA, 72)})
    return(return_vec)
}

itermax <- 100
seed_start <- 436

## these change
nullpi_seq   <- c(0.5, 1)
Nsamp_seq    <- c(10, 20)
ncontrol_seq <- c(100)

par_vals <- expand.grid(list((1 + seed_start):(itermax + seed_start),
                             nullpi_seq, Nsamp_seq, ncontrol_seq))
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
args_val$log2foldsd   <- 0.8
args_val$tissue       <- "skin"
args_val$path         <- "../../../data/gtex_tissue_gene_reads_v6p/"
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

one_rep(par_list[[1]], args_val)

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


save(sout, file = "general_sims2.Rd")





mse_mat <- cbind(par_vals, sout[, 1:24])
auc_mat <- cbind(par_vals, sout[, 25:48])
pi0_mat <- cbind(par_vals, sout[, 49:72])
write.csv(mse_mat, file = "mse_mat2.csv", row.names = FALSE)
write.csv(auc_mat, file = "auc_mat2.csv", row.names = FALSE)
write.csv(pi0_mat, file = "pi0_mat2.csv", row.names = FALSE)





## from par_list[[1070]], chosen by a random seed
library(coda)
bout <- vicar::ruvb(Y = Y, X = X, k = num_sv, ctl = control_genes,
                    return_mcmc = TRUE)

mcmc_b <- mcmc(t(bout$betahat_post[1, , ,drop = TRUE]))
gout <- geweke.diag(mcmc_b)
qqnorm(gout$z)
abline(0, 1)

traceplot(mcmc_b)

library(ggplot2)
eout <- effectiveSize(mcmc_b)
qplot(eout, bins = 20, fill = I("white"), color = I("black")) + theme_bw()

out$betahat_post
