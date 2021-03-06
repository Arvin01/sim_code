library(Rmpi)
library(snow)

one_rep <- function(new_params, current_params) {
    return_vec <- tryCatch(expr = {
        source("../code/data_generators.R")
        source("../code/adjustment_methods.R")
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


        num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

        rmax <- min(sum(control_genes), nrow(X) - ncol(X))
        ## num_sv <- max(cate::est.confounder.num(~ Treatment, X.data = as.data.frame(X),
        ##                                        Y = Y, rmax = rmax,
        ##                                        nRepeat = 100, bcv.plot = FALSE)$r, 1)


        start.time <- proc.time()
        method_list            <- list()
        method_list$ols        <- ols(Y = Y, X = X)
        method_list$ruv2       <- ruv2(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
        method_list$ruv2_rsvar <- ruv2(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
        method_list$ruv3       <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes,
                                       multiplier = FALSE)
        ## method_list$ruv3_rsvar <- ruv3(Y = Y, X = X, num_sv = num_sv,
        ##                                control_genes = control_genes,
        ##                                multiplier = TRUE)
        method_list$ruv4       <- ruv4(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
        method_list$ruv4_rsvar <- ruv4_rsvar_ebayes(Y = Y, X = X, num_sv = num_sv,
                                                    control_genes = control_genes)
        ## method_list$ruvem   <- ruvem(Y = Y, X = X, num_sv = num_sv,
        ##                              control_genes = control_genes)
        method_list$catenc     <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                          control_genes = control_genes,
                                          calibrate = TRUE)
        method_list$ruvb       <- ruvb_bfa_gs(Y = Y, X = X, num_sv = num_sv,
                                              control_genes = control_genes)



        ## false_signs <- sign(beta_true[!control_genes]) != sign(method_list$ruvb$betahat)
        ## sorder <- order(method_list$ruvb$svalues)
        ## qplot(1:(args_val$Ngene - args_val$ncontrols),
        ##       c(method_list$ruvb$svalues)[sorder], col = false_signs[sorder])

        ## fsr <- cumsum(false_signs[sorder]) / (1:(args_val$Ngene - args_val$ncontrols))

        ## plot(c(method_list$ruvb$svalues)[sorder], fsr)
        ## abline(0, 1)

        get_mse <- function(args, beta_true, control_genes) {
            if (length(args$betahat) == length(control_genes)) {
                mean((args$betahat[!control_genes] - beta_true[!control_genes]) ^ 2)
            } else {
                mean((args$betahat - beta_true[!control_genes]) ^ 2)
            }
        }

        get_auc <- function(args, which_null, control_genes) {
            if (sum(which_null) == length(which_null)) {
                return(NA)
            }
            if (length(args$pvalues) == length(control_genes)) {
                pROC::roc(predictor = args$pvalues[!control_genes],
                          response = which_null[!control_genes])$auc
            } else {
                pROC::roc(predictor = c(args$pvalues), response = which_null[!control_genes])$auc
            }
        }

        get_coverage <- function(args, beta_true, control_genes) {
            if(length(args$lower) == length(control_genes)) {
                mean(args$lower[!control_genes] < beta_true[!control_genes] &
                     args$upper[!control_genes] > beta_true[!control_genes])
            } else {
                mean(args$lower < beta_true[!control_genes] &
                     args$upper > beta_true[!control_genes])
            }
        }

        mse_vec <- sapply(method_list, get_mse, beta_true = beta_true,
                          control_genes = control_genes)
        auc_vec <- sapply(method_list, get_auc, which_null = which_null,
                          control_genes = control_genes)
        cov_vec <- sapply(method_list, get_coverage, beta_true = beta_true,
                          control_genes = control_genes)

        return_vec <- c(mse_vec, auc_vec, cov_vec)
        xtot.time <- proc.time() - start.time
        return(return_vec)

    }, error = function(e){rep(NA, 24)})
    return(return_vec)
}

itermax <- 200
seed_start <- 2222

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(3, 5, 10, 20)
ncontrol_seq <- c(10, 100)

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
args_val$tissue       <- "muscle"
args_val$path         <- "../../../data/gtex_tissue_gene_reads/"
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## one_rep(par_list[[21]], args_val)

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
mse_mat <- cbind(par_vals, sout[, 1:8])
auc_mat <- cbind(par_vals, sout[, 9:16])
cov_mat <- cbind(par_vals, sout[, 17:24])
write.csv(mse_mat, file = "mse_mat2.csv", row.names = FALSE)
write.csv(auc_mat, file = "auc_mat2.csv", row.names = FALSE)
write.csv(cov_mat, file = "cov_mat2.csv", row.names = FALSE)
