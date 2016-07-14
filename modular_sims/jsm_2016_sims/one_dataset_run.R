library(ggplot2)
library(reshape2)
source("../code/data_generators.R")
source("../code/adjustment_methods.R")
source("../code/summary_stat_methods.R")
source("../code/evaluations.R")

set.seed(117)
dout <- pois_thin(Nsamp = 10, nullpi = 1, path = "../../../data/gtex_tissue_gene_reads/",
                  ncontrol = 100, Ngene = 10000)
Y <- dout$Y
X <- dout$X
num_sv <- dout$num_sv
control_genes <- dout$control_genes
method_list               <- list()
method_list$ols           <- ols(Y = Y, X = X)
method_list$ruv2          <- ruv2(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
method_list$ruv4          <- ruv4(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
method_list$ruvinv        <- ruvinv(Y = Y, X = X, control_genes = control_genes)
method_list$vruv4         <- vruv4(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = control_genes)
method_list$sva           <- sva(Y = Y, X = X, num_sv = num_sv)
method_list$cate_nc_nocal <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes, calibrate = FALSE)
method_list$cate_nc_cal   <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes, calibrate = TRUE)
method_list$cate_rr_nocal <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = FALSE)
method_list$cate_rr_cal   <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = TRUE)



## pos <- 888
## tout <- t.test(Y[as.logical(X[, 2]), pos], Y[!as.logical(X[, 2]), pos])
## tout$p.value
## method_list$ols$pvalues[pos]

pdf(file = "ols_p.pdf")
qplot(method_list$ols$pvalues, bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("T-test P-values") +
    ylab("Counts") +
    theme(text = element_text(size = 20))
dev.off()

sum(p.adjust(p = method_list$ols$pvalues, method = "BH") < 0.1)

qplot(method_list$ruv4$pvalues, bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("RUV4 P-values") +
    ylab("Count")

qplot(method_list$cate_nc_cal$pvalues, bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Cate_nc_cal P-values") +
    ylab("Count")

qplot(method_list$cate_nc_nocal$pvalues, bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Cate_nc_nocal P-values") +
    ylab("Count")

pdf(file = "sva_p.pdf")
qplot(method_list$sva$pvalues, bins = 30, col = I("black"), fill = I("#FABB69")) +
    theme_bw() +
    xlab("SVA P-values") +
    ylab("Counts") +
    theme(text = element_text(size = 20))
dev.off()

qplot(method_list$vruv4$pvalues, bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Vicar P-values") +
    ylab("Count")


pmat <- cbind(method_list$ols$pvalues,
              method_list$sva$pvalues,
              method_list$ruv2$pvalues,
              method_list$ruv4$pvalues,
              method_list$ruvinv$pvalues,
              method_list$cate_nc_nocal$pvalues,
              method_list$cate_nc_cal$pvalues,
              method_list$cate_rr_nocal$pvalues,
              method_list$cate_rr_cal$pvalues)

colnames(pmat) <- c("OLS", "SVA", "RUV2", "RUV4", "RUVinv", "CATEnc", "Cal CATEnc", "CATErr",
                    "Cal CATErr")

pdf(file = "all_null.pdf", height = 5, width = 9)
longdat <- melt(pmat, id.vars = NULL)
colnames(longdat) <- c("Var1", "Method", "P-values")
ggplot(data = longdat, mapping = aes(x = `P-values`, fill = Method, color = I("black"))) +
    geom_histogram(bins = 30) +
    facet_wrap(~ Method) +
    theme_bw()  +
    theme(text = element_text(size = 20)) +
    ylab("Counts")
dev.off()




pmat <- cbind(method_list$ruv2$pvalues,
              method_list$ruv4$pvalues,
              method_list$ruvinv$pvalues)
colnames(pmat) <- c("RUV2", "RUV4", "RUVinv")

pdf(file = "ruv.pdf", height = 6, width = 9)
longdat <- melt(pmat, id.vars = NULL)
colnames(longdat) <- c("Var1", "Method", "P-values")
ggplot(data = longdat, mapping = aes(x = `P-values`, fill = Method, color = I("black"))) +
    geom_histogram(bins = 30) +
    facet_wrap(~ Method) +
    theme_bw()  +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    ylab("Counts")
dev.off()



pmat <- cbind(method_list$cate_nc_nocal$pvalues,
              method_list$cate_nc_cal$pvalues,
              method_list$cate_rr_nocal$pvalues,
              method_list$cate_rr_cal$pvalues)

colnames(pmat) <- c("CATEnc", "Calibrated CATEnc", "CATErr",
                    "Calibrated CATErr")

pdf(file = "cate.pdf", height = 6.5, width = 10)
longdat <- melt(pmat, id.vars = NULL)
colnames(longdat) <- c("Var1", "Method", "P-values")
ggplot(data = longdat, mapping = aes(x = `P-values`, fill = Method, color = I("black"))) +
    geom_histogram(bins = 30) +
    facet_grid(.~ Method) +
    theme_bw()  +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 20, angle = 90),
          legend.position = "none") +
    ylab("Counts")
dev.off()



### Same as above but nullpi = 0.5 -----------------------------------------------------------
library(ggplot2)
source("../code/data_generators.R")
source("../code/adjustment_methods.R")
source("../code/summary_stat_methods.R")
source("../code/evaluations.R")


set.seed(121)
dout <- pois_thin(Nsamp = 10, nullpi = 0.5, path = "../../../data/gtex_tissue_gene_reads/",
                  ncontrol = 500, Ngene = 10000)
Y <- dout$Y
X <- dout$X
num_sv <- dout$num_sv
control_genes <- dout$control_genes
which_null <- as.logical(dout$which_null)
method_list               <- list()
method_list$ols           <- ols(Y = Y, X = X)
method_list$ruv2          <- ruv2(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
method_list$ruv4          <- ruv4(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
method_list$ruvinv        <- ruvinv(Y = Y, X = X, control_genes = control_genes)
method_list$vruv4         <- vruv4(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = control_genes)
method_list$sva           <- sva(Y = Y, X = X, num_sv = num_sv)
method_list$cate_nc_nocal <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes, calibrate = FALSE)
method_list$cate_nc_cal   <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                     control_genes = control_genes, calibrate = TRUE)
method_list$cate_rr_nocal <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = FALSE)
method_list$cate_rr_cal   <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = TRUE)



## pos <- 888
## tout <- t.test(Y[as.logical(X[, 2]), pos], Y[!as.logical(X[, 2]), pos])
## tout$p.value
## method_list$ols$pvalues[pos]

qplot(method_list$ols$pvalues[which_null], bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("OLS Null P-values") +
    ylab("Count")

qplot(method_list$ruv4$pvalues[which_null], bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("RUV4 Null P-values") +
    ylab("Count")

qplot(method_list$cate_nc_cal$pvalues[which_null], bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Calibrated Cate Null P-values") +
    ylab("Count")

qplot(method_list$cate_rr_cal$pvalues[which_null], bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Calibrated Cate Null P-values") +
    ylab("Count")

qplot(method_list$cate_nc_cal$pvalues[!which_null], bins = 30, col = I("black"),
      fill = I("#52EF87")) +
    theme_bw() +
    xlab("Calibrated Cate Non-null P-values") +
    ylab("Count")

qplot(method_list$cate_nc_nocal$pvalues[which_null], bins = 30,
      col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Uncalibrated Cate Null P-values") +
    ylab("Count")

qplot(method_list$sva$pvalues[which_null], bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Cate_nc_nocal Null P-values") +
    ylab("Count")

qplot(method_list$vruv4$pvalues[which_null], bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Vicar Null P-values") +
    ylab("Count")

qplot(method_list$vruv4$pvalues[!which_null], bins = 30, col = I("black"), fill = I("#52EF87")) +
    theme_bw() +
    xlab("Vicar Non-null P-values") +
    ylab("Count")


pmat <- cbind(method_list$cate_nc_nocal$pvalues[which_null],
              method_list$cate_nc_cal$pvalues[which_null],
              method_list$cate_rr_nocal$pvalues[which_null],
              method_list$cate_rr_cal$pvalues[which_null])

colnames(pmat) <- c("CATEnc", "Calibrated CATEnc", "CATErr",
                    "Calibrated CATErr")

pdf(file = "cate_nonnull.pdf", height = 6, width = 9)
longdat <- melt(pmat, id.vars = NULL)
colnames(longdat) <- c("Var1", "Method", "Null P-values")
ggplot(data = longdat, mapping = aes(x = `Null P-values`, fill = Method, color = I("black"))) +
    geom_histogram(bins = 30) +
    facet_grid(.~ Method) +
    theme_bw()  +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(size = 20, angle = 90),
          legend.position = "none") +
    ylab("Counts")
dev.off()
