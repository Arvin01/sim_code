source("data_generators/datamaker_only_counts.R")

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 1
args_val$tissue       <- "muscle"
args_val$path         <- "../../../data/gtex_tissue_gene_reads/"
args_val$Ngene        <- 10000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0
args_val$Nsamp        <- Nsamp
args_val$nullpi       <- nullpi

if (nullpi != 1) {
    args_val$poisthin <- TRUE
}

d_out                  <- datamaker_counts_only(args_val)
which_null             <- d_out$meta$null
control_genes          <- as.logical(which_null)
nnull                  <- sum(control_genes)
control_genes[control_genes][sample(1:nnull, size = nnull - ncontrol)] <- FALSE
beta_true              <- rep(0, length = args_val$Ngene)
if (nullpi != 1) {
    beta_true[!which_null] <- d_out$meta$true_log2foldchange
}
X                      <- as.matrix(model.matrix(~d_out$input$condition))
colnames(X)            <- c("Intercept", "Treatment")
Y                      <- t(log2(as.matrix(d_out$input$counts + 1)))
num_sv                 <- sva::num.sv(t(Y), mod = X, method = "be")
