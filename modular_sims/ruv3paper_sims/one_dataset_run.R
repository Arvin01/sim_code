library(ggplot2)
library(reshape2)
library(gridExtra)
source("../code/data_generators.R")
source("../code/adjustment_methods.R")
source("../code/summary_stat_methods.R")
source("../code/evaluations.R")

"#52EF87"
fill_color <- "grey90"
seed_seq <- 117:119
for (index in 1:length(seed_seq)) {
    current_seed <- seed_seq[index]
    set.seed(current_seed)
    dout <- pois_thin(Nsamp = 10, nullpi = 1, path = "../../../data/gtex_tissue_gene_reads/",
                      ncontrol = 300, Ngene = 10000)
    Y <- dout$Y
    X <- dout$X
    olsout <- ols(Y = Y, X = X)

    pmat_temp <- data.frame(pvalues = olsout$pvalues)
    pmat_temp$seed <- current_seed

    if (index == 1) {
        pmat <- pmat_temp
    } else {
        pmat <- rbind(pmat, pmat_temp)
    }
}

pdf(file = "all_null.pdf", height = 2, width = 6, family = "Times", colormodel = "cmyk")
longdat <- melt(pmat[pmat$seed == seed_seq[1], ], id.vars = "seed")
p1 <- ggplot(data = longdat, mapping = aes(x = value, color = I("black"), fill = I(fill_color))) +
    geom_histogram(bins = 20) +
    ylab("Counts") +
    xlab("P-values") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank())

longdat <- melt(pmat[pmat$seed == seed_seq[2], ], id.vars = "seed")
p2 <- ggplot(data = longdat, mapping = aes(x = value, color = I("black"), fill = I(fill_color))) +
    geom_histogram(bins = 20) +
    ylab("") +
    xlab("P-values") +
    theme_bw() +
    theme(strip.background = element_rect(fill = fill_color),
          strip.text.y = element_text(size = 9))

longdat <- melt(pmat[pmat$seed == seed_seq[3], ], id.vars = "seed")
p3 <- ggplot(data = longdat, mapping = aes(x = value, color = I("black"), fill = I(fill_color))) +
    geom_histogram(bins = 20) +
    ylab("") +
    xlab("P-values") +
    theme_bw() +
    theme(strip.background = element_rect(fill = fill_color),
          strip.text.y = element_text(size = 9))

gridExtra::grid.arrange(p1, p2, p3, ncol = 3, nrow = 1, padding = 0)
dev.off()


### Same as above but nullpi = 0.5 -----------------------------------------------------------
rm(list = ls())

source("../code/data_generators.R")
source("../code/adjustment_methods.R")
source("../code/summary_stat_methods.R")
source("../code/evaluations.R")

"#52EF87"
fill_color <- "grey90"

get_p <- function(args) {
    return(args$pvalues)
}

seed_seq <- c(117, 119)
for (index in 1:length(seed_seq)) {
    current_seed <- seed_seq[index]
    set.seed(current_seed)
    dout <- pois_thin(Nsamp = 10, nullpi = 0.5, path = "../../../data/gtex_tissue_gene_reads/",
                      ncontrol = 300, Ngene = 10000)
    Y <- dout$Y
    X <- dout$X
    num_sv <- dout$num_sv
    control_genes <- dout$control_genes
    method_list               <- list()
    method_list$ols           <- ols(Y = Y, X = X)
    ## method_list$sva           <- sva(Y = Y, X = X, num_sv = num_sv)
    method_list$ruv2          <- ruv2(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
    method_list$vruv2         <- vruv2(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
    method_list$ruv4          <- ruv4(Y = Y, X = X, num_sv = num_sv, control_genes = control_genes)
    method_list$cate_nc_nocal <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                         control_genes = control_genes, calibrate = FALSE)
    method_list$cate_nc_cal   <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                         control_genes = control_genes, calibrate = TRUE)
    method_list$vruv4         <- vruv4(Y = Y, X = X, num_sv = num_sv,
                                       control_genes = control_genes)
    method_list$ruvinv        <- ruvinv(Y = Y, X = X, control_genes = control_genes)
    method_list$vruvinv       <- vruvinv(Y = Y, X = X, control_genes = control_genes)
    ## method_list$cate_rr_nocal <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = FALSE)
    ## method_list$cate_rr_cal   <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = TRUE)

    pmat_temp <- as.data.frame(sapply(method_list, get_p))
    pmat_temp$seed <- current_seed
    pmat_temp$which_null <- as.logical(dout$which_null)

    if (index == 1) {
        pmat <- pmat_temp
    } else {
        pmat <- rbind(pmat, pmat_temp)
    }
}


pdf(file = "half_null.pdf", height = 8, width = 6, family = "Times", colormodel = "cmyk")
colnames(pmat) <-  c("OLS", "RUV2", "vRUV2", "RUV4", "CATEnc",
                     "CalCATEnc", "vRUV4", "RUVinv", "vRUVinv",
                     "seed", "which_null")
longdat <- melt(pmat[pmat$seed == seed_seq[1] & pmat$which_null, ],
                id.vars = c("seed", "which_null"))
p1 <- ggplot(data = longdat, mapping = aes(x = value, color = I("black"), fill = I(fill_color))) +
    geom_histogram(bins = 20) +
    facet_grid(variable~., scales = "free") +
    ylab("Counts") +
    xlab("P-values") +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank())

longdat <- melt(pmat[pmat$seed == seed_seq[2] & pmat$which_null, ],
                id.vars = c("seed", "which_null"))
p2 <- ggplot(data = longdat, mapping = aes(x = value, color = I("black"), fill = I(fill_color))) +
    geom_histogram(bins = 20) +
    facet_grid(variable~., scales = "free") +
    ylab("") +
    xlab("P-values") +
    theme_bw() +
    theme(strip.background = element_rect(fill = fill_color),
          strip.text.y = element_text(size = 9))

gridExtra::grid.arrange(p1, p2, ncol = 2, nrow = 1)
dev.off()
