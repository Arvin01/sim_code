#########################
## Synopsis: plot simulations from ruv3paper_sims.R
########################
library(reshape2)
library(ggplot2)
library(dplyr)

mse_mat <- read.csv("mse_mat8.csv")
auc_mat <- read.csv("auc_mat8.csv")
cov_mat <- read.csv("cov_mat8.csv")


## Coverage Plot -------------------------------------------------
longdat <- melt(data = cov_mat, id.vars = 1:5)
ggplot(data = longdat, mapping = aes(y = value, x = variable)) +
    geom_boxplot() +
    facet_grid(Nsamp ~ nullpi + ncontrols) +
    geom_hline(yintercept = 0.95) +
    theme_bw()



## AUC plot ------------------------------------------------------
nullpi_seq <- unique(auc_mat$nullpi)
nsamp_seq <- unique(auc_mat$Nsamp)
dummy_dat <- expand.grid(nullpi_seq, nsamp_seq)
colnames(dummy_dat) <- c("nullpi", "Nsamp")
med_mat <- matrix(NA, nrow = length(nullpi_seq) * length(nsamp_seq), ncol = ncol(auc_mat) - 5)
for (index in 6:ncol(auc_mat)) {
    form1 <- as.formula(paste(colnames(auc_mat)[index], "~ nullpi + Nsamp"))
    out1 <- aggregate(form1, FUN = median, na.rm = TRUE,
                      data = auc_mat)
    med_mat[, index - 5] <- out1[, 3]
}
dummy_dat <- cbind(expand.grid(nullpi_seq[nullpi_seq != 1], nsamp_seq), apply(med_mat, 1, max))
colnames(dummy_dat) <- c("nullpi", "Nsamp", "max_med")

longdat <- filter(melt(data = auc_mat, id.vars = 1:5), nullpi != 1)
ggplot(data = longdat, mapping = aes(y = value, x = variable)) +
    geom_boxplot() +
    facet_grid(Nsamp ~ nullpi + ncontrols) +
    geom_hline(data = dummy_dat, mapping = aes(yintercept = max_med)) +
    theme_bw()

temp <- filter(auc_mat, Nsamp == 5, nullpi == 0.5)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
temp <- filter(auc_mat, Nsamp == 5, nullpi == 0.9)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
temp <- filter(auc_mat, Nsamp == 10, nullpi == 0.5)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
temp <- filter(auc_mat, Nsamp == 10, nullpi == 0.9)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)


longdat <- melt(data = mse_mat, id.vars = 1:5)
ggplot(data = longdat, mapping = aes(y = value, x = variable)) +
    geom_boxplot() +
    facet_grid(Nsamp ~ nullpi + ncontrols) +
    ylim(c(0, 0.5)) +
    theme_bw()

## Linked FA ------------------------------------------------------------
mse_mat <- read.csv("mse_mat.csv")
auc_mat <- read.csv("auc_mat.csv")
cov_mat <- read.csv("cov_mat.csv")

temp <- filter(auc_mat, Nsamp == 5, nullpi == 0.5, ncontrols == 10)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
temp <- filter(auc_mat, Nsamp == 5, nullpi == 0.5, ncontrols == 100)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)

temp <- filter(auc_mat, Nsamp == 5, nullpi == 0.9, ncontrols == 10)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
temp <- filter(auc_mat, Nsamp == 5, nullpi == 0.9, ncontrols == 100)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)

temp <- filter(auc_mat, Nsamp == 10, nullpi == 0.5, ncontrols == 10)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
temp <- filter(auc_mat, Nsamp == 10, nullpi == 0.5, ncontrols == 100)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)


temp <- filter(auc_mat, Nsamp == 10, nullpi == 0.9, ncontrols == 10)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
temp <- filter(auc_mat, Nsamp == 10, nullpi == 0.9, ncontrols == 100)
t.test(temp$ruvb, temp$ruv2, paired = TRUE)
