#########################
## Synopsis: plot simulations from ruv3paper_sims.R
########################
library(reshape2)
library(ggplot2)
library(dplyr)

mse_mat <- read.csv("mse_mat2.csv")
auc_mat <- read.csv("auc_mat2.csv")
cov_mat <- read.csv("cov_mat2.csv")

colnames(cov_mat)[6:ncol(cov_mat)] <- c("OLS", "RUV2", "RUV2c", "RUV3",
                                        "RUV4", "RUV4c", "CATE", "RUVB")

colnames(auc_mat)[6:ncol(auc_mat)] <- c("OLS", "RUV2", "RUV2c", "RUV3",
                                        "RUV4", "RUV4c", "CATE", "RUVB")
## Coverage Plot -------------------------------------------------

longdat <- melt(data = cov_mat, id.vars = 1:5)
longdat$Nsamp <- longdat$Nsamp * 2
pdf(file = "coverage.pdf", height = 8, width = 6.5, family = "Times")
p <- ggplot(data = longdat, mapping = aes(y = value, x = variable)) +
    geom_boxplot(outlier.size = 0.5, size = 0.5) +
    facet_grid(nullpi + ncontrols ~ Nsamp) +
    geom_hline(yintercept = 0.95, lty = 2) +
    xlab("Method") + ylab("Coverage") + ylim(0.57, 1) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(strip.background = element_rect(fill="white"))
print(p)
dev.off()



## AUC plot ------------------------------------------------------
nullpi_seq <- unique(auc_mat$nullpi)


nsamp_seq <- unique(auc_mat$Nsamp)
ncontrol_seq <- unique(auc_mat$ncontrols)
dummy_dat <- expand.grid(nullpi_seq, nsamp_seq, ncontrol_seq)
colnames(dummy_dat) <- c("nullpi", "Nsamp", "ncontrols")
med_mat <- matrix(NA, nrow = (length(nullpi_seq) - 1) * length(nsamp_seq) * length(ncontrol_seq),
                  ncol = ncol(auc_mat) - 5)
for (index in 6:ncol(auc_mat)) {
    form1 <- as.formula(paste(colnames(auc_mat)[index], "~ nullpi + Nsamp + ncontrols"))
    out1 <- aggregate(form1, FUN = median, na.rm = TRUE,
                      data = auc_mat)
    med_mat[, index - 5] <- out1[, 4]
}
dummy_dat <- cbind(expand.grid(nullpi_seq[nullpi_seq != 1], nsamp_seq, ncontrol_seq),
                   apply(med_mat, 1, max))
colnames(dummy_dat) <- c("nullpi", "Nsamp", "ncontrols", "max_med")

longdat <- filter(melt(data = auc_mat, id.vars = 1:5), nullpi != 1)
longdat$Nsamp <- longdat$Nsamp * 2
dummy_dat$Nsamp <- dummy_dat$Nsamp * 2
pdf(file = "auc.pdf", height = 8, width = 6.5, family = "Times")
p <- ggplot(data = longdat, mapping = aes(y = value, x = variable)) +
    geom_boxplot(outlier.size = 0.5, size = 0.5) +
    facet_grid(nullpi + ncontrols ~ Nsamp) +
    geom_hline(data = dummy_dat, mapping = aes(yintercept = max_med), lty = 2) +
    xlab("Method") + ylab("AUC") +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(strip.background = element_rect(fill="white"))
print(p)
dev.off()


auc_mat$bm2 <- auc_mat$RUVB - auc_mat$RUV2
auc_mat$threem2 <-auc_mat$RUV3 - auc_mat$RUV2
auc_mat$bmc <- auc_mat$RUVB - auc_mat$CATE
auc_mat$bm3 <- auc_mat$RUVB - auc_mat$RUV3


temp <- filter(auc_mat, nullpi != 1)

pdf(file = "bm2.pdf", family = "Times", height = 3, width = 3)
p1 <- ggplot(data = temp, mapping = aes(y = bm2, x = as.factor(Nsamp))) +
    facet_grid(nullpi ~ ncontrols) +
    geom_boxplot(outlier.size = 0.3, size = 0.2) +
    geom_hline(yintercept = 0, lty = 2) +
    xlab("Sample Size") + ylab("RUVB - RUV2") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"))
print(p1)
dev.off()

pdf(file = "tm2.pdf", family = "Times", height = 3, width = 3)
p2 <- ggplot(data = temp, mapping = aes(y = threem2, x = as.factor(Nsamp))) +
    facet_grid(nullpi ~ ncontrols) +
    geom_boxplot(outlier.size = 0.2, size = 0.2) +
    geom_hline(yintercept = 0, lty = 2) +
    xlab("Sample Size") + ylab("RUV3 - RUV2") +
    theme_bw() +
    theme(strip.background = element_rect(fill="white"))
print(p2)
dev.off()


ggplot(data = temp, mapping = aes(y = bm3, x = as.factor(Nsamp))) +
    facet_grid(nullpi ~ ncontrols) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty = 2) +
    xlab("Sample Size") + ylab("RUVB - RUV2")

ggplot(data = temp, mapping = aes(y = bmc, x = as.factor(Nsamp))) +
    facet_grid(nullpi ~ ncontrols) +
    geom_boxplot() +
    geom_hline(yintercept = 0, lty = 2) +
    xlab("Sample Size") + ylab("RUVB - RUV2")
