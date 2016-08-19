library(ggplot2)
library(reshape2)

## ASH results ----------------------------------------------------------------
pi0hat_ash <- read.csv("ash_pi0hat_mat.csv")

names(pi0hat_ash) <- c("Seed", "Proportion Null", "Sample Size",
                       "Number Controls", "poisthin", "OLS", "RUV2",
                       "RUVinv", "CATEnc", "vRUV4", "RUV3")

pi0hat_ash$`Sample Size` <- pi0hat_ash$`Sample Size` * 2

longdat <- melt(pi0hat_ash,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(pi0hat_ash)[-(1:5)])


pdf(file = "ash_all.pdf", height = 6, width = 10)
longdat_null <- longdat[longdat$`Proportion Null` == 1, ]
ggplot(longdat_null, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 1, color = "red", lty = 2) +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)

longdat_half <- longdat[longdat$`Proportion Null` == 0.9, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 0.9, color = "red", lty = 2) +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)

longdat_half <- longdat[longdat$`Proportion Null` == 0.5, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 0.5, color = "red", lty = 2) +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)
dev.off()


## qvalue results ---------------------------------------------------------
rm(list = ls())
pi0hat_qv <- read.csv("qvalue_pi0hat_mat.csv")
names(pi0hat_qv) <- c("Seed", "Proportion Null", "Sample Size",
                       "Number Controls", "poisthin", "OLS", "RUV2",
                       "RUVinv", "CATEnc", "vRUV4", "RUV3")

pi0hat_qv$`Sample Size` <- pi0hat_qv$`Sample Size` * 2


longdat <- melt(pi0hat_qv,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(pi0hat_qv)[-(1:5)])

pdf(file = "qv_all.pdf", height = 6, width = 10)
longdat_null <- longdat[longdat$`Proportion Null` == 1, ]
ggplot(longdat_null, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 1, color = "red", lty = 2) +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)


longdat_half <- longdat[longdat$`Proportion Null` == 0.9, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 0.9, color = "red", lty = 2) +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)

longdat_half <- longdat[longdat$`Proportion Null` == 0.5, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 0.5, color = "red", lty = 2) +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)

dev.off()


## ASH AUC ----------------------------------------------------------------------
rm(list = ls())
auc_ash <- read.csv("ash_auc_mat.csv")

names(auc_ash) <- c("Seed", "Proportion Null", "Sample Size",
                       "Number Controls", "poisthin", "OLS", "RUV2",
                       "RUVinv", "CATEnc", "vRUV4", "RUV3")

auc_ash$`Sample Size` <- auc_ash$`Sample Size` * 2

longdat <- melt(auc_ash,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(auc_ash)[-(1:5)])


pdf(file = "ash_all_auc.pdf", height = 6, width = 10)
longdat_half <- longdat[longdat$`Proportion Null` == 0.9, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("AUC") +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)

longdat_half <- longdat[longdat$`Proportion Null` == 0.5, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("AUC") +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)
dev.off()



## qvalue auc ---------------------------------------------------------
rm(list = ls())
auc_pvalues <- read.csv("qvalue_auc_mat.csv")
names(auc_pvalues) <- c("Seed", "Proportion Null", "Sample Size",
                       "Number Controls", "poisthin", "OLS", "RUV2",
                       "RUVinv", "CATEnc", "vRUV4", "RUV3")

auc_pvalues$`Sample Size` <- auc_pvalues$`Sample Size` * 2


longdat <- melt(auc_pvalues,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(auc_pvalues)[-(1:5)])

pdf(file = "pvalues_all_auc.pdf", height = 6, width = 10)
longdat_half <- longdat[longdat$`Proportion Null` == 0.9, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("AUC") +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)

longdat_half <- longdat[longdat$`Proportion Null` == 0.5, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("AUC") +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)

dev.off()



## mse ---------------------------------------------------------
rm(list = ls())
mse_mat <- read.csv("mse_mat.csv")
names(mse_mat) <- c("Seed", "Proportion Null", "Sample Size",
                       "Number Controls", "poisthin", "OLS", "RUV2",
                       "RUVinv", "CATEnc", "vRUV4", "RUV3", "RUVsoft")

mse_mat$`Sample Size` <- mse_mat$`Sample Size` * 2


longdat <- melt(mse_mat,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(mse_mat)[-(1:5)])

pdf(file = "all_mse.pdf", height = 6, width = 10)
longdat_half <- longdat[longdat$`Proportion Null` == 1, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("MSE") +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)


longdat_half <- longdat[longdat$`Proportion Null` == 0.9, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("MSE") +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)

longdat_half <- longdat[longdat$`Proportion Null` == 0.5, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid("Controls" + `Number Controls`~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("MSE") +
    geom_vline(xintercept = 9.5, color = "blue", lty = 2, alpha = 1/2)

dev.off()
