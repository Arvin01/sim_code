library(ggplot2)
library(reshape2)

## ASH results ----------------------------------------------------------------
pi0hat_ash <- read.csv("pi0hat_ash.csv")
pi0hat_ruvinflate <- read.csv("pi0hat_ash_inflate.csv")
pi0hat_ash <- cbind(pi0hat_ash, pi0hat_ruvinflate[, -(1:5)])

names(pi0hat_ash) <- c("Seed", "Sample Size", "Proportion Null",
                       "Number Controls", "poisthin", "OLS", "RUV2",
                       "RUV4", "RUVinv", "VICAR", "SVA", "CATEnc",
                       "CATErr", "RUV4Inf", "RUV4InfEbayes", "RUV4Weird")

pi0hat_ash <- pi0hat_ash[, c(1:6, 11, 7:9, 12:13, 10, 14:16)]
pi0hat_ash$`Sample Size` <- pi0hat_ash$`Sample Size` * 2

longdat <- melt(pi0hat_ash,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(pi0hat_ash)[-(1:5)])


pdf(file = "ash_all_pi0hat.pdf", height = 6, width = 10)
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
dev.off()


## qvalue results ---------------------------------------------------------
rm(list = ls())
pi0hat_qv <- read.csv("pi0hat_qvalue.csv")
pi0hat_ruvinflate <- read.csv("pi0hat_qvalue_inflate.csv")
pi0hat_qv <- cbind(pi0hat_qv, pi0hat_ruvinflate[, -(1:5)])
names(pi0hat_qv) <- c("Seed", "Sample Size", "Proportion Null",
                       "Number Controls", "poisthin", "OLS", "RUV2",
                       "RUV4", "RUVinv", "VICAR", "SVA", "CATEnc", "Cal CATEnc",
                       "CATErr", "Cal CATErr", "RUV4Inf", "RUV4InfEbayes", "RUV4Weird")

pi0hat_qv <- pi0hat_qv[, c(1:6, 11, 7:9, 12:15, 10, 16:18)]
pi0hat_qv$`Sample Size` <- pi0hat_qv$`Sample Size` * 2


longdat <- melt(pi0hat_qv,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(pi0hat_qv)[-(1:5)])

pdf(file = "qv_pi0hat.pdf", height = 6, width = 10)
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
dev.off()
