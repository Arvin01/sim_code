library(ggplot2)
library(reshape2)

## ASH results ----------------------------------------------------------------
pi0hat_ash <- read.csv("pi0hat_ash.csv")

names(pi0hat_ash) <- c("Seed", "Sample Size", "Proportion Null",
                       "Number Controls", "poisthin", "OLS", "vRUV4",
                       "SVA", "CATEnc", "CATErr", "SUCCOTASH")

pi0hat_ash <- pi0hat_ash[, c(1:6, 8, 7, 9:11)]
pi0hat_ash$`Sample Size` <- pi0hat_ash$`Sample Size` * 2

longdat <- melt(pi0hat_ash,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(pi0hat_ash)[-(1:5)])


pdf(file = "./fig/ash_all.pdf", height = 6, width = 10)
longdat_null <- longdat[longdat$`Proportion Null` == 1, ]
ggplot(longdat_null, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid(.~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 1, color = "red", lty = 2) +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)
dev.off()

pdf(file = "./fig/ash_all_nine.pdf", height = 6, width = 10)
longdat_half <- longdat[longdat$`Proportion Null` == 0.9, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid(.~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 0.9, color = "red", lty = 2) +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)
dev.off()


pdf(file = "./fig/ash_all_half.pdf", height = 6, width = 10)
longdat_half <- longdat[longdat$`Proportion Null` == 0.5, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid(.~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("Estimated Proportion Null") +
    geom_hline(yintercept = 0.5, color = "red", lty = 2) +
    geom_vline(xintercept = 7.5, color = "blue", lty = 2, alpha = 1/2)
dev.off()


## ASH AUC ----------------------------------------------------------------
rm(list = ls())
auc_ash <- read.csv("auc_ash.csv")

names(auc_ash) <- c("Seed", "Sample Size", "Proportion Null",
                       "Number Controls", "poisthin", "OLS", "vRUV4",
                       "SVA", "CATEnc", "CATErr", "SUCCOTASH")

auc_ash <- auc_ash[, c(1:6, 8, 7, 9:11)]
auc_ash$`Sample Size` <- auc_ash$`Sample Size` * 2

longdat <- melt(auc_ash,
                id.vars = c("Seed", "Sample Size", "Proportion Null", "Number Controls"),
                measure.vars = colnames(auc_ash)[-(1:5)])

dummydat <- aggregate(value ~ variable + `Sample Size` + `Proportion Null`,
                      data = longdat, FUN = median)


pdf(file = "./fig/ash_all_nine_auc.pdf", height = 6, width = 10)
dummydat_temp <- aggregate(data = dummydat[dummydat$`Proportion Null` == 0.9, ],
                           value ~ `Sample Size`, FUN = max)
longdat_half <- longdat[longdat$`Proportion Null` == 0.9, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid(.~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("AUC") +
    geom_hline(data = dummydat_temp, mapping = aes(yintercept = value), col = 2, lty = 2)
dev.off()


pdf(file = "./fig/ash_all_half_auc.pdf", height = 6, width = 10)
dummydat_temp <- aggregate(data = dummydat[dummydat$`Proportion Null` == 0.5, ],
                           value ~ `Sample Size`, FUN = max)
longdat_half <- longdat[longdat$`Proportion Null` == 0.5, ]
ggplot(longdat_half, mapping = aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    facet_grid(.~ "Sample Size" + `Sample Size`) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "none",
          text = element_text(size = 20)) +
    xlab("Methods") +
    ylab("AUC") +
    geom_hline(data = dummydat_temp, mapping = aes(yintercept = value), col = 2, lty = 2)
dev.off()
