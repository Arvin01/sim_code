## Plots from mouthwash_sims.R
library(stringr)
library(ggplot2)
library(dplyr)

## auc matrix
auc_mat <- read.csv("auc_mat2.csv")
names(auc_mat) <- str_replace(names(auc_mat), pattern = "auc_", replacement = "")
auc_mat <- filter(select(auc_mat, -current_seed, -poisthin, -ncontrols), nullpi != 1)

mean_na <- function(x) mean(x, na.rm = TRUE)
auc_summ <- auc_mat %>% group_by(nullpi, Nsamp) %>% summarize_each(funs = funs(mean_na))

as.matrix(auc_summ[1,])
as.matrix(auc_summ[2,])

## pi0
pi0_mat <- read.csv("pi0_mat2.csv")
names(pi0_mat) <- str_replace(names(pi0_mat), "pi0_", "")
pi0_mat <- select(pi0_mat, -current_seed, -poisthin, -ncontrols)
pi0_sum <- pi0_mat %>% group_by(nullpi, Nsamp) %>% summarize_each(funs(mean_na))

as.matrix(pi0_sum[1,])
as.matrix(pi0_sum[2,])
as.matrix(pi0_sum[3,])
as.matrix(pi0_sum[4,])

## mse
mse_mat <- read.csv("mse_mat2.csv")
names(mse_mat) <- str_replace(names(mse_mat), "mse_", "")
mse_mat <- select(mse_mat, -current_seed, -poisthin, -ncontrols)
mse_sum <- mse_mat %>% group_by(nullpi, Nsamp) %>% summarize_each(funs(mean_na))

as.matrix(mse_sum[1,])
as.matrix(mse_sum[2,])
as.matrix(mse_sum[3,])
as.matrix(mse_sum[4,])
