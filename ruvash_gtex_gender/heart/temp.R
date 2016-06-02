
## Comparison to SUCCOTASH
Now run SUCCOTASH and see if we get anything very different.
```{r, cache = TRUE}
succ_out <- succotashr::succotash(Y = t(heart$Y), X = heart$X, k = k)
```

```{r}
succ_onsex <- onsex[order(succ_out$lfdr)]
```

SUCCOTASH says that `r sum(succ_onsex[1:num_sig])` of the top `r num_sig` are on a sex chromosome. The ordering between SUCCOTASH and RUVASH is also very similar for the top few genes. The rankings between SUCCOTASH and RUVASH are much more similar for all of the genes than between RUVASH and RUV4.

```{r}
qplot(rank(c(succ_out$lfdr))[lfdr_order[1:num_sig]], rank(ruvash_out$lfdr)[lfdr_order[1:num_sig]],
     xlab = "SUCCOTASH Rank", ylab = "RUVASH Rank", main = "Rankings",
     color = onsex[lfdr_order[1:num_sig]]) +
   guides(color=guide_legend(title="On Sex?")) +
    ggtitle(paste("Rankings for top", num_sig, "genes"))

qplot(rank(c(succ_out$lfdr)), rank(ruvash_out$lfdr),
     xlab = "SUCCOTASH Rank", ylab = "RUVASH Rank", main = "Rankings",
     color = onsex) +
   guides(color=guide_legend(title="On Sex?")) +
    ggtitle("Rankings for all genes")
```

SUCCOTASH says that `r sum(succ_out$lfdr < 0.5)` of the genes have a posterior probability less than a half of being zero.


The estimate of $\pi_0$ for SUCCOTASH is `r succ_out$pi0` and that of RUVASH is `r ruvash_out$fitted.g$pi[1]`. When I run qvalue with the p-values from RUV4 I get `r qout$pi0`. Which is very different from the estimates of RUVASH and SUCCOTASH.


The estimates of the variance scaling parameters are very similar. SUCCOTASH inflates the variance by `r succ_out$scale_val`, while RUVASH inflates the variance by `r ruvash_out$ruv$multiplier`.
