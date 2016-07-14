get_mse <- function(betahat, beta_true) {
    mse <- mean((betahat - beta_true) ^ 2)
    return(mse)
}
