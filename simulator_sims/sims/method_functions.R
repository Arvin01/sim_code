## @knitr methods

R.utils::sourceDirectory("./adjustment_methods")

fit_ols <- new_method("ols", "Ordinary Least Squares",
                       method = function(model, draw) {
                         list(fit = ols(draw))
                       })
