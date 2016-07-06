if (length(unique(which_null)) == 1) {
    auc <- NA
} else {
    auc <- pROC::roc(predictor = lfdr, response = which_null)$auc
}
score <- auc
