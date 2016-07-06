simulate:
    exec: pois_thin.R
    seed: R(1:2)
    params:
        Nsamp: 5
        nullpi: 0.9
        ncontrol: 100
    return: Y, X, num_sv, which_null, beta_true, control_genes

adjust:
    exec: ols.R, sva.R, ruv2.R, ruv4.R, ruvinv.R, vruv4.R, cate_nc.R, cate_rr.R
    params:
        Y: $Y
        X: $X
        num_sv: $num_sv
        control_genes: $control_genes
        exec[7:8]:
            calibrate: 0, 1
    return: betahat, sebetahat, df, pvalues

sum_stat:
    exec: fit_ash.R, fit_qvalue.R
    .alias: ash, qvalue
    params:
        betahat: $betahat
        sebetahat: $sebetahat
        df: $df
        model: "ET"
        pvalues: $pvalues
    return: betahat, lfdr, pi0hat

auc:
    exec: auc.R
    .alias: auc
    params:
        which_null: $which_null
        lfdr: $lfdr
    return: auc

mse:
    exec: mse.R
    .alias: mse
    params:
        beta_true: $beta_true
        betahat: $betahat
    return: mse


DSC:
    run: simulate *
         adjust *
         sum_stat *
         (auc, mse)
    exec_path: data_generators, adjustment_methods, summary_stat_methods, evaluations
    R_libs: limma, sva, qvalue, ruv, pROC, dcgerard/vicar, stephens999/ashr
    output: dsc_result