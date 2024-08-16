Trans_label_est_split <- function(fold, tau,  X_source_matrix, B_source_hat, beta_target_half,
                                  X_target, y_target, N_target, N_source,
                                  beta_hat_ini_project, beta_target_full, which_beta_to_use){
  L = ncol(B_source_hat)
  ### estimate b(l), Sigma_Q, Gamma_beta
  Sigma_Q_hat = t(X_target) %*% X_target / nrow(X_target)
  Gamma_beta_hat = t(B_source_hat) %*% Sigma_Q_hat %*% B_source_hat

  ### estimate initialization of gamma_l, s_Q, beta
  ini_all = est_ini(fold = fold, sample = 'split', nfold = NA, X_target, y_target,
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_half, L)

  s2_Q_hat_ini_target = ini_all$ini_target$s2_Q
  gamma_hat_ini_target = ini_all$ini_target$gamma

  s2_Q_hat_ini_linear = ini_all$ini_linear_comb$s2_Q
  gamma_hat_ini_linear = ini_all$ini_linear_comb$gamma

  s2_Q_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$s2_Q
  gamma_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$gamma

  s2_Q_hat_min = min(s2_Q_hat_ini_linear, s2_Q_hat_ini_target) + tau

  if(which_beta_to_use == 'target_only'){
    beta_hat_in_obj = beta_target_full
  }else{
    beta_hat_in_obj = beta_target_half
  }

  test_rid = ini_all$test_rid
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]
  beta_hat_in_constraint = beta_target_half

  ### estiamte beta
  coef_min_s2Q_robust_1 =
    Trans_label_robust_A_robust(X_target_test, y_target_test,
                                B_source_hat, Sigma_Q_hat, beta_hat_in_constraint,
                                beta_hat_in_obj,  s2_Q_hat_min, beta_hat_ini_project)

  coef_min_s2Q_robust_A_base = list(`gamma_hat` = coef_min_s2Q_robust_1$gamma_hat,
                                    `beta_hat` = coef_min_s2Q_robust_1$beta_hat,
                                    `success` = coef_min_s2Q_robust_1$success)
  result = list(
    `min_s2Q_robust_A_base_split` = coef_min_s2Q_robust_A_base,
    `ini_split` = ini_all
  )
  return(result)
}



Trans_label_est_split_ave <-  function(X_source_matrix, tau,  B_source_hat,
                                       beta_target_half_1, beta_target_half_2,
                                       X_target, y_target, N_target, N_source,
                                       beta_hat_ini_project, beta_target_full, which_beta_to_use){
  L = ncol(B_source_hat)
  ### estimate b(l), Sigma_Q, Gamma_beta
  Sigma_Q_hat = t(X_target) %*% X_target / nrow(X_target)
  Gamma_beta_hat = t(B_source_hat) %*% Sigma_Q_hat %*% B_source_hat

  ### estimate initialization of gamma_l, s_Q, beta
  ##!!!! notice that fold contains 0 and 1, not 1 and 2.
  ini_all = est_ini(fold = 1, sample = 'split', nfold = NA, X_target, y_target,
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_half_2, L)

  s2_Q_hat_ini_target = ini_all$ini_target$s2_Q
  gamma_hat_ini_target = ini_all$ini_target$gamma

  s2_Q_hat_ini_linear = ini_all$ini_linear_comb$s2_Q
  gamma_hat_ini_linear = ini_all$ini_linear_comb$gamma

  s2_Q_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$s2_Q
  gamma_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$gamma

  s2_Q_hat_min = min(s2_Q_hat_ini_linear, s2_Q_hat_ini_target) + tau
  # s2_Q_hat_min = s2_Q_hat_ini_linear_robust + tau

  if(which_beta_to_use == 'target_only'){
    beta_hat_in_obj = beta_target_full
  }else{
    beta_hat_in_obj = beta_target_half_2
  }

  test_rid = ini_all$test_rid
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]
  beta_hat_in_constraint = beta_target_half_2

  ### estiamte beta
  coef_min_s2Q_robust_1 =
    Trans_label_robust_A_robust(X_target_test, y_target_test,
                                B_source_hat, Sigma_Q_hat, beta_hat_in_constraint,
                                beta_hat_in_obj,  s2_Q_hat_min, beta_hat_ini_project)

  ### use the second half!
  ### estimate initialization of gamma_l, s_Q, beta
  ini_all = est_ini(fold = 0, sample = 'split', nfold = NA, X_target, y_target,
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_half_1, L)

  s2_Q_hat_ini_target = ini_all$ini_target$s2_Q
  gamma_hat_ini_target = ini_all$ini_target$gamma

  s2_Q_hat_ini_linear = ini_all$ini_linear_comb$s2_Q
  gamma_hat_ini_linear = ini_all$ini_linear_comb$gamma

  s2_Q_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$s2_Q
  gamma_hat_ini_linear_robust = ini_all$ini_robust_linear_comb$gamma

  s2_Q_hat_min = min(s2_Q_hat_ini_linear, s2_Q_hat_ini_target) + tau
  # s2_Q_hat_min = s2_Q_hat_ini_linear_robust + tau

  if(which_beta_to_use == 'target_only'){
    beta_hat_in_obj = beta_target_full
  }else{
    beta_hat_in_obj = ini_all$ini_target$beta
  }
  test_rid = ini_all$test_rid
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]

  ### estiamte beta
  coef_min_s2Q_robust_2 =
    Trans_label_robust_A_robust(X_target_test, y_target_test,
                                B_source_hat, Sigma_Q_hat, ini_all$ini_target$beta,
                                beta_hat_in_obj,  s2_Q_hat_min, beta_hat_ini_project)


  coef_min_s2Q_robust_A_robust = list(`gamma_hat` = (coef_min_s2Q_robust_1$gamma_hat + coef_min_s2Q_robust_2$gamma_hat)/2,
                                      `beta_hat` = (coef_min_s2Q_robust_1$beta_hat + coef_min_s2Q_robust_2$beta_hat)/2,
                                      `success` = (coef_min_s2Q_robust_1$success + coef_min_s2Q_robust_2$success)/2)
  result = list(
    `min_s2Q_robust_A_robust_split` = coef_min_s2Q_robust_A_robust,
    `ini_split` = ini_all
  )
  return(result)
}






