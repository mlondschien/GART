
est_ini_w_target <- function(X_target_train, y_target_train,
                             X_target_test, y_target_test, beta_target_half, L){
  n_train = length(y_target_train)
  n_test = length(y_target_test)
  s2_Q_hat = mean((y_target_test - X_target_test %*% beta_target_half)^2)
  gamma_ini = rep(1/L, L)
  return(list(`beta` = beta_target_half, `s2_Q` = s2_Q_hat, `gamma` = gamma_ini))
}

est_ini_w_linear_comb <- function(B_source_hat, beta_target_half,
                                  X_target_train, y_target_train,
                                  X_target_test, y_target_test, L){
  n_train = length(y_target_train)
  n_test = length(y_target_test)
  ### Use gamma to estimate beta_hat_ini!
  gammaHat <- CVXR::Variable(ncol(B_source_hat))
  objective <- CVXR::Minimize(mean((X_target_train %*% B_source_hat %*% gammaHat -
                                y_target_train)^2))
  problem <- CVXR::Problem(objective, constraints =
                       list(gammaHat >= 1e-10, gammaHat <= 1, sum(gammaHat)==1))
  result <- CVXR::solve(problem)
  if(result$status %in% c('optimal', 'optimal_inaccurate')){
    gamma_ini = as.vector(result$getValue(gammaHat))
    beta_hat_ini = B_source_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_test - X_target_test %*% beta_hat_ini)^2)
  }else{
    gamma_ini = rep(1/L, L)
    beta_hat_ini = B_source_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_test - X_target_test %*% beta_hat_ini)^2)
  }
  return(list(`beta` = beta_hat_ini, `s2_Q` = s2_Q_hat, `gamma` = gamma_ini))
}


est_ini_w_robust_linear_comb <- function(B_source_hat, beta_target_half,
                                         X_target_train, y_target_train,
                                         X_target_test, y_target_test, L){
  n_train = length(y_target_train)
  n_test = length(y_target_test)

  gammaHat <- CVXR::Variable(ncol(B_source_hat)+1)
  B_hat = cbind(beta_target_half, B_source_hat)

  objective <- CVXR::Minimize(mean((X_target_test %*% B_hat %*% gammaHat -
                                y_target_test)^2))
  problem <- CVXR::Problem(objective, constraints =
                       list(gammaHat >= 1e-10, gammaHat <= 1, sum(gammaHat)==1))
  result <- CVXR::solve(problem)

  if(result$status == 'optimal'){
    gamma_ini = as.vector(result$getValue(gammaHat))
    beta_hat_ini = B_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_train - X_target_train %*% beta_hat_ini)^2)
  }else{
    gamma_ini = rep(1/L, L)
    beta_hat_ini = B_hat %*% gamma_ini
    s2_Q_hat = mean((y_target_train - X_target_train %*% beta_hat_ini)^2)
  }
  return(list(`beta` = beta_hat_ini, `s2_Q` = s2_Q_hat, `gamma` = gamma_ini))
}



est_ini <- function(fold=1, sample = c('split', 'cv'), nfold, X_target, y_target,
                    Sigma_Q_hat, Gamma_beta_hat, B_source_hat, beta_target_half, L){
  N_target = nrow(X_target)

  test_rid = (1:N_target)[((1:N_target) %% 2) == fold]
  train_rid = setdiff(1:N_target, test_rid)
  X_target_train = X_target[train_rid, ]
  y_target_train = y_target[train_rid]
  X_target_test = X_target[test_rid, ]
  y_target_test = y_target[test_rid]
  ini_target = est_ini_w_target(X_target_train, y_target_train,
                                X_target_test, y_target_test, beta_target_half, L)
  ini_linear_comb = est_ini_w_linear_comb(B_source_hat, beta_target_half,
                                          X_target_train, y_target_train,
                                          X_target_test, y_target_test, L)
  ini_robust_linear_comb = est_ini_w_robust_linear_comb(B_source_hat, beta_target_half,
                                                        X_target_train, y_target_train,
                                                        X_target_test, y_target_test, L)
  return(list(`test_rid` = test_rid,
              `ini_target` = ini_target,
              `ini_linear_comb` = ini_linear_comb,
              `ini_robust_linear_comb` = ini_robust_linear_comb
  )
  )
}
