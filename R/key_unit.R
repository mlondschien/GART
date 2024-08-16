Trans_label_robust_A_robust <- function(X_target_test, y_target_test,
                                        B_source_hat, Sigma_Q_hat, beta_hat_in_constraint,
                                        beta_hat_in_obj, s2_Q_hat_min, beta_hat_ini_project){
  B_hat_old = as.matrix(cbind(beta_hat_in_constraint, B_source_hat))
  B_hat_new = as.matrix(cbind(beta_hat_in_obj, B_source_hat))
  B_hat_new_fun = as.matrix(cbind(beta_hat_in_obj, B_source_hat)) - as.vector(beta_hat_ini_project)

  Gamma_beta_hat_new = t(B_hat_new_fun) %*% Sigma_Q_hat %*% B_hat_new_fun

  gammaHat <- CVXR::Variable(ncol(B_hat_new))
  objective <- CVXR::Minimize(CVXR::quad_form(gammaHat, Gamma_beta_hat_new))
  constraints0 = mean((X_target_test %*% B_hat_old %*% gammaHat - y_target_test)^2) - s2_Q_hat_min <= 0

  problem = CVXR::Problem(objective, constraints =
                      list(constraints0, gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1))
  result <- CVXR::solve(problem)
  if(result$status %in% c('optimal', 'optimal_inaccurate')){
    gamma_hat = as.vector(result$getValue(gammaHat))
    beta_hat = B_hat_new %*% gamma_hat
    success = 1
  }else{
    gamma_hat = beta_hat = NA
    success = 0
    print('Failed!')
  }
  return(list(`gamma_hat` = gamma_hat,
              `beta_hat` = beta_hat,
              `success` = success))
}




