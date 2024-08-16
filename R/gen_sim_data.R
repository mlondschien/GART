X_cov_source <- function(p){
  cov_ma = sapply(1:p, function(i){
    0.6^abs(i-(1:p))
  })
  diag(cov_ma) = diag(cov_ma) + stats::rnorm(p, mean = 0.1, sd = 0.01)
  # cov_ma = diag(x = rep(1, p))
  return(cov_ma)
}

X_cov_target <- function(p){
  cov_ma = sapply(1:p, function(i){
    0.7^abs(i-(1:p))
  })
  # cov_ma = diag(x = rep(1, p))
  return(cov_ma)
}



gen_coeff_bv <- function(alpha = 0.01){
  Beta_matrix =  matrix(c(0.3, 0.1, 0.5, -0.2, -0.7, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1, 0, 0,
                          0.2, 0.05, 0.4, -0.15, -0.6, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, -1, 0,0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, -1, 0,
                          0.3, 0.2, 0.6, -0.3, -0.6, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1,0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 1,
                          -0.2, -0.1, -0.5, 0.2, 0.6, 0, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0, 0, 0, 0, 0.2, 0, 0, 0, 0.3, 0, 0, 0, 0.4, 0, 0, 0), ncol = 4)

  target_beta = rowMeans(Beta_matrix[,-1]) + alpha * rep(c(-1, 1, -1, -1, 1), 7)
  return(list(Beta_matrix, target_beta))
}



#' Generate simulation data to run GART
#' @param alpha A numeric parameter controling the distance between the source convex hull and the target model.
#' @return A list of two elements including the input data for GART estimation and a validation dataset for GART evaluation. Note that for illustration purpose, the validation population shares the same linear model as the training target population.
#' @export
simu_data <- function(alpha = 0.07){
  coeff_matrix = gen_coeff_bv(alpha)
  X_lst = list()
  y_lst = list()
  ratio_target_n = 0.01
  N_vec = c(20000, 10000, 15000, 20000)
  L = length(N_vec)
  N_target = round(mean(N_vec)*ratio_target_n)
  s_Q=1
  s_P=.5
  for(l in 1:L){
    N = N_vec[l]
    beta = coeff_matrix[[1]][,l]
    p = length(beta)
    X_l = MASS::mvrnorm(n = N, mu = stats::rnorm(p, mean = stats::rexp(1), sd = 0.01),
                        Sigma = X_cov_source(p))
    epsilon = matrix(stats::rnorm(N, mean = 0, sd = s_P), N, 1)
    y_l = X_l %*% beta + epsilon
    X_lst = c(X_lst, list(X_l))
    y_lst = c(y_lst, list(y_l))
  }

  Sigma_Q = X_cov_target(p)
  X_target = MASS::mvrnorm(n = N_target, mu = rep(0, p),
                           Sigma = Sigma_Q)

  epsilon = matrix(stats::rnorm(N_target, mean = 0, sd = s_Q), N_target, 1)
  beta_target = c(`intercept` = 1, coeff_matrix[[2]])
  y_target = cbind(1, X_target) %*% beta_target + epsilon

  N_valid = 5000
  X_target_valid = MASS::mvrnorm(n = N_valid, mu = rep(0, p),
                                 Sigma = Sigma_Q)
  epsilon_valid = matrix(stats::rnorm(N_valid, mean = 0, sd = s_Q), N_valid, 1)
  beta_valid = c(1, coeff_matrix[[2]])
  y_target_valid = cbind(1, X_target_valid) %*% beta_valid + epsilon_valid

  target_data = list(`X` = X_target, `y` = y_target, `beta` = beta_target)
  source_data = list(`X` = X_lst, `y` = y_lst, `beta` = coeff_matrix[[1]])
  valid_data = list(`X` = X_target_valid, `y` = y_target_valid, `beta` = beta_valid)
  return(list(`GART_est_input` = list(`target` = target_data, `source` = source_data),
              `GART_eval_input` = valid_data))
}



