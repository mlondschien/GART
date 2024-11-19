compare_est <- function(B_source_hat, X_source, y_source, X_target, y_target){
  L = length(X_source)

  print("target only")
  ## target only estimator
  out = glmnet::cv.glmnet(X_target, y_target, standardize=F)
  out = glmnet::glmnet(X_target, y_target, lambda = out$lambda.min, intercept = T, standardize = F)
  beta_target = stats::coef(out)

  print("linear source estimator")
  ## best linear source estimator
  Sigma_Q_hat = t(X_target) %*% X_target / nrow(X_target)
  Gamma_beta_hat = t(B_source_hat) %*% Sigma_Q_hat %*% B_source_hat
  n_target = length(y_target)

  gammaHat <- CVXR::Variable(ncol(B_source_hat))
  objective <- CVXR::Minimize(mean((X_target %*% B_source_hat %*% gammaHat -
                                y_target)^2))
  problem <- CVXR::Problem(objective, constraints =
                       list(gammaHat >= 1e-10, sum(gammaHat)==1))
  result <- CVXR::solve(problem)

  gamma_linear_source = as.vector(result$getValue(gammaHat))
  beta_linear_source = c(beta_target[1], B_source_hat %*% gamma_linear_source)

  print("maximin")
  ## maximin
  gammaHat <- CVXR::Variable(ncol(B_source_hat))
  objective <- CVXR::Minimize(CVXR::quad_form(gammaHat, Gamma_beta_hat))
  problem <- CVXR::Problem(objective, constraints =
                       list(gammaHat >= 1e-10, sum(gammaHat)==1))
  result <- CVXR::solve(problem)
  gamma_maximin = as.vector(result$getValue(gammaHat))
  beta_maximin = c(beta_target[1], B_source_hat %*% gamma_maximin)

  print("TransGLM")
  ## TransGLM
  # source_data = sapply(1:L, function(i){
  #   list(`x` = X_source[[i]], `y` = y_source[[i]])
  # }, simplify = FALSE)
  # target_data = list(`x` = X_target, `y` = y_target)
  # out_transglm = glmtrans::glmtrans(target = target_data,  source = source_data, family = "gaussian",
  #                         intercept = TRUE, detection.info = FALSE)
  # beta_transglm = out_transglm$beta
  beta_transglm = rep(0, length(beta_target))

  print("TransLasso")
  ## TransLasso
  n.vec = c(nrow(X_target), sapply(X_source, FUN = nrow))
  X_source_matrix = do.call(rbind, X_source)
  y_source_vec = do.call(rbind, y_source)
  prop.re1 <- tryCatch(Trans.lasso(rbind(X_target, X_source_matrix),
                                   c(y_target, y_source_vec), n.vec,
                                   I.til = 1:round(n_target/2), l1 = T), error = function(e) NULL)

  prop.re2 <- tryCatch(Trans.lasso(rbind(X_target, X_source_matrix),
                                   c(y_target, y_source_vec), n.vec,
                                   I.til = (round(n_target/2)+1):n_target, l1=T), error = function(e) NULL)
  if(!is.null(prop.re1) & !is.null(prop.re2)){
    beta.prop <- c(beta_target[1], (prop.re1$beta.hat + prop.re2$beta.hat) / 2)
    gamma.prop = (prop.re1$theta.hat + prop.re2$theta.hat) / 2
  }else{
    beta.prop = rep(NA, length(beta_target))
    gamma.prop = NA
  }
  result = data.frame(
    `target_only` = as.vector(beta_target),
    `linear_source_beta` = as.vector(beta_linear_source),
    `maximin` = as.vector(beta_maximin),
    `trans_lasso` = as.vector(beta.prop),
    `trans_glm` = as.vector(beta_transglm))
  return(result)
}


GART_eval <- function(beta_est, data_valid){
  if(nrow(beta_est) == ncol(data_valid$X) + 1){
    X_target_valid = cbind(1, data_valid$X)
  }else{
    X_target_valid = data_valid$X
  }
  y_target_valid = data_valid$y
  beta_valid = data_valid$beta

  predict_mse0 = as.vector(y_target_valid) - X_target_valid %*% as.matrix(beta_est)
  predict_mse = apply(predict_mse0, 2, function(x){
    mean(x^2)
  })
  R2 = apply(predict_mse0, 2, function(x){
    1 - mean(x^2) / stats::var(y_target_valid)
  })

  if(!is.null(beta_valid)){
    if(length(beta_valid) == nrow(beta_est)-1){
      beta_valid = c(0, beta_valid)
    }
    epsilon = X_target_valid %*% (as.matrix(beta_est) - as.vector(beta_valid))
    mse = apply(epsilon, 2, function(x){
      mean(x^2)
    })
    beta_diff2_norm = sqrt(apply(beta_est - as.vector(beta_valid), 2, function(x){
      sqrt(sum(x^2))
    }))
    result = list(`beta_diff2_norm` = beta_diff2_norm,
                  `mse` = mse,
                  `predict_mse` = predict_mse,
                  `R2` = R2)
  }else{
    result = list(`predict_mse` = predict_mse,
                  `R2` = R2)
  }
  return(result)
}

