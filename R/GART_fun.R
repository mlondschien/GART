#' Run GART main function
#'
#' @param data A list including two elements named `target` and `source`. Both contain two sublists named `X` and `y`, corresponding to the training data of source sites and the target site. Note that `X` and `y` in `source` should include the same number of sublists if multiple source sites exist.
#' @param base_names A vector of names representing different initial estimators for GART. Possible choices include 'target' (target only estimator), 'zero', 'convex' (convex combination of sources and target), 'weighted' (weighted average of sources and target), 'source' (convex combination of sources).
#' @param tau A number constroling the size of constraint set in GART. Default is 1/n_target.
#' @param is_benchmark If \code{TRUE}, calculate benchmark methods including target only estimator, source mixture, maximin, transLasso and transGLM. Default is \code{FALSE}.
#' @param is_valid If \code{TRUE}, evaluate model performance based on `data_valid`. Default is \code{FALSE}.
#' @param data_valid A list including three elements for the validation data if the true coefficient generating validation data is known, namely `X` (an N by p matrix), `y` (a vector of length N) and `beta` (a vector of length p or p+1 (with intercept)). If the true coefficient is unknown, only the first two elements are enough.
#' @returns A list containing estimation and/or evaluation results by running GART (and other benchmark methods if `is_benchmark=TRUE`).
#' @export
GART <- function(data, base_names = 'weighted', tau = NULL,
                 is_benchmark = F, is_valid = F, data_valid = NULL){
  est_model = GART_est(data, base_names, tau = tau,
                       is_benchmark = is_benchmark)
  if(is_benchmark){
    beta_est = cbind(est_model$coef_beta, est_model$other_model)
  }else{
    beta_est = as.matrix(est_model$coef_beta)
  }
  result = list(`est` = est_model)
  if(is_valid){
    message('Valid model performance...')
    eval = GART_eval(beta_est, data_valid)
    result = list(`est` = est_model,
                  `eval` = eval)
  }
  return(result)
}

GART_est <- function(data, base_names, tau = NULL, is_benchmark = F){
  X_target = data$target$X; y_target = data$target$y
  X_source = data$source$X; y_source = data$source$y
  N_source = sapply(X_source, nrow); N_target = nrow(X_target)
  if(is.null(tau)){
    tau = 1/length(y_target)
  }

  L = length(X_source)
  p = ncol(X_target)
  trans_intercept = F

  ### estimate b(l), Sigma_Q, Gamma_beta
  message("Estimate source coef...")
  X_source_matrix = Reduce(rbind, X_source)
  y_source_vec = unlist(y_source)
  B_source_hat = sapply(1:L, function(l){
    est = glmnet::cv.glmnet(X_source[[l]], y_source[[l]], standardize=F)
    est = glmnet::glmnet(X_source[[l]], y_source[[l]],
                         lambda = est$lambda.min, standardize=F)
    b_l = stats::coef(est)
    return(as.vector(b_l))
  })
  if(!trans_intercept){
    B_source_intercept = B_source_hat[1, ]
    B_source_hat = B_source_hat[-1, ]
  }

  est = glmnet::cv.glmnet(X_target, y_target, standardize=F)
  est = glmnet::glmnet(X_target, y_target, lambda = est$lambda.min, standardize=F)
  b_l = stats::coef(est)
  if(!trans_intercept){
    beta_target_full_intercept = b_l[1]
    b_l = b_l[-1]
  }
  beta_target_full = b_l
  beta_est = gamma_est = c()

  message('Estimate GART coef...')
  # for(which_beta_to_use in c('target_only', 'ini_beta')){
  for(which_beta_to_use in c('target_only')){
    base_est = GART_base_est = beta_target_half_lst = c()
    for(fold in c(0, 1)){
      test_rid = (1:N_target)[((1:N_target) %% 2) == fold]
      train_rid = setdiff(1:N_target, test_rid)
      est = glmnet::cv.glmnet(X_target[train_rid, ], y_target[train_rid], standardize=F)
      est = glmnet::glmnet(X_target[train_rid, ], y_target[train_rid], lambda = est$lambda.min, standardize=F)
      b_l = stats::coef(est)
      if(!trans_intercept){
        b_l = b_l[-1]
      }
      beta_target_half = b_l
      beta_target_half_lst = c(beta_target_half_lst, list(beta_target_half))

      gammaHat <- CVXR::Variable(ncol(B_source_hat)+1)
      B_hat_1 = cbind(beta_target_half, B_source_hat)
      objective <- CVXR::Minimize(mean((X_target[test_rid, ] %*% B_hat_1 %*% gammaHat -
                                    y_target[test_rid])^2))
      GART_step1 = step1 = c()
      ini_bound = 1
      for(base in base_names){
        problem = switch(base,
                         convex = CVXR::Problem(objective, constraints =
                                            list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1)),
                         source = CVXR::Problem(objective, constraints =
                                            list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1, gammaHat[1] == 0)),
                         mixST = CVXR::Problem(objective, constraints =
                                           list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1, gammaHat[1] == 0)),
                         mixT0 = CVXR::Problem(objective, constraints =
                                           list(gammaHat[1]>=0, gammaHat[2:(L+1)] >= -ini_bound, gammaHat <= ini_bound)),
                         target = CVXR::Problem(objective, constraints =
                                            list(gammaHat >= 0, gammaHat <= 1, sum(gammaHat)==1, gammaHat[1] == 1)),
                         weighted = CVXR::Problem(objective, constraints =
                                              list(gammaHat[1]>=0, gammaHat[2:(L+1)] >= -ini_bound, gammaHat <= ini_bound)),
                         zero = CVXR::Problem(objective, constraints =
                                          list(gammaHat == 0)))
        result1 <- CVXR::solve(problem)
        if(result1$status %in% c('optimal', 'optimal_inaccurate')){
          gamma_step1_1 = as.vector(result1$getValue(gammaHat))
        }else{
          # print('Failed!')
          gamma_step1_1 = c(1, rep(0, L))
        }
        beta_hat_ini_project = B_hat_1 %*% gamma_step1_1
        if(base == 'mixT0'){
          gamma_weight = gamma_step1_1

          sample_mixT = sample(1:nrow(X_target), size = round((1-tau)*N_target))
          sample_0 = sample(1:nrow(X_target), size = round((tau)*N_target))
          X_ini = rbind(X_target[sample_mixT, ], X_target[sample_0, ])
          Y_ini = c(X_target[sample_mixT, ] %*% B_hat_1 %*% gamma_weight,
                    X_target[sample_0, ] %*% rep(0, p))

          beta = CVXR::Variable(p)
          model = CVXR::Problem(CVXR::Minimize(mean((Y_ini - X_ini %*% beta)^2)))
          result <- CVXR::solve(model)
          beta_hat_ini_project = result$getValue(beta)
          mean((beta_hat_ini_project - data$beta_valid)^2)
        }
        if(base == 'mixST'){
          gamma_source = gamma_step1_1

          gammaHat <- CVXR::Variable(ncol(B_source_hat)+1)
          B_hat_1 = cbind(beta_target_half, B_source_hat)
          objective <- CVXR::Minimize(mean((X_target[test_rid, ] %*% B_hat_1 %*% gammaHat -
                                        y_target[test_rid])^2))
          problem = CVXR::Problem(objective, constraints =
                              list(gammaHat[1]>=0, gammaHat[2:(L+1)] >= -ini_bound, gammaHat <= ini_bound))
          result = CVXR::solve(problem)
          gamma_weight = result$getValue(gammaHat)

          sample_mixST = sample(1:nrow(X_target), size = round((1-tau)*N_target))
          sample_mixsource = sample(1:nrow(X_target), size = round((tau)*N_target))
          X_ini = rbind(X_target[sample_mixST, ], X_target[sample_mixsource, ])
          Y_ini = c(X_target[sample_mixST, ] %*% B_hat_1 %*% gamma_weight,
                    X_target[sample_mixsource, ] %*% B_source_hat %*% gamma_source[-1])


          beta = CVXR::Variable(p)
          model = CVXR::Problem(CVXR::Minimize(mean((Y_ini - X_ini %*% beta)^2)))
          result <- CVXR::solve(model)
          beta_hat_ini_project = result$getValue(beta)
          mean((beta_hat_ini_project - data$beta_valid)^2)
        }
        if(base == 'target'){
          beta_hat_ini_project = beta_target_full
        }
        result = list(`gamma` = gamma_step1_1, `beta` = beta_hat_ini_project)
        step1 = c(step1, list(result))
        result = Trans_label_est_split(fold = fold, tau, X_source_matrix, B_source_hat, beta_target_half,
                                       X_target, y_target, N_target, N_source,
                                       beta_hat_ini_project, beta_target_full, which_beta_to_use)
        GART_step1 = c(GART_step1, list(result))
        # if(base == base_names[1]){
        #   print(paste0('s2_Q in target only:', round(result$ini_split$ini_target$s2_Q, 3)))
        #   print(paste0('s2_Q in linear comb:', round(result$ini_split$ini_linear_comb$s2_Q, 3)))
        #   print(paste0('s2_Q in robust linear comb:', round(result$ini_split$ini_robust_linear_comb$s2_Q, 3)))
        # }
      }
      names(GART_step1) = names(step1) = base_names
      GART_base_est = c(GART_base_est, list(GART_step1))
      base_est = c(base_est, list(step1))
    }

    GART_base_ave_later_lst = base_lst = c()
    for(base in base_names){
      GART_base_avg = base_avg = list()
      GART_base_avg$gamma_hat = (GART_base_est[[1]][[base]]$min_s2Q_robust_A_base_split$gamma_hat +
                                      GART_base_est[[2]][[base]]$min_s2Q_robust_A_base_split$gamma_hat) / 2
      GART_base_avg$beta_hat = (GART_base_est[[1]][[base]]$min_s2Q_robust_A_base_split$beta_hat +
                                     GART_base_est[[2]][[base]]$min_s2Q_robust_A_base_split$beta_hat) / 2
      base_avg$gamma_hat = (base_est[[1]][[base]]$gamma + base_est[[2]][[base]]$gamma) / 2
      base_avg$beta_hat = (base_est[[1]][[base]]$beta + base_est[[2]][[base]]$beta) / 2
      GART_base_ave_later_lst = c(GART_base_ave_later_lst, list(GART_base_avg))
      base_lst = c(base_lst, list(base_avg))
    }
    names(base_lst) = paste0('base_', base_names)
    names(GART_base_ave_later_lst) = paste0('GART_ave_later_base_', base_names, '_', which_beta_to_use)

    GART_base_ave_first_lst = c()
    beta_target_half_1 = beta_target_half_lst[[1]]
    beta_target_half_2 = beta_target_half_lst[[2]]
    for(base in names(base_lst)){
      beta_hat_ini_project = base_lst[[base]]$beta
      result = Trans_label_est_split_ave(X_source_matrix, tau, B_source_hat, beta_target_half_1, beta_target_half_2,
                                         X_target, y_target, N_target, N_source,
                                         beta_hat_ini_project, beta_target_full, which_beta_to_use)
      junk = list(`beta_hat` = result$min_s2Q_robust_A_robust_split$beta_hat,
                  `gamma_hat` = result$min_s2Q_robust_A_robust_split$gamma_hat)
      GART_base_ave_first_lst = c(GART_base_ave_first_lst, list(junk))
    }
    # names(GART_base_ave_first_lst) = paste0('GART_ave_first_base_', base_names, '_', which_beta_to_use)
    names(GART_base_ave_first_lst) = paste0('GART_', base_names)

    beta_est = cbind(beta_est, cbind(
      sapply(GART_base_ave_first_lst, function(x) x$beta_hat)
      # sapply(GART_base_ave_later_lst, function(x) x$beta_hat)
    ))
    gamma_est = cbind(gamma_est, cbind(
      sapply(GART_base_ave_first_lst, function(x) x$gamma_hat)
      # sapply(GART_base_ave_later_lst, function(x) x$gamma_hat)
    ))
  }

  if(is_benchmark){
    message('Add benchmark methods (maximin, transLasso, transGLM, etc)...')
    add_benchmark = compare_est(B_source_hat, X_source, y_source, X_target, y_target)
    beta_est = rbind(`intercept` = beta_target_full_intercept, beta_est)
    other_model = add_benchmark
  }else{
    beta_est = rbind(`intercept` = beta_target_full_intercept, beta_est)
    other_model = data.frame(`target_only` = c(beta_target_full_intercept, beta_target_full))
  }

  return(list(`coef_beta` = beta_est, `weight_gamma` = gamma_est,
              `other_model` = other_model))
}
