rank_simu<- function(d, q_star, n, phi, glm_family, q_max = 10, case = 'Case1'){
  if (glm_family$family=='binomial'){
    glm_weights = 10; meanshift = 0
  }else{
    glm_weights = matrix(rep(1, n * d), nrow = n, ncol = d)
    if(glm_family$family == 'poisson'){meanshift = 1
    if (glm_family$link=='sqrt') {meanshift = 5}
    }else{
      meanshift = 1 }
  }
  if (case == 'Case1' | case == 'Case2'){
    W_star <- matrix(rnorm(q_star * d, 0, sqrt(1)), nrow = q_star, ncol = d)
    Z_star <- matrix(rnorm(n * q_star, 0, sqrt(1)), nrow = n, ncol = q_star)
  }else if (case == 'Case3' | case == 'Case4'){
    W_star <- diag(rep(1,d))[1:q_star,]
    Z_star <- matrix(runif(q_star * n, -n/50, n/50), nrow = n, ncol = q_star)
  }else if (case == 'Case5' | case == 'Case6'){
    W_star <- diag(rep(1,d))[1:q_star,]
    Z_star <- matrix(runif(q_star * n, -0.1, 0.1), nrow = n, ncol = q_star)
  }else{
    stop('Define the right cases')
  }
  
  eta_star <- crossprod(t(Z_star), W_star) + meanshift
  
  param_star = dmf_identify(eta_star, q_star)
  ratio_ = c( round(log(diag(crossprod(param_star$L))[q_star])/log(n),4), glm_family$family)
  
  
  mu_star <- glm_family$linkinv(eta_star)
  
  Y <- generate_Y(glm_family, eta_star, glm_weights = glm_weights, phi = phi)
  
  
  effect_rank_act = act_rank(Y, glm_family, glm_weights)
  effect_rank_onaski = onatski_correction(Y/glm_weights, family = glm_family, q_max, weights =glm_weights, max_iter = 20)$q_hat
  
  return (list(q_act = effect_rank_act,
               q_onaski = effect_rank_onaski,
               ratio_ = ratio_))
  
}



onatski_correction <- function(x, family = gaussian(), q_max, weights =1,
                               offset= zeros(x), max_iter =10){
  n <- nrow(x); p <- ncol(x); d = ncol(x)
  if (n < p) stop("fewer observations than predictors")
  if (length(weights) == 1)
    weights <- rep(weights, n * p)
  else {
    if (is.vector(weights)) {
      if (length(weights) != n)
        stop("inconsistent number of weights")
      else
        weights <- rep(weights, p)
    } else {
      if (nrow(weights) != n && ncol(weights) != p)
        stop("inconsistent number of weights")
    }
  }
  valid <- (weights > 0) & (!is.na(x))
  if (any(apply(!valid, 1, all)) || any(apply(!valid, 2, all)))
    stop("full row or column with zero weights or missing values")
  
  invlink_Y = family$linkfun(x/weights)
  invlink_Y[invlink_Y == Inf] = max(invlink_Y[which(invlink_Y < Inf)])
  invlink_Y[invlink_Y == -Inf] = min(invlink_Y[which(invlink_Y > -Inf)])
  if (sum(is.na(invlink_Y))!=0) {
    warning("NA produced in g(Y) transformation, imputing with mean")
    invlink_Y[is.na(invlink_Y)] = mean(invlink_Y, na.rm =TRUE)
  }
  fit_result = svd(invlink_Y, nu = p, nv = p)
  fit_result$L =sweep(fit_result$u, 2, fit_result$d, '*')
  fit_result$V = fit_result$v
  
  eta = tcrossprod(fit_result$L, fit_result$V)
  cov_matrix = cov(eta)
  lambda_hats = eigen(cov_matrix)$value
  tol = Inf
  j_update = q_max + 1
  while (tol>=1){
    j = j_update
    y_reg = lambda_hats[j:(j+4)]
    x_reg = (j + seq(-1, 3, 1))^(2/3)
    delta_temp = as.numeric(2* abs(coef(lm(y_reg ~ x_reg))))[2]
    flag = which(-diff(lambda_hats)>= delta_temp)
    if(length(flag) == 0){q_hat = 0}else{q_hat = tail(flag[flag<=q_max],1)}
    j_update = q_hat + 1
    tol = abs(j_update -j)
  }
  
  return (list(q_hat = q_hat, delta = delta_temp, L =fit_result$L,
               V = fit_result$V, deviance = fit_result$deviance, family = fit_result$family))
  
}