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
  effect_rank_onaski = onatski_gaussapprox(Y/glm_weights, family = glm_family, q_max, weights =glm_weights, max_iter = 20)$q_hat
  #effect_rank_onaski = onatski_rank(Y/glm_weights, family = glm_family, q_max, weights =glm_weights, max_iter = 20, fit_rank  = q_max)$q_hat
  
  return (list(q_act = effect_rank_act,
               q_onaski = effect_rank_onaski,
               ratio_ = ratio_))
  
}



