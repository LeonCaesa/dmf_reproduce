setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)


case = 'Case3'
n <- 500
p <- 50
q_star <- 6
#q_star <- 15
#q_max = p-5
phi = 2
glm_family = negative.binomial(phi)
#glm_family = Gamma('log')
#glm_family = gaussian()
#glm_family = binomial()
d = p

if (glm_family$family=='binomial'){
  #glm_weights = matrix(rpois(n * d, lambda = n), nrow=n, ncol = d)
  glm_weights = 10
  meanshift = 0
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
  #Z_star <- matrix(runif(q_star * n, -n/10, n/10), nrow = n, ncol = q_star)  
  Z_star <- matrix(runif(q_star * n, -n/50, n/50), nrow = n, ncol = q_star)  
}else if (case == 'Case5' | case == 'Case6'){
  W_star <- diag(rep(1,d))[1:q_star,]
  Z_star <- matrix(runif(q_star * n, -0.1, 0.1), nrow = n, ncol = q_star)
}else{
  stop('Define the right cases')
}

# 
# onatski_improved = function(x, family, q_max, fit_rank,  weights = 1, offset = zeros(x), 
#                             max_iter =10){
#   n <- nrow(x); p <- ncol(x)
#   if (p < q_max + 5 ) stop('decrease rmax')
#   
#   fit_result = dmf(x, family, rank= fit_rank, weights, offset) 
#   
#   eta = tcrossprod(fit_result$L, fit_result$V)
#   cov_matrix = cov(eta)
#   lambda_hats = eigen(cov_matrix)$value
#   
#   tol = 1e3
#   j_update = q_max + 1
#   while (tol>=1){
#     j = j_update
#     
#     y_reg = lambda_hats[j:(j+4)]
#     x_reg = (j + seq(-1, 3, 1))^(2/3)
#     delta_temp = as.numeric(2* abs(coef(lm(y_reg ~ x_reg))))[2]
#     flag = which((-diff(lambda_hats)>= delta_temp))
#     if(length(flag) == 0){q_hat = 0}else{q_hat = tail(flag[flag<=q_max],1)}
#     j_update = q_hat + 1
#     tol = abs(j_update -j)
#   }
#   return (list(q_hat = q_hat, delta = delta_temp, L =fit_result$L,
#                V = fit_result$V, deviance = fit_result$deviance, family = fit_result$family))
# }

eta_star <- crossprod(t(Z_star), W_star) + meanshift

Y <- generate_Y(glm_family, eta_star, glm_weights = glm_weights, phi = phi)

eta_diff = rep(0, p-1)

# dmf_star =  dmf(Y, glm_family, q_star)
# eta_qstar = tcrossprod(dmf_star$L, dmf_star$V)

# for (q_fit in (p-5):2){
for (q_fit in p:2){

    dmf_fit = dmf(Y/glm_weights, glm_family, q_fit, weights = glm_weights)
    eta_hat = tcrossprod(dmf_fit$L,dmf_fit$V)  
    cov_hat = cov(eta_hat)
    lambda_hats = eigen(cov_hat)$value
    plot(lambda_hats, main = q_fit)
    abline(v = q_star)
    abline(v = q_fit)
    
    #eta_diff[q_fit-1] = mean( (eta_hat - eta_star)^2)
    
    
      
    #eta_diff[q_fit-1] = mean( (eta_hat - eta_qstar)^2)

    #rank_result = onatski_improved(Y/glm_weights, glm_family, q_max = p-5, glm_weights, max_iter =20, fit_full = TRUE)
    #rank_result = onatski_improved(Y/glm_weights, glm_family, q_max = 10, fit_rank = q_fit, glm_weights, max_iter =20)
    #effect_rank = rank_result$q_hat
    
    #print(c(q_fit, effect_rank, effect_rank == q_fit))
    
    
    #print( c(q_fit, eta_diff[q_fit-1], dmf_fit$deviance))
    #break
}

plot(2:p, eta_diff)
abline(v = q_star)
# eta_hat - glm_family$linkfun(Y)

