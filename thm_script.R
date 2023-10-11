setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)
library(rstiefel)

#p_lists =  seq(100, 400, by = 100)
#p_lists =  seq(10, 40, by = 10)
#n_lists =  seq(1000, 5000, by = 1000)
#n_lists =  seq(500, 1500, by = 200)
n_lists =  seq(200, 1000, by = 200)
#n_lists =  seq(10000, 50000, by = 10000)
#p_lists =  seq(10, 100, by = 20)
#p_lists =  seq(10, 70, by = 20)
#p_lists =  seq(100, 500, by = 100)
p_lists =  seq(100, 200, by = 20)
#p_lists =  seq(500, 1000, by = 100)


#phi = 50
phi = 5


family_list = c()
# family_list[[1]] = negative.binomial(phi)
# family_list[[2]] = poisson()
# family_list[[3]] = gaussian(0.3)
family_list[[1]] = binomial() 

q_star <- 6
n_runs = 10
#power_ = 0.3
power_ = 1
glm_weights = 10
control = glm.control(epsilon = 1e-8, maxit = 300)

set.seed(1)
# 
# n = 500
# p = 50
# glm_family = family_list[[1]]
result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results'
for (glm_family in family_list){
  family_dir = paste(result_dir, glm_family$family, sep = '/')
  dir.create(file.path(family_dir), showWarnings = FALSE)
    for (n in n_lists){
        for(p in p_lists){
          # create model dir
            np_dir = paste(family_dir, paste(n, p, sep = '_'), sep = '/')
            dir.create(file.path(np_dir), showWarnings = FALSE)
            
            for (run in 1:n_runs){
                #dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
                # create model dir
                #param_star = gene_trueLV(q_star, n, p, mu_v =0, sd_v =1 , mu_l= 0, sd_l = 1)
                #param_star = gene_trueLV(q_star, n, p, mu_v =0, sd_v =0.3 , mu_l= 0, sd_l = 0.3)
                #param_star = gene_trueLV(q_star, n, p, mu_v =1, sd_v =0.3 , mu_l= 1, sd_l = 0.3)
                #param_star = gene_trueLV(q_star, n, p, mu_v =1, sd_v =1 , mu_l= 1, sd_l = 1)
                #param_star = gene_trueLV(q_star, n, p, mu_v =1, sd_v =0.3 , mu_l= 1, sd_l = 0.3)
                #param_star = gene_trueLV(q_star, n, p, mu_v =1, sd_v =0.3 , mu_l= 1, sd_l = 0.3)
                
                param_star = gene_trueLV(q_star, n, p, power_)
                # crossprod(param_star$L_star)
                # crossprod(param_star$V_star)
              
                eta_star = tcrossprod(param_star$L_star,param_star$V_star)
                
                Y <- generate_Y(glm_family, eta_star, glm_weights = glm_weights, phi = phi)
                param_star$Y = Y/glm_weights
                
                # [DMF fit]
                dmf_result = dmf(param_star$Y, family = glm_family, rank = q_star, weights =  glm_weights)
                # [Save fitted result]
                param_savedir = paste(np_dir, paste('run',run, '_param_star.RData', sep = ''), sep = '/')
                dmf_savedir = paste(np_dir, paste('run',run, '_dmf.RData', sep = ''), sep = '/')
                save(param_star, file = param_savedir)
                save(dmf_result, file = dmf_savedir)
            }# end of n runs
        }# end of p loop
    }# end of n loop
}# end of family loop

#https://rdrr.io/cran/rstiefel/man/rustiefel.html
# eigen_list
# n_list -> np^2, np^3.... np^3
# For n in n_list
#     Sample L, V ~ otrhognal,
#     Mutiply L column by (p.already p^2, pgrids = p+1, ...p+p)
#     identify L*, V*
#     Compute eigen(L^TL^T)  

    






