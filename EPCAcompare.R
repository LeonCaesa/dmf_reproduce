setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)
library(rstiefel)


n_lists =  seq(200, 1000, by = 200)
p_lists =  seq(100, 200, by = 20)


phi = 5
family_list = c()
family_list[[1]] = binomial() 
 
# family_list[[1]] = negative.binomial(phi)
# family_list[[2]] = poisson()
# family_list[[3]] = gaussian()
# family_list[[4]] = binomial() 

q_star <- 6
n_runs = 10
power_ = 1
glm_weights = 10


set.seed(1)
time_df <- data.frame( n= numeric(0), p = numeric(0), 
                      family = character(0), model = character(0),
                      time = numeric(0), 
                      repeats = numeric())


result_dir = '/projectnb/dmfgrp/dmf_revision/comp_results'
for (glm_family in family_list){
  family_dir = paste(result_dir, glm_family$family, sep = '/')
  dir.create(file.path(family_dir), showWarnings = FALSE)
  for (n in n_lists){
    for(p in p_lists){
      # create model dir
      np_dir = paste(family_dir, paste(n, p, sep = '_'), sep = '/')

      if (!file.exists(np_dir)){
        dir.create(file.path(np_dir), showWarnings = FALSE)
      for (run in 1:n_runs){
        #dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
        
        param_star = gene_trueLV(q_star, n, p, power_)
        eta_star = tcrossprod(param_star$L_star,param_star$V_star)
        Y <- generate_Y(glm_family, eta_star, glm_weights = glm_weights, phi = phi)
        param_star$Y = Y/glm_weights
        
        # [DMF fit]
        start.time <- Sys.time()
        dmf_result = dmf(param_star$Y, family = glm_family, rank = q_star, weights =  glm_weights)
        end.time <- Sys.time()
        time_diff = as.numeric(end.time-start.time)
        time_df[nrow(time_df)+1, ] <- c(n, p, glm_family$family, 'dmf', time_diff, run)

        # [EPCA fit]
        start.time <- Sys.time()
        irls_result = IRLS_Matrix(param_star$Y, q_star, glm_family, glm_weights = matrix(rep(glm_weights, n * p), nrow = n), verbose = FALSE)
        end.time <- Sys.time()
        time_diff = as.numeric(end.time-start.time)
        time_df[nrow(time_df)+1, ] <- c(n, p, glm_family$family, 'EPCA', time_diff, run)
        
        print(tail(time_df,2))

        # [Save fitted result]
        param_savedir = paste(np_dir, paste('run',run, '_param_star.RData', sep = ''), sep = '/')
        dmf_savedir = paste(np_dir, paste('run',run, '_dmf.RData', sep = ''), sep = '/')
        irls_savedir = paste(np_dir, paste('run',run, '_irls.RData', sep = ''), sep = '/')
        
        save(param_star, file = param_savedir)
        save(dmf_result, file = dmf_savedir)
        save(irls_result, file = irls_savedir)
      }# end of n runs
        
    } # end of if argunent
      
    }# end of p loop
    
  }# end of n loop
  
}# end of family loop

# paste(np_dir, paste('run',run, '_param_star.RData', sep = ''), sep = '/')
# save(time_df, file =  "/projectnb/dmfgrp/dmf_revision/comp_results/time_BinomExp.RData")
