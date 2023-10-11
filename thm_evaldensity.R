setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)
result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/lowerCvgThd/'

n_runs = 10
family_list = c()
#phi = 3
#phi = 50
phi = 6
family_list[[1]] = negative.binomial(phi)
family_list[[2]] = poisson()
family_list[[3]] = gaussian()


# [eval results with random sampled residuals]
p_lists =  seq(10, 40, by = 10)
#p_lists =  seq(100, 400, by = 100)
n_lists =  seq(1000, 5000, by =1000)




density_size = 512 * 2 

residL_df <- data.frame(matrix(ncol = density_size +4, nrow = 0))
colnames(residL_df ) <- c(paste0("sample", c(1:density_size)), 'n', 'p', 'run', 'family')
residV_df <- data.frame(matrix(ncol =density_size + 4, nrow = 0))
colnames(residV_df ) <- c(paste0("sample", c(1:density_size)), 'n', 'p', 'run', 'family')

resideta_df <- data.frame(matrix(ncol = density_size + 4, nrow = 0))
colnames(resideta_df ) <- c(paste0("sample", c(1:density_size)), 'n', 'p', 'run', 'family')



for (glm_family in family_list){
  family_dir = paste(result_dir, glm_family$family, sep = '/')
  for (n in n_lists){
    for(p in p_lists){
      # create model dir
      np_dir = paste(family_dir, paste(n, p, sep = '_'), sep = '/')
      
      
      for (run in 1:n_runs){
        param_savedir = paste(np_dir, paste('run',run, '_param_star.RData', sep = ''), sep = '/')
        dmf_savedir = paste(np_dir, paste('run',run, '_dmf.RData', sep = ''), sep = '/')
        
        # [load model from dir]
        load(param_savedir)
        load(dmf_savedir)
        # [compute the residuals]
        resid_L= c(param_star$L - dmf_result$L)
        resid_V = c(param_star$V - dmf_result$V)
        eta_star = tcrossprod(param_star$L, param_star$V)
        eta_hat = tcrossprod(dmf_result$L, dmf_result$V)
        resid_eta = c(eta_star - eta_hat)
        
        
        addL_rowidx = nrow(residL_df)+1
        L_DenEsti = density(resid_L, density_size)
        residL_df[addL_rowidx, ] <- c(L_DenEsti$x, L_DenEsti$y ,n, p , run, glm_family$family)
        
        addV_rowidx = nrow(residV_df)+1
        V_DenEsti = density(resid_V, density_size)
        residV_df[addV_rowidx, ] <- c(V_DenEsti$x, V_DenEsti$y, n, p , run, glm_family$family)
        
        addeta_rowidx = nrow(resideta_df)+1
        eta_DenEsti = density(resid_eta, density_size)
        resideta_df[addeta_rowidx, ] <- c(eta_DenEsti$x, eta_DenEsti$y, n, p , run, glm_family$family)
        
      }# end of n runs
    }# end of p loop
  }# end of n loop
}# end of family loop

library(tidyverse)

resideta_long <-  resideta_df %>%gather(resid_x, value_x, sample1:sample512)
resideta_long <-  resideta_long %>%gather(resid_y, value_y, sample513:sample1024)
resideta_long$value_x = as.numeric(resideta_long$value_x)
resideta_long$value_y = as.numeric(resideta_long$value_y)


residL_long <-  residL_df %>%gather(resid_x, value_x,  sample1:sample512)
residL_long <-  residL_long %>%gather(resid_y, value_y, sample513:sample1024)
residL_long$value_x = as.numeric(residL_long$value_x)
residL_long$value_y = as.numeric(residL_long$value_y)

residV_long <-  residV_df %>%gather(resid_x, value_x,  sample1:sample512)
residV_long <-  residV_long %>%gather(resid_y, value_y,  sample513:sample1024)
residV_long$value_x = as.numeric(residV_long$value_x)
residV_long$value_y = as.numeric(residV_long$value_y)

ggplot(residV_long) + geom_point(aes(x = resid_x, y = resid_y, colour = n)) + facet_wrap(~family +p, nrow = length(family_list))
ggplot(residL_long) + geom_points(aes(x = resid_x, y = resid_y, colour = n)) + facet_wrap(~family +p, nrow = length(family_list))
ggplot(resideta_long) + geom_points(aes(x = resid_x, y = resid_y, colour = n)) + facet_wrap(~family +p, nrow = length(family_list))


# 
# base_family = c(family_list[[1]]$family)
# base_p = c(100)
# base_n = c(5000)
# range_density = filter(residL_long, n %in% base_n, p %in% base_p, family %in% base_family )
# l_lim <- range(density(range_density$value)$x)
# range_density = filter(residV_long, n %in% base_n, p %in% base_p, family %in% base_family )
# v_lim <- range(density(range_density$value)$x)
# range_density = filter(resideta_long, n %in% base_n, p %in% base_p, family %in% base_family )
# eta_lim <- range(density(range_density$value)$x)

# l_lim = c(-100,100)
# v_lim = c(-0.1,0.1)
# eta_lim = c(-3,3)

l_lim = range(residL_long$value)
v_lim = range(residV_long$value)
eta_lim =  range(resideta_long$value)

ggplot(residL_long) + geom_density(aes(x = value, colour = n)) + 
  xlim(l_lim) + facet_wrap(~ family +p, nrow = length(family_list)) 
ggplot(residV_long) + geom_density(aes(x = value, colour = n)) + 
  xlim(v_lim) + facet_wrap(~family +p, nrow = length(family_list))

ggplot(resideta_long) + geom_density(aes(x = value, colour = n)) + 
  xlim(range(leta$x)) + facet_wrap(~family +p, nrow = length(family_list))

# [some investigation]
# 
# eta_star = tcrossprod(param_star$L_star, param_star$V_star)
# eta_hat =  tcrossprod(dmf_result$L, dmf_result$V)
# plot(density(eta_star - eta_hat))
# qqplot(param_star$L_star, dmf_result$L)
# qqplot(param_star$V_star, dmf_result$V)
# qqplot(eta_star,  eta_hat)
# 
# q_star = 6
# 
# slv_star = svd(eta_star, nu = q_star, nv = q_star)
# negative_flag = c(1:q_star)[slv_star$u[1,]<0]
# slv_star$u[,negative_flag] = matrix(-slv_star$u[,negative_flag])
# slv_star$v[,negative_flag] = matrix(-slv_star$v[,negative_flag])
# 
# plot(density(slv_star$v - param_star$V_star))
# 
# 
