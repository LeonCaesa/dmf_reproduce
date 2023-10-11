setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(reshape2)
library(tidyverse)
library(MASS)
library(rstiefel)
library(lintr)
library(reshape2)
library(gridExtra)

phi = 1
p = 200
q_star = 5
power_= 1
n = 2000
G = 1000
nrun = 200

gof_simu <- function(p, q_star, n, phi, G, glm_family, power_ = 1){

  param_star = gene_trueLV(q_star, n, p, power_)
  eta_star = tcrossprod(param_star$L, param_star$V)
  Y <- generate_Y(glm_family, eta_star, glm_weights = 0, phi = phi)

  if(grepl('^Negative Binomial', glm_family$family)){
    phi_estimate = mean(Y)^2/(sd(Y)^2 - mean(Y))
    glm_family_estimate = negative.binomial(phi_estimate)
    glm_family_naive = poisson()
  }else if(glm_family$family == 'Gamma'){
    phi_estimate = mean(Y)^2/(sd(Y)^2)
    glm_family_estimate = glm_family
    glm_family_naive = gaussian()
  }else if(glm_family$family == 'poisson'){
    phi_estimate = mean(Y)^2/(sd(Y)^2)
    glm_family_estimate = negative.binomial(phi_estimate)
    glm_family_naive = negative.binomial(1e5)}
  # [dmf fit]
  lv_true = dmf(Y, glm_family, rank = q_star)
  lv_estimate = dmf(Y, glm_family_estimate, rank = q_star)
  lv_naive = dmf(Y, glm_family_naive, rank = q_star)

  # [residual compute]
  estimate_result = family_test(Y, lv_estimate, G)
  true_result = family_test(Y, lv_true, G)
  naive_result = family_test(Y, lv_naive, G)

  # [p value calcs]
  p_estimate = 1- pchisq(estimate_result$chisq_stat, G-1)
  p_true = 1- pchisq(true_result$chisq_stat, G-1)
  p_naive = 1- pchisq(naive_result$chisq_stat, G-1)


  # [plotting]
  norm_vec1 = true_result$norm_vec
  norm_vec2 = estimate_result$norm_vec
  norm_vec3 = naive_result$norm_vec

  norm_df = data.frame(cbind(norm_vec1, norm_vec2, norm_vec3))
  colnames(norm_df) = c('P(NegBinom True Phi)', 'P(NegBinom Esti Phi)', 'P(Poisson)')
  long_norm_df = melt(norm_df)
  g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) +
    geom_density(alpha = 0.5)+
    scale_fill_brewer(palette="Dark2") +
    theme_classic() +
    stat_function(fun = dnorm, aes(colour = 'Standard Normal'), linetype = "dashed")+
    labs(y="Denstiy Value", x = "Y") + xlim(-5,5)
  g1 = ggplot(long_norm_df) +
    stat_qq(aes(sample = value , color = variable)) +
    scale_color_brewer(palette="Dark2")+ theme_classic() +
    geom_abline(slope =1, intercept =0)+
    labs(y="Quantile(Test)", x = "Quantile(Norm)")

  grid.arrange(g1,g2, ncol =2)

  return(data.frame(p_true = p_true, p_estimate = p_estimate, p_naive = p_naive))
}
glm_family = negative.binomial(phi)
gof_simu(p, q_star, n, phi, G, glm_family, power_ = 1)




p_result <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(p_result) <- c('p_true', 'p_estimate', 'p_naive')
family_list = c()
family_list[[1]] = negative.binomial(phi)
family_list[[2]] = poisson()

result_dir = '/projectnb/dmfgrp/dmf_revision/family_results/'

for (glm_family in family_list){
    for (i in 1:nrun){
    print(i)
    p_result[nrow(p_result)+1, ] = gof_simu(p, q_star, n, phi = 1, G, glm_family)
    }
    folder_dir = paste(result_dir, glm_family$family, sep = '/')
    dir.create(file.path(folder_dir), showWarnings = FALSE)
    file_name = paste(paste(n,p, phi, power_, G, sep = '_'), '.RData', sep = '')
    save(p_result, file = paste(folder_dir, file_name, sep = '/'))
}

