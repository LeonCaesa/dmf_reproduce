setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(reshape2)
library(tidyverse)
library(MASS)
library(rstiefel)
library(lintr)
library(reshape2)
library(gridExtra)
#lintr::lint('thm_facilitate.r')

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
#family_list[[1]] = poisson()

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



##### [end of the script]
#
# glm_family = negative.binomial(phi)
# param_star = gene_trueLV(q_star, n, p, power_)
# eta_star = tcrossprod(param_star$L_star,param_star$V_star)
# Y <- generate_Y(glm_family, eta_star, glm_weights = 1, phi = phi)
#
#
# # W_star <- matrix(rep(1, q_star*p), nrow = q_star, ncol = p)
# # Z_star <- matrix(rep(0.5, q_star*p), nrow = n, ncol = q_star)
# # eta_star <- crossprod(t(Z_star), W_star)
# # L_star = dmf_identify(eta_star,q_star)$L_star
# # V_star = dmf_identify(eta_star,q_star)$V_star
#
# # mu_star <- glm_family$linkinv(eta_star)
# # glm_weights = matrix(rep(1,n*p), nrow=n, ncol = p)
# # Y <- generate_Y(glm_family, eta_star, glm_weights = glm_weights, phi = phi)
#
#
# if(grepl('^Negative Binomial', glm_family$family)){
#   phi_estimate = mean(Y)^2/(sd(Y)^2 - mean(Y))
#   glm_family_estimate = negative.binomial(phi_estimate)
#   glm_family_naive = poisson()
# }else if(glm_family$family == 'Gamma'){
#   phi_estimate = mean(Y)^2/(sd(Y)^2)
#   glm_family_estimate = glm_family
#   glm_family_naive = gaussian()
# }else if(glm_family$family == 'poisson'){
#   phi_estimate = mean(Y)^2/(sd(Y)^2)
#   glm_family_estimate = negative.binomial(phi_estimate)
#   glm_family_naive = negative.binomial(1e5)}
#
# lv_true = dmf(Y, glm_family, rank = q_star)
# lv_estimate = dmf(Y, glm_family_estimate, rank = q_star)
# lv_naive = dmf(Y, glm_family_naive, rank = q_star)
#
# estimate_result = family_test(Y, lv_estimate, G)
# true_result = family_test(Y, lv_true, G)
# naive_result = family_test(Y, lv_naive, G)
#
#
# p_estimate = 1- pchisq(estimate_result$chisq_stat, G-1)
# p_true = 1- pchisq(true_result$chisq_stat, G-1)
# p_naive = 1- pchisq(naive_result$chisq_stat, G-1)
#
#
# lv_plot = lv_naive
# plot(density(Y - lv_plot$family$linkinv(tcrossprod(lv_plot$L, lv_plot$V))))
#
# norm_vec1 = true_result$norm_vec
# norm_vec2 = estimate_result$norm_vec
# norm_vec3 = naive_result$norm_vec
#
# norm_df = data.frame(cbind(norm_vec1, norm_vec2, norm_vec3))
# colnames(norm_df) = c('P(NegBinom True Phi)', 'P(NegBinom Esti Phi)', 'P(Poisson)')
#
#
# long_norm_df = melt(norm_df)
#
#
# g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) +
#   geom_density(alpha = 0.5)+
#   scale_fill_brewer(palette="Dark2") +
#   theme_classic() +
#   stat_function(fun = dnorm, aes(colour = 'Standard Normal'), linetype = "dashed")+
#   labs(y="Denstiy Value", x = "Y") + xlim(-5,5)
#
#
# g1 = ggplot(long_norm_df) +
#   stat_qq(aes(sample = value , color = variable)) +
#   scale_color_brewer(palette="Dark2")+ theme_classic() +
#   geom_abline(slope =1, intercept =0)+
#   labs(y="Quantile(Test)", x = "Quantile(Norm)")
#
#
# grid.arrange(g1,g2, ncol =2)
# print(data.frame(p_true = p_true, p_estimate = p_estimate, p_naive = p_naive))



# [script experiments]

#
#
# gof_simu <- function(p, q_star, n, phi, G, glm_family){
#
#     param_star = gene_trueLV(q_star, n, p, power_)
#     eta_star = tcrossprod(param_star$L_star,param_star$V_star)
#     Y <- generate_Y(glm_family$family, eta_star, glm_weights = 0, phi = phi)
#
#
#     if(grepl('^Negative Binomial', glm_family$family)){
#       phi_estimate = mean(Y)^2/(sd(Y)^2 - mean(Y))
#       glm_family_estimate = negative.binomial(phi_estimate)
#       glm_family2 = poisson()
#     }else if(glm_family$family == 'Gamma'){
#       phi_estimate = mean(Y)^2/(sd(Y)^2)
#       glm_family_estimate = glm_family
#       glm_family2 = gaussian()
#     }else if(glm_family$family == 'poisson'){
#       phi_estimate = mean(Y)^2/(sd(Y)^2)
#       glm_family_estimate = negative.binomial(phi_estimate)
#       glm_family2 = negative.binomial(1e5)}
#
#     lv_true = dmf(Y, glm_family, rank = q_star)
#     lv_estimate = dmf(Y, glm_family_estimate, rank = q_star)
#     lv_naive = dmf(Y, glm_family2, rank = q_star)
#
#
#     p_estimate = 1- pchisq(estimate_result$chisq_stat, G-1)
#     p_true = 1- pchisq(true_result$chisq_stat, G-1)
#     p_naive = 1- pchisq(naive_result$chisq_stat, G-1)
#
#
#     if(grepl('^Negative Binomial', glm_family$family)){
#       p_naive = 1- pchisq(family_test(Y, lv_naive, G, chisq_stat = TRUE), G-1)
#     }else{
#       p_naive = 1- pchisq(family_test(Y, lv_naive, G, chisq_stat = TRUE), G-1)}
#   return(data.frame(p_true = p_true, p_estimate = p_estimate, p_naive = p_naive))
# }
#
# nrun = 200
# # nrun = 10
# p_result <- data.frame(matrix(ncol = 3, nrow = 0))
# colnames(p_result) <- c('p_true', 'p_estimate', 'p_naive')
#
#
# family_list = c()
# family_list[[1]] = negative.binomial(phi)
# family_list[[2]] = poisson()
# #family_list[[1]] = poisson()
#
# result_dir = '/projectnb/dmfgrp/dmf_revision/family_results/'
#
# for (glm_family in family_list){
#     for (i in 1:nrun){
#       print(i)
#       p_result[nrow(p_result)+1, ] = gof_simu(p, q_star, n, phi = 1, G, glm_family)
#     }
#     folder_dir = paste(result_dir, glm_family$family, sep = '/')
#     dir.create(file.path(folder_dir), showWarnings = FALSE)
#     file_name = paste(paste(n,p, phi, power_, G, sep = '_'), '.RData', sep = '')
#     save(p_result, file = paste(folder_dir, file_name, sep = '/'))
# }


#[plotting and evaluation]

pvalue_df = p_result
#colnames(pvalue_df) = c('NegBinom (True Phi)', 'NegBinom (Esti Phi)', 'Poisson')
colnames(pvalue_df) = c('P(Poisson)', 'P(NegBinom Esti Phi)', 'P(NegBinom True Phi)')
long_pvalue = melt(pvalue_df)
load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/GOF_Experiments/Sensitivity_Analysis.RData")
ggplot(long_pvalue, aes(x= variable, y = value)) + geom_boxplot()

# ggplot(long_pvalue, aes(x=variable, y=value)) +
#   geom_violin(trim=FALSE, fill="gray")+ geom_boxplot(fill='white', width = 0.1)  + geom_point(size = 0.2)+ stat_summary(fun.y=median, geom="point", size=2, color="red")+
#   scale_fill_brewer(palette="Dark2")+ theme_classic() +ylab('P-value')+ylim(c(0,1))


png('/projectnb/dmfgrp/dmf_revision/figure/GOF_OPSensi.png', units="in", width=6, height=4, res=300)
long_pvalue = filter(long_pvalue, variable %in% c('P(Poisson)', 'P(NegBinom Esti Phi)'))

library(latex2exp)
#recode(long_pvalue$variable, p_true = TeX('Negbinom($phi$)'), p_estiamte = TeX('Negbinom$\\hat{phi}$'), p_naive = 'Poisson')


ggplot(long_pvalue, aes(x=variable, y=value)) +
  geom_boxplot(fill='white',width = 0.8)  + geom_point(size = 0.2, position = 'jitter')+ stat_summary(fun.y=median, geom="point", size= 2, color="red")+
  scale_fill_brewer(palette="Dark2")+ theme_light() +ylab('P-value') + ylim(c(0,1)) + xlab('DMF Fit Family (repeated 200 times)')
dev.off()


png('/projectnb/dmfgrp/dmf_revision/figure/gof_sensitivity2.png', units="in", width=6, height=4, res=300)
ggplot(long_pvalue, aes(x=variable, y=value)) +
  geom_boxplot(fill='white',width =0.8) + scale_x_discrete(labels=c("p_true" = TeX('Negbinom($phi$)'), 'p_estimate' = TeX('Negbinom ($\\hat{phi}$)'),
                                                                    "p_naive" = 'Poisson'))+ geom_point(size = 0.2, position = 'jitter')+
  stat_summary(fun.y=median, geom="point", size= 2, color="red") +
  scale_fill_brewer(palette="Dark2")+ theme_light(base_size = 12) +ylab('P-value')+ylim(c(0,1)) + xlab('DMF Fit Family (repeated 10k times)')
dev.off()
