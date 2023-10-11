if(!exists("foo", mode="function")) source("../dmf/R/dmf.R")
if(!exists("foo", mode="function")) source('./casestudy_utils.R')

library(MASS)
library(NMF)

data_dir = '/projectnb2/dmfgrp/dmf_revision/data'
q_hat = 20



# [leukemia]
data_name = 'leukemia'
factor_families = c();
factor_families[[1]] = negative.binomial(0.1);
factor_families[[2]] = poisson('identity');
factor_families[[3]] = poisson();
weights_list = c(1,1,1) 
save_L  = FALSE
save_dir = paste('../data/', data_name, sep = '')
load(paste(data_dir, paste(data_name, '.RData', sep = ''), sep= '/'));
load(paste(data_dir, paste(data_name, '_label.RData', sep = ''), sep= '/'));

for (idx in 1:length(factor_families)){
  glm_family = factor_families[[idx]]
  glm_weight = weights_list[idx]
  if (glm_family$link == 'identity' & glm_family$family ==  "poisson"){
    fit_nmf = nmf(Y + 0.001, q_hat, nrun=5, seed = 'nndsvd')
    dmf_fit = c();dmf_fit$V = t(fit_nmf@fit@H); dmf_fit$L = fit_nmf@fit@W; dmf_fit$family = glm_family
    glm_family$family = 'nmf'
  }else{dmf_fit = dmf(Y/glm_weight, glm_family, q_hat, weights = glm_weight)  }
  dmf_identified = center_identify(dmf_fit)
  if (save_L){save_df = data.frame(cbind(dmf_identified$L, y))}else{save_df = data.frame(cbind(dmf_identified$V, y))}
  save(save_df, file = paste(save_dir, paste(glm_family$family, '_q', q_hat, '.RData', sep = ''), sep = '/')) 
}




# [cbcl]
data_name = 'cbcl'
factor_families = c();
factor_families[[1]] = binomial(); factor_families[[2]] = negative.binomial(0.462);
factor_families[[3]] = poisson('identity'); factor_families[[4]] = poisson();
weights_list = c(255,1,1,1) 
save_L  = TRUE
save_dir = paste('../data/', data_name, sep = '')
load(paste(data_dir, paste(data_name, '.RData', sep = ''), sep= '/'));load(paste(data_dir, paste(data_name, '_label.RData', sep = ''), sep= '/'));

for (idx in 1:length(factor_families)){
  glm_family = factor_families[[idx]]
  glm_weight = weights_list[idx]
  if (glm_family$link == 'identity' & glm_family$family ==  "poisson"){
      fit_nmf = nmf(Y + 0.001, q_hat, nrun=5, seed = 'nndsvd')
      dmf_fit = c();dmf_fit$V = t(fit_nmf@fit@H); dmf_fit$L = fit_nmf@fit@W; dmf_fit$family = glm_family
      glm_family$family = 'nmf'
  }else{dmf_fit = dmf(Y/glm_weight, glm_family, q_hat, weights = glm_weight)  }
  dmf_identified = center_identify(dmf_fit)
  if (save_L){save_df = data.frame(cbind(dmf_identified$L, y))}else{save_df = data.frame(cbind(dmf_identified$V, y))}
  save(save_df, file = paste(save_dir, paste(glm_family$family, '_q', q_hat, '.RData', sep = ''), sep = '/')) 
}



# [network]
data_name = 'network'
factor_families = c();
factor_families[[1]] = binomial();
weights_list = c(1) 
save_L  = TRUE

save_dir = paste('../data/', data_name, sep = '')
load(paste(data_dir, paste(data_name, '.RData', sep = ''), sep= '/'));load(paste(data_dir, paste(data_name, '_label.RData', sep = ''), sep= '/'));

for (idx in 1:length(factor_families)){
  glm_family = factor_families[[idx]]
  glm_weight = weights_list[idx]
  if (glm_family$link == 'identity' & glm_family$family ==  "poisson"){
    fit_nmf = nmf(Y + 0.001, q_hat, nrun=5, seed = 'nndsvd')  
    dmf_fit = c();dmf_fit$V = t(fit_nmf@fit@H); dmf_fit$L = fit_nmf@fit@W; dmf_fit$family = glm_family
    glm_family$family = 'nmf'
  }else{dmf_fit = dmf(Y/glm_weight, glm_family, q_hat, weights = glm_weight)  }
  dmf_identified = center_identify(dmf_fit)
  if (save_L){save_df = data.frame(cbind(dmf_identified$L, y))}else{save_df = data.frame(cbind(dmf_identified$V, y))}
  save(save_df, file = paste(save_dir, paste(glm_family$family, '_q', q_hat, '.RData', sep = ''), sep = '/')) 
}
D = diag(rowSums(Y))
L = D-Y
D2 = diag(1/sqrt(rowSums(Y)))
L_sym = tcrossprod(crossprod(D2,L), D2)
eigen_results = eigen(L_sym)
save_df = data.frame(cbind(eigen_results$vectors[,(ncol(eigen_results$vectors)-q_hat+1):ncol(eigen_results$vectors)]),y)
save(save_df, file = paste(save_dir, paste('spectrum', '_q', q_hat, '.RData', sep = ''), sep = '/')) 



# [nlp]
data_name = 'nlp'
factor_families = c();
factor_families[[1]] = negative.binomial(0.006);
factor_families[[2]] = poisson('identity'); factor_families[[3]] = poisson();
weights_list = c(1,1,1) 
save_L  = TRUE

save_dir = paste('../data/', data_name, sep = '')
load(paste(data_dir, paste(data_name, '.RData', sep = ''), sep= '/'));load(paste(data_dir, paste(data_name, '_label.RData', sep = ''), sep= '/'));
for (idx in 1:length(factor_families)){
  glm_family = factor_families[[idx]]
  glm_weight = weights_list[idx]
  if (glm_family$link == 'identity' & glm_family$family ==  "poisson"){
    fit_nmf = nmf(Y+0.001, q_hat, nrun=5, seed = 'nndsvd')  
    dmf_fit = c();dmf_fit$V = t(fit_nmf@fit@H); dmf_fit$L = fit_nmf@fit@W; dmf_fit$family = glm_family
    glm_family$family = 'nmf'
  }else{dmf_fit = dmf(Y/glm_weight, glm_family, q_hat, weights = glm_weight)}
  dmf_identified = center_identify(dmf_fit)
  if (save_L){save_df = data.frame(cbind(dmf_identified$L, y))}else{save_df = data.frame(cbind(dmf_identified$V, y))}
  save(save_df, file = paste(save_dir, paste(glm_family$family, '_q', q_hat, '.RData', sep = ''), sep = '/')) 
}


