setwd('/projectnb2/dmfgrp/dmf_revision/')
library(MASS)
library(tidyverse)
source('./dmf/R/dmf.R')
source('./power_util.R')


# set.seed(1)
# n = 100
# #p = 20
# p = 50
# q_star = 5
# std_eta= 0.1
# mean_eta =0
# G = 200
# nrepeats = 1
# add_name = 'test'
argv <- commandArgs(TRUE)
if (length(argv) > 0){
  n <- as.numeric( argv[1] )
  p <- as.numeric( argv[2] )
  q_star <- as.numeric( argv[3] )
  std_eta  <- as.numeric( argv[4] )
  mean_eta  <- as.numeric( argv[5] )
  G  <- as.numeric( argv[6] )
  nrepeats <- as.numeric( argv[7] )
  add_name <- as.character(argv[8])
}

# [Logging the experiment info]
paste("n:",n)
paste("p:",p)
paste("q_star:",q_star)
paste("std_eta:",std_eta)
paste("mean_eta:",mean_eta)
paste("G:", G)
paste("Number of Simulations:",nrepeats)
paste("add_name:", add_name)


save_dir = '/projectnb/dmfgrp/dmf_revision/power_script/power_result'
save_name = paste('n', n,
                  '_p', p,
                  '_q', q_star,
                  '_sd', std_eta,
                  '_mu', mean_eta,
                  '_G', G,
                  '_S', nrepeats,
                  add_name,'.RData', sep = '')
save_file = paste(save_dir, save_name, sep = '/')




phi = 1
family_list = c()
family_list[[1]] = Gamma('log')
family_list[[2]] = binomial()
family_list[[3]] = negative.binomial(5)
family_list[[4]] = poisson('log')
name_list = c('Gamma', 'Binom', 'Negbinom', 'Poisson')


weight_list = c(1,90,1,1)
phi_list = c(1, 1, 5, 1)
domain_list = c('c', 'i', 'i', 'i')



summary_table = data.frame(matrix(ncol = 11))
colnames(summary_table) = c('n', 'p', 'q_star', 'std_eta', 'mean_eta', 'G', 'repeats',
                            'add_name', 'true_family', 'comp_family', 'pvalue')

# [simulate data]
for (repeats in 1:nrepeats){
  for(true_idx in 1:4){

    true_family = family_list[[true_idx]]
    true_weight = weight_list[true_idx]
    V_star <- matrix(rnorm(q_star * p, 0, std_eta), nrow = p, ncol = q_star)
    L_star <- matrix(rnorm(n * q_star, 0, std_eta), nrow = n, ncol = q_star)
    eta_star <- tcrossprod(L_star, V_star) + mean_eta


    phi = phi_list[true_idx] # default dispersion for non-gaussian and non-negbinom
    Y = generate_Y(true_family, eta_star, phi = phi, glm_weights = true_weight)

    lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight, init_svd = TRUE, control = glm.control(maxit = 200))
    true_result = family_test(Y/true_weight, lv_true, G, dispersion = phi, weights = true_weight)
    p_true = 1- pchisq(true_result$chisq_stat, length(true_result$norm_vec)-1)

    summary_table[nrow(summary_table)+1, ] = c(n, p, q_star, std_eta, mean_eta, length(true_result$norm_vec), repeats, add_name,
                                               name_list[true_idx], name_list[true_idx], p_true)
    print(tail(summary_table,1))
    for(esti_idx in 1:4){
      estimate_family <- family_list[[esti_idx]]
      correct_Y = correct_domain(Y, domain_list[esti_idx])
      if(true_family$family != 'binomial' && estimate_family$family == 'binomial'){correct_weight = max(correct_Y)}else{correct_weight = true_weight}

      if (estimate_family$family != true_family$family){
        lv_estimate = dmf(correct_Y/correct_weight, estimate_family, rank = q_star, weights = correct_weight, control = glm.control(maxit = 200))
        estimate_result = family_test(correct_Y/correct_weight, lv_estimate, G , weights = correct_weight)
        p_estimate = 1- pchisq(estimate_result$chisq_stat, length(estimate_result$norm_vec) - 1)#todo: adjust number of G

        summary_table[nrow(summary_table)+1, ] = c(n, p, q_star, std_eta, mean_eta, length(estimate_result$norm_vec),
                                                   repeats, add_name, name_list[true_idx], name_list[esti_idx], p_estimate)

        print(tail(summary_table,1))
      }
    }#end of esti_idx
  }# end of true_idx
}# end of repeats
summary_table = summary_table[-1,]
save(summary_table, file = save_file)