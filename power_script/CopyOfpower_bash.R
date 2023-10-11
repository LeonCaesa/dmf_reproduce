setwd('/projectnb2/dmfgrp/dmf_revision/')
library(MASS)
library(tidyverse)
source('~/dmf/R/dmf.R')
source('/projectnb/dmfgrp/dmf_revision/power_util.R')


# set.seed(1)

#
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
# family_list[[1]] = gaussian()
family_list[[1]] = Gamma('log')
family_list[[2]] = binomial()
# family_list[[3]] = negative.binomial(phi)
family_list[[3]] = negative.binomial(5)
family_list[[4]] = poisson('log')
# name_list = c('Gaussian', 'Binom', 'Negbinom', 'Poisson')
name_list = c('Gamma', 'Binom', 'Negbinom', 'Poisson')

# weight_list = c(1,100,1,1)
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
    # if(true_family$family == 'binomial'){
    #   true_weight = matrix(rnegbin(n*p, true_weight, 1), nrow = n)
    #   true_weight[true_weight==0] = 1}

    V_star <- matrix(rnorm(q_star * p, 0, std_eta), nrow = p, ncol = q_star)
    L_star <- matrix(rnorm(n * q_star, 0, std_eta), nrow = n, ncol = q_star)
    eta_star <- tcrossprod(L_star, V_star) + mean_eta

    # [added to reset Gaussian variance]
    phi = phi_list[true_idx] # default dispersion for non-gaussian and non-negbinom
    # if(true_family$family == 'gaussian'){phi = phi * std_eta * 150}
    # if(grepl("Negative Binomial", true_family$family)){phi = 20}

    Y = generate_Y(true_family, eta_star, phi = phi, glm_weights = true_weight)
    if(true_family$family == 'gaussian'){Y = Y - min(Y)} # todo: positive gaussian?

    lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight, init_svd = TRUE, control = glm.control(maxit = 200))
    #lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight, control = glm.control(maxit = 200))
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


# #[test script]
set.seed(1)
n_list = c(600)
q_star = 5
std_eta= 0.1
mean_eta = 0
nrepeats = 1000
add_name = 'test'




family_list = c()
#family_list[[1]] = gaussian()
family_list[[1]] = Gamma('log')
family_list[[2]] = binomial()
# family_list[[3]] = negative.binomial(40)
family_list[[3]] = negative.binomial(5)
family_list[[4]] = poisson('log')
#name_list = c('Gaussian', 'Binom', 'Negbinom', 'Poisson')
name_list = c('Gamma', 'Binom', 'Negbinom', 'Poisson')
weight_list = c(1, 90, 1, 1)
domain_list = c('c', 'i', 'i', 'i')
phi_list = c(1, 1, 5, 1)

summary_table = data.frame(matrix(ncol = 12))
colnames(summary_table) = c('n', 'p', 'q_star', 'std_eta', 'mean_eta', 'G', 'repeats',
                            'add_name', 'true_family', 'comp_family', 'pvalue', 'chistat')

n = 500
p = 200
#for (n in n_list){
for (n in c(200)){
    p_list = n * seq(0.2, 0.7, by =0.1)
  # p_list = n * c(0.2, 0.7)
    for (p in c(80)){
   #for (p in p_list){

      # for( std_eta in seq(0.1, 0.4, by = 0.1)){
      for( std_eta in c(0.3)){
      #for( std_eta in seq(0.27,0.31, by = 0.01)){
       # for( std_eta in c(0.5, 1, 1.5)){

          for (repeats in 1:nrepeats){
          G = n*p/50
          for(true_idx in c(1,2,3,4)){
            true_family = family_list[[true_idx]]
            true_weight = weight_list[true_idx]


            V_star <- matrix(rnorm(q_star * p, 0, std_eta), nrow = p, ncol = q_star)
            L_star <- matrix(rnorm(n * q_star, 0, std_eta), nrow = n, ncol = q_star)
            eta_star <- tcrossprod(L_star, V_star) + mean_eta

            phi = phi_list[true_idx] # default dispersion for non-gaussian and non-negbinom
            # if(true_family$family == 'gaussian'){phi = phi * std_eta * 150}
            # if(true_family$family == 'Gamma'){phi = 1}

            Y = generate_Y(true_family, eta_star, phi = phi, glm_weights = true_weight)
            if(true_family$family == 'gaussian'){phi = 1; Y = Y - min(Y)} # todo: positive gaussian?

            # lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight, control = glm.control(maxit = 200))
            lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight, init_svd = TRUE, control = glm.control(maxit = 200))
            true_result = family_test(Y/true_weight, lv_true, G,  dispersion = phi, weights = true_weight)
            p_true = 1- pchisq(true_result$chisq_stat, length(true_result$norm_vec)-1)
            #print(p_true)

            summary_table[nrow(summary_table)+1, ] = c(n, p, q_star, std_eta, mean_eta, length(true_result$norm_vec), repeats, add_name,
                                                       name_list[true_idx], name_list[true_idx], p_true, true_result$chisq_stat)

            # for(esti_idx in c(1,2,3,4)){
            # # for(esti_idx in 3){
            #   #esti_idx = 3
            #   estimate_family <- family_list[[esti_idx]]
            #   correct_Y = correct_domain(Y, domain_list[esti_idx])
            #   if(true_family$family != 'binomial' && estimate_family$family == 'binomial'){correct_weight = max(correct_Y) }else{correct_weight = true_weight}

              # if (estimate_family$family != true_family$family){
              #   lv_estimate = dmf(correct_Y/correct_weight, estimate_family, rank = q_star, weights = correct_weight,  control = glm.control(maxit = 200))
              #   estimate_result = family_test(correct_Y/correct_weight, lv_estimate, G, weights = correct_weight)
              #   p_estimate = 1- pchisq(estimate_result$chisq_stat, length(estimate_result$norm_vec) - 1)#todo: adjust number of G
              #
              #   summary_table[nrow(summary_table)+1, ] = c(n, p, q_star, std_eta, mean_eta, length(estimate_result$norm_vec),
              #                                              repeats, add_name, name_list[true_idx], name_list[esti_idx], p_estimate)
              #
              #     }
            }#end of esti_idx
            # print(c(n,p,std_eta, phi, repeats))
          # }# end of true_idx
          print(c(n,p,std_eta, phi, repeats))
        } #end of repeats
      } # end of std_eta
    } # end of p
}# end of n
summary_table = summary_table[-1,]
summary_table$pvalue =round( as.numeric(summary_table$pvalue), 2)
print(summary_table)
#
# summary_table$n = as.numeric(summary_table$n)
# summary_table$p = as.numeric(summary_table$p)
# summary_table$pvalue = as.numeric(summary_table$pvalue)

summary_table$chistat = as.numeric(summary_table$chistat)
summary_table$pvalue = 1- mapply(pchisq, summary_table$chista, df = as.numeric(summary_table$G) - 1)
summary_table$pvalue = round( as.numeric(summary_table$pvalue), 2)
#qqplot(summary_table$chistat, rchisq(length(summary_table$chistat), as.numeric(summary_table$G)) -5)

summary_table$power = summary_table$pvalue<=0.05
summary_table %>% group_by(true_family, comp_family) %>% summarise(p_ = mean(power))



#
#
# ggplot(filter(summary_table)) + geom_boxplot(aes(x =  as.factor(n * p) ,
#   y = pvalue, colour = comp_family)) + facet_wrap(~std_eta + true_family)
#
# ggplot(filter(summary_table, n ==600, p ==240)) + geom_boxplot(aes(x =  as.factor(std_eta) ,
#                                                  y = pvalue, colour = comp_family)) + facet_wrap(~ true_family)

# #
# sigma2_list = seq(0, 1000, by =100)
# variance_list = rep(0, length(sigma2_list))
# for (sigma2_idx in 1:length(sigma2_list)){
#   sigma2 = sigma2_list[sigma2_idx]
#   variance_list[sigma2_idx] = 3 * sqrt(sigma2) + 9/20 * sigma2
# }
# plot(sigma2_list, sigma2_list/variance_list)
