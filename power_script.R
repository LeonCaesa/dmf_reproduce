setwd('/projectnb2/dmfgrp/dmf_revision/')
source('~/dmf/R/dmf.R')
source('/projectnb/dmfgrp/dmf_revision/power_util.R')
library(reshape2)
library(tidyverse)
library(MASS)
library(rstiefel)
library(gridExtra)
#lintr::lint('thm

set.seed(1)
phi = 1
p = 20
q_star = 5
#power_= 1
power_= 1.2
n = 1000
nrepeats = 50
G = 500


family_list = c()
family_list[[1]] = gaussian()
family_list[[2]] = Gamma('log')
family_list[[3]] = negative.binomial(phi)
family_list[[4]] = poisson('log')
family_list[[5]] = poisson('sqrt')
family_list[[6]] = binomial()
family_list[[7]] = binomial('cloglog')

domain_list = c('c', 'c', 'i', 'i', 'i', 'i', 'i')

# offset_list = c(5, 5, 5, 5, 5, 0, 0)
offset_list = c(0, 0, 0, 0, 0, 0, 0)

# working desired
# offset_list = c(0, # gaussian
#                 5, # Gamma # worked
#                 0, # Negbinom
#                 4, # Poisson(log) #worked except negbin
#                 5, # Poisson(sqrt) # worked
#                 0, # Binomial(logit)
#                 0  # Binomial(cloglog)
# )

weights_list = c(1, 1, 1, 1, 1,100, 100)

name_list = c('Gaussian', 'Gamma(log)', 'Negbinom(log)', 'Poisson(log)', 'Poisson(sqrt)', 'Binomial(logit)',  'Binomial(cloglog)')

# [dmf fit]
summary_table = data.frame(matrix(ncol = 7))


for (repeats in 1:nrepeats){
  for(true_idx in 1:7){
      true_family <- family_list[[true_idx]]
      true_weight <- weights_list[true_idx]
      if(true_family$family == 'binomial'){
        true_weight = matrix(rnegbin(n*p, true_weight, theta = 50), nrow = n)
        true_weight[true_weight==0] = 1}

      if(true_family$family == 'gaussian'){next}

      # V_star <- matrix(rnorm(q_star * p, 0, sqrt(0.005/q_star)), nrow = p, ncol = q_star)
      # L_star <- matrix(rnorm(n * q_star, 1, sqrt(0.005/q_star)), nrow = n, ncol = q_star)
      # eta_star <- tcrossprod(L_star, V_star) + offset_list[true_idx]
      param_star = switch (true_family$family,
        # Gamma = gene_trueLV(q_star, n, p, power_ = 0.5),
        # binomial = gene_trueLV(q_star, n, p, power_ = 1.2),
        # poisson = gene_trueLV(q_star, n, p, power_ = 0.5),
        gene_trueLV(q_star, n, p, power_)
      )

      eta_star = tcrossprod(param_star$L, param_star$V) + offset_list[true_idx]

      Y = generate_Y(true_family, eta_star, phi = phi, glm_weights = true_weight)



      lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight)
      true_result = family_test(Y/true_weight, lv_true, G, dispersion = phi, weights = true_weight)
      #qqplot(true_result$norm_vec, rnorm(length(true_result$norm_vec)))
      p_true = 1- pchisq(true_result$chisq_stat, G-1)


      for (estimate_idx in 1:7){
          estimate_family <- family_list[[estimate_idx]]
          correct_Y = correct_domain(Y, domain_list[estimate_idx])
          if(true_family$family != 'binomial'){correct_weight = max(correct_Y)
          }else{correct_weight = true_weight}

          if ((estimate_family$family != true_family$family) | ((estimate_family$link != true_family$link))){
            lv_estimate <- switch(estimate_family$family,
              binomial = dmf(correct_Y/correct_weight, estimate_family, rank = q_star, weights = correct_weight),
              dmf(correct_Y, estimate_family, rank = q_star, weights = 1)
            )
            estimate_result <- switch(estimate_family$family,
              binomial = family_test(correct_Y/correct_weight, lv_estimate, G , weights = correct_weight),
              family_test(correct_Y, lv_estimate, G , weights = 1)
              )
            p_estimate = 1- pchisq(estimate_result$chisq_stat, G-1)

            summary_table[nrow(summary_table) + 1,] = c(name_list[true_idx],
                                                        name_list[estimate_idx],
                                                         p_true, p_estimate,
                                                        G, offset_list[true_idx], repeats)

            print(c(paste(true_family$family, true_family$link, sep = '-'),
                    paste(estimate_family$family, estimate_family$link, sep = '-')
                    ))
            print(round(c(p_true, p_estimate, repeats), 4))
          } # end of self exclusion
      } # end of estimate idx
  }# end of true idx
}# end of repeats



save_name = '/projectnb/dmfgrp/dmf_revision/family_results/power_result6.RData'
summary_table = summary_table[-1,]
colnames(summary_table) = c('True', 'Estimate', 'P_true', 'P_esti', 'G', 'offset', 'Repeats')
save(summary_table, file = save_name)



# [test script]
setwd('/projectnb2/dmfgrp/dmf_revision/')
source('/projectnb/dmfgrp/dmf_revision/power_util.R')
source('~/dmf/R/dmf.R')
library(reshape2)
library(tidyverse)
library(MASS)
library(rstiefel)
library(reshape2)
library(gridExtra)
domain_list = c('c', 'c', 'i', 'i', 'i', 'i', 'i')

#set.seed(1)
phi = 1
p = 20
q_star = 5
phi = 1
power_ = 0.5
n = 1000
family_list = c()
family_list[[1]] = gaussian('log')
family_list[[2]] = Gamma('log')
family_list[[3]] = negative.binomial(phi)
family_list[[4]] = poisson('log')
family_list[[5]] = poisson('sqrt')
family_list[[6]] = binomial()
family_list[[7]] = binomial('cloglog')
# true_idx = 4; estimate_idx = 5; # possion(log) versus poisson(sqrt)
true_idx = 4; estimate_idx = 7; # possion(log) versus binomial(cloglog)
# true_idx = 4; estimate_idx = 2; # possion(log) versus Gamma(log)
# true_idx = 5; estimate_idx = 4; # possion(sqrt) versus poisson(log)
# true_idx = 5; estimate_idx = 7; # possion(sqrt) versus binomial(cloglog)
# true_idx = 2; estimate_idx = 3; # Gamma(log) versus Negbinom(log)
# true_idx = 2; estimate_idx = 4; # Gamma(log) versus poisson(log)
# true_idx = 7; estimate_idx = 6; # binomial(cloglog) versus binomial(logit)
# true_idx = 3; estimate_idx = 4; # negbin(log) versus poisson(log)

offset_test = 4; true_weight = 1; G = 500
true_family <- family_list[[true_idx]]
estimate_family <- family_list[[estimate_idx]]

if(true_family$family == 'binomial'){
  true_weight = matrix(rnegbin(n*p, true_weight, 1), nrow = n)
  true_weight[true_weight==0] = 1}

param_star = gene_trueLV(q_star, n, p, power_)
eta_star = tcrossprod(param_star$L, param_star$V) + offset_test


Y = generate_Y(true_family, eta_star, phi = phi, glm_weights = true_weight)


lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight)
true_result = family_test(Y/true_weight, lv_true, G, dispersion = phi, weights = true_weight)
qqplot(true_result$norm_vec, rnorm(length(true_result$norm_vec)))
p_true = 1- pchisq(true_result$chisq_stat, G-1)

correct_Y = correct_domain(Y, domain_list[estimate_idx])
if(true_family$family != 'binomial'){correct_weight = max(correct_Y)
}else{correct_weight = true_weight}

if ((estimate_family$family != true_family$family) | ((estimate_family$link != true_family$link))){
  lv_estimate <- switch(estimate_family$family,
                        binomial = dmf(correct_Y/correct_weight, estimate_family, rank = q_star, weights = correct_weight),
                        dmf(correct_Y, estimate_family, rank = q_star, weights = 1)
  )
  estimate_result <- switch(estimate_family$family,
                            binomial = family_test(correct_Y/correct_weight, lv_estimate, G , weights = correct_weight),
                            family_test(correct_Y, lv_estimate, G , weights = 1)
  )
  p_estimate = 1- pchisq(estimate_result$chisq_stat, G-1)
}
print(c(p_true, p_estimate))



