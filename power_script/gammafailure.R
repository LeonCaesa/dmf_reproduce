set.seed(1)
library(MASS)
library(tidyverse)
source('~/dmf/R/dmf.R')
source('/projectnb/dmfgrp/dmf_revision/power_util.R')

std_eta = 0.4
n = 600
p = 420
q_star = 5
mean_eta = 0
phi = 1
true_family = Gamma('log')
true_weight = 1
G = n*p/50
V_star <- matrix(rnorm(q_star * p, 0, std_eta), nrow = p, ncol = q_star)
L_star <- matrix(rnorm(n * q_star, 0, std_eta), nrow = n, ncol = q_star)
eta_star <- tcrossprod(L_star, V_star) + mean_eta

Y = generate_Y(true_family, eta_star, phi = phi, glm_weights = true_weight)

#lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight, init_svd = TRUE, control = glm.control(maxit = 200))
lv_true = dmf(Y/true_weight, true_family, rank = q_star, weights = true_weight, control = glm.control(maxit = 200))
true_result = family_test(Y/true_weight, lv_true, G, dispersion = phi, weights = true_weight)
p_true = 1- pchisq(true_result$chisq_stat, length(true_result$norm_vec)-1)
print(p_true)

