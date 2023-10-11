setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)
library(tidyverse)



#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/BinomialExp/large_np'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/BinomialExp/large_p'


#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/phi5_largenp/'

#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio0.3/small_np'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio0.3/large_np'

# result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/large_np/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/large_n/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/large_p/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/smallest_np'

#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/huge_n/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/hhuge_n/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/middle_n/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/ratio1.0/small_n/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/WOTorthexp/meanshiftedlargep/'
#result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/WOTorthexp/meanshiftedsmallp/'
result_dir = '/projectnb/dmfgrp/dmf_revision/thm_results/WOTorthexp/meanshiftedmiddle/'

family_names = list.files(result_dir)
phi = as.numeric(regmatches(family_names, gregexpr("[0-9.]+", family_names))[[1]])
n_runs = 10
q_star = 6

family_list = c()
family_list[[1]] = negative.binomial(phi)
family_list[[2]] = poisson()
# family_list[[3]] = binomial()



residsL_nsample = 1000
residsV_nsample = 50
residseta_nsample = 1000


# [initalize the residual df]
residL_df <- data.frame(matrix(ncol = residsL_nsample + 5, nrow = 0))
colnames(residL_df ) <- c(paste0("residsL_sample", c(1:residsL_nsample)), 'n', 'p', 'run', 'family', 'ratio')
residV_df <- data.frame(matrix(ncol = residsV_nsample + 5, nrow = 0))
colnames(residV_df ) <- c(paste0("residsV_sample", c(1:residsV_nsample)), 'n', 'p', 'run', 'family', 'ratio')
resideta_df <- data.frame(matrix(ncol = residseta_nsample + 5, nrow = 0))
colnames(resideta_df ) <- c(paste0("residseta_sample", c(1:residseta_nsample)), 'n', 'p', 'run', 'family', 'ratio')


for (glm_family in family_list){
    family_dir = paste(result_dir, glm_family$family, sep = '/')
    #family_dir = paste(result_dir, sep = '/')
    file_names = list.files(family_dir)
    n_lists = as.numeric(gsub("(.*)_(.*)", "\\1", file_names))
    p_lists =  as.numeric(gsub("(.*)_(.*)", "\\2", file_names))
  checkC2_list = rep(0, length(p_lists))
  checkVeij_list = rep(0, length(p_lists))
  for (idx in 1:length(n_lists)){
      n = n_lists[idx]
      p = p_lists[idx]
      # create model dir
      np_dir = paste(family_dir, paste(n, p, sep = '_'), sep = '/')

      resideij_df = array(0, c(n_runs, n, p))
      sdedV_df =  array(0, c(n_runs, p, q_star))
      vareij_df = array(0, c(n_runs, n, p))
      
      for (run in 1:n_runs){
        param_savedir = paste(np_dir, paste('run',run, '_param_star.RData', sep = ''), sep = '/')
        dmf_savedir = paste(np_dir, paste('run',run, '_dmf.RData', sep = ''), sep = '/')
        
        # [load model from dir]
        load(param_savedir)
        load(dmf_savedir)
        
        #ratio = log(diag(crossprod(param_star$L_star))[q_star]/p)/log(n)
        ratio = log(diag(crossprod(param_star$L_star))[q_star])/log(n)
        
        
        # [compute the residuals]
        resid_L= c(param_star$L - dmf_result$L)
        resid_V = c(param_star$V - dmf_result$V)
        eta_star = tcrossprod(param_star$L, param_star$V)
        eta_hat = tcrossprod(dmf_result$L, dmf_result$V)
        resid_eta = c(eta_star - eta_hat)
        
        # [verify the conditions]
        resideij_df[run, ,] =  param_star$Y - glm_family$linkinv(eta_star) 
        sdedV_df[run, ,] = dmf_result$V
        vareij_df[run,,] = glm_family$variance(glm_family$linkinv(eta_star) )

        addL_rowidx = nrow(residL_df)+1
        residL_df[addL_rowidx, ] <- c(sample(resid_L, residsL_nsample), n, p , run, glm_family$family, ratio)
        
        addV_rowidx = nrow(residV_df)+1
        residV_df[addV_rowidx, ] <- c(sample(resid_V, residsV_nsample), n, p , run, glm_family$family, ratio)
        
        addV_rowidx = nrow(resideta_df)+1
        resideta_df[addV_rowidx, ] <- c(sample(resid_eta, residseta_nsample), n, p , run, glm_family$family, ratio )
        
    }# end of n_runs
      
  v_sd = apply(sdedV_df, c(2,3), sd)    
  eij_sd = apply(resideij_df, c(2,3), sd)
  var_eij = apply(vareij_df, c(2,3), mean)
  
  
  sdedV_df = sweep(sdedV_df, c(1,2), v_sd, '/')
  q_test = 1; n_test = 1
  checkC2_list[idx] = mean(rowSums(resideij_df[,n_test,] * sdedV_df[,,q_test]))
  checkVeij_list[idx] = mean(rowSums(var_eij))
  
  # plot(p_lists, checkC2_list, main = glm_family$family)
  # points(p_lists, 1/p_lists, col = 'red')
  # 
  # plot(p_lists, checkVeij_list, col ='blue')
  # points(p_lists, 1/p_lists, col = 'red')
  
  }# end of n,p loop
  
  plot_df = data.frame(cbind(p_lists, n_lists, checkC2_list, 1/p_lists))
  colnames(plot_df) = c('p', 'n', 'C2', '1/p')
  ggplot(plot_df) + geom_point(aes(x = p ,y = C2, colour = as.factor(n))) + 
    geom_point(aes(x = p_lists, y= 1/p_lists), shape=23) + ggtitle(glm_family$family)
  

  plot(p_lists, checkC2_list, main = glm_family$family)
  points(p_lists, 1/p_lists, col = 'red')

  plot(p_lists, checkVeij_list, col ='blue')
  points(p_lists, 1/p_lists, col = 'red')
  
}# end of family loop


# 
resideta_df = resideta_df %>% filter(! n %in% c(1000), ! p %in% c(90))
residL_df = residL_df %>% filter(! n %in% c(1000), ! p %in% c(90))
residV_df = residV_df %>% filter(! n %in% c(1000), ! p %in% c(90))

resideta_long <-  resideta_df %>%gather(resid_, value,  1:residseta_nsample)
resideta_long$value = as.numeric(resideta_long$value)
resideta_long$p = as.numeric(resideta_long$p);resideta_long$n = as.factor(as.numeric(resideta_long$n))
# 
residL_long <-  residL_df %>%gather(resid_, value, 1:residsL_nsample)
residL_long$value = as.numeric(residL_long$value)
residL_long$p = as.numeric(residL_long$p);residL_long$n = as.factor(as.numeric(residL_long$n))
# 
residV_long <-  residV_df %>% gather(resid_, value, 1:residsV_nsample)
residV_long$value = as.numeric(residV_long$value)
residV_long$p = as.numeric(residV_long$p);residV_long$n = as.factor(as.numeric(residV_long$n))

base_family = c(family_list[[1]]$family)
# base_p = c(p_lists[length(p_lists)])
# base_n = c(n_lists[length(n_lists)])
base_p = c(p_lists[1])
base_n = c(n_lists[1])

range_density = filter(residL_long, n %in% base_n, p %in% base_p, family %in% base_family )
l_lim <- range(density(range_density$value)$x)
range_density = filter(residV_long, n %in% base_n, p %in% base_p, family %in% base_family )
v_lim <- range(density(range_density$value)$x)
range_density = filter(resideta_long, n %in% base_n, p %in% base_p, family %in% base_family )
eta_lim <- range(density(range_density$value)$x)


l_lim = as.numeric(quantile(residL_long$value, c(0.001, 0.999)))
v_lim = as.numeric(quantile(residV_long$value, c(0.01, 0.99)))
eta_lim = as.numeric(quantile(resideta_long$value, c(0.01, 0.99)))

# l_lim = range(residL_long$value)
# v_lim = range(residV_long$value)
# eta_lim =  range(resideta_long$value)

ggplot(residL_long) + geom_density(aes(x = value, colour = n)) +
  xlim(l_lim) + facet_wrap(~ family +p, nrow = length(family_list),
                           labeller = labeller(
                             .multi_line = FALSE)) + theme(legend.position="bottom")


ggplot(residV_long) + geom_density(aes(x = value, colour = n)) +
  xlim(v_lim) + facet_wrap(~family +p, nrow = length(family_list),
                           labeller = labeller(
                             .multi_line = FALSE)) + theme(legend.position="bottom")

#png('/projectnb/dmfgrp/dmf_revision/figure/Final/Eta_neffect.png', units="in", width=9, height=4, res=300)
# png('/projectnb/dmfgrp/dmf_revision/figure/Final/Eta_peffect.png', units="in", width=9, height=4, res=300)
ggplot(resideta_long) + geom_density(aes(x = value, colour = as.factor(p))) + xlab(TeX('$e_{np}^{eta}$')) + ylab(TeX('KDE | $e_{np}^{eta}$')) +
  xlim(range(eta_lim)) + facet_wrap(~family , ncol= length(family_list),
                                    labeller = labeller(
    .multi_line = FALSE)) + theme(legend.position="bottom")
# dev.off()


# png('/projectnb/dmfgrp/dmf_revision/figure/largeN/ConsisL_n5000p70.png', units="in", width=9, height=8, res=300)
#png('/projectnb/dmfgrp/dmf_revision/figure/Final/L_neffect.png', units="in", width=9, height=8, res=300)
# png('/projectnb/dmfgrp/dmf_revision/figure/Final/L_peffect.png', units="in", width=9, height=8, res=300)
ggplot(residL_long) + geom_density(aes(x = value)) +
  xlab(TeX('$e_{np}^{Lambda}$')) + ylab(TeX('KDE | $e_{np}^{Lambda}$')) +
  xlim(l_lim) + facet_wrap(~ family + p, nrow = length(family_list),
                           labeller = labeller(
                             .multi_line = FALSE)) + theme(legend.position="bottom")
# dev.off()

# png('/projectnb/dmfgrp/dmf_revision/figure/largeN/ConsisV_n5000p70.png', units="in", width=9, height=8, res=300)
# png('/projectnb/dmfgrp/dmf_revision/figure/Final/V_neffect.png', units="in", width=9, height=8, res=300)
# png('/projectnb/dmfgrp/dmf_revision/figure/Final/V_peffect.png', units="in", width=9, height=8, res=300)
ggplot(residV_long) + geom_density(aes(x = value)) +
  xlab(TeX('$e_{np}^{V}$')) + ylab(TeX('KDE | $e_{np}^{V}$')) +
  xlim(v_lim) + facet_wrap(~family +p , nrow = length(family_list),
                           labeller = labeller(
                             .multi_line = FALSE)) + theme(legend.position="bottom")
# dev.off()

# [some investigation]
plot(c(density(as.numeric(residL_long$ratio))))








