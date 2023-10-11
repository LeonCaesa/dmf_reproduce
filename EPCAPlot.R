setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)
library(tidyverse)

# [speed comparision]
load( "/projectnb/dmfgrp/dmf_revision/comp_results/time_exp.RData")
time_agg = time_df
load( "/projectnb/dmfgrp/dmf_revision/comp_results/time_exp2.RData")
time_agg = rbind(time_df, time_agg)
time_agg = filter(time_agg, !(family  %in% c('binomial')))
load("/projectnb/dmfgrp/dmf_revision/comp_results/time_BinomExp.RData")
time_agg = rbind(time_df, time_agg)


time_agg$time = round(as.numeric(time_agg$time), 4) 
time_agg$n = as.factor(as.numeric(time_agg$n))
time_agg$p = as.factor(as.numeric(time_agg$p))
# ggplot(time_agg) + geom_boxplot(aes(x = n, y = time, colour = model)) +
#   theme_light() +
#   facet_wrap(~family + p , scales = "free",
#              labeller = labeller(.multi_line = FALSE)) + theme(legend.position="bottom")


png('/projectnb/dmfgrp/dmf_revision/figure/Reply/speedcomp.png', 
    units="in", width=18, height=10, res=300)
ggplot(filter(time_agg, !(family %in% c('gaussian')))) + 
  geom_boxplot(aes(x = p, y = time, colour = model)) +
  theme_light() +
  facet_wrap(~family + n, scales = "free", nrow = 3,
             labeller = labeller(.multi_line = FALSE)) + theme(legend.position="bottom")
dev.off()

# [likelihood comparision]
deivance_df <- data.frame( n= numeric(0), p = numeric(0), 
                       family = character(0), model = character(0),
                       deviance_= numeric(0), 
                       repeats = numeric())
n_lists =  seq(200, 1000, by = 200)
p_lists =  seq(100, 200, by = 20)
nruns = 10
family_list = c()
phi = 5
family_list[[1]] = negative.binomial(phi)
family_list[[2]] = poisson()
family_list[[3]] = gaussian()
family_list[[4]] = binomial()
result_dir = '/projectnb/dmfgrp/dmf_revision/comp_results'
for (glm_family in family_list){
  family_dir = paste(result_dir, glm_family$family, sep = '/')
  for (n in n_lists){
    for(p in p_lists){
      # create model dir
      np_dir = paste(family_dir, paste(n, p, sep = '_'), sep = '/')
      for (run in 1:nruns){
  
        dmf_savedir = paste(np_dir, paste('run',run, '_dmf.RData', sep = ''), sep = '/')
        irls_savedir = paste(np_dir, paste('run',run, '_irls.RData', sep = ''), sep = '/')
        load(dmf_savedir);load(irls_savedir);
        #irls_result$deviance[length(irls_result$deviance)]
        deivance_df[nrow(deivance_df)+1, ] = c(n,p,glm_family$family, 'dmf',dmf_result$deviance, run)
        deivance_df[nrow(deivance_df)+1, ] = c(n,p,glm_family$family, 'EPCA', irls_result$deviance[length(irls_result$deviance)], run)
      }# end of run
    }# end of p
  }#end of n
}# end of family
deivance_df$deviance_ = round(as.numeric(deivance_df$deviance_), 4) 
deivance_df$avg_dev = deivance_df$deviance/as.numeric(deivance_df$n)/as.numeric(deivance_df$p)
deivance_df$n = as.factor(as.numeric(deivance_df$n))
deivance_df$p = as.factor(as.numeric(deivance_df$p))


png('/projectnb/dmfgrp/dmf_revision/figure/Reply/devcomp.png', 
    units="in", width=18, height=10, res=300)
ggplot(filter(deivance_df, !(family %in% c('gaussian')))) + geom_boxplot(aes(x = p, y = avg_dev, colour = model)) +
  theme_light() + ylab('Avg Deviance | eta')+
  facet_wrap(~family + n, scales = "free", nrow = 3,
             labeller = labeller(.multi_line = FALSE)) + theme(legend.position="bottom")
dev.off()
