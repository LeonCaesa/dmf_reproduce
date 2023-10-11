setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)
library(tidyverse)
library(reshape2)
library(latex2exp)

n = 500
p = 100

np_name = paste(n, p, sep = '_')
np_dir = paste("/projectnb/dmfgrp/dmf_revision/rank_result/", np_name , "/RankCase", sep = '')
ratio_dir = paste("/projectnb/dmfgrp/dmf_revision/rank_result/", np_name , "/RatioCase", sep = '')



q_star = c(6,6,6,6,6,6)
for (i in 1:6){
  load_name = paste(np_dir, i, ".RData", sep = '')
  ratio_name = paste(ratio_dir, i, ".RData", sep = '')
  load(load_name); load(ratio_name)
  if (exists('rank_agg')){
    rank_table$case = paste('Case', i, sep = '')
    rank_agg = rbind(rank_agg, rank_table)
    ratio_agg_case$case = paste('Case', i, sep = '')
    ratio_agg = rbind(ratio_agg, ratio_agg_case)
  }else{
    rank_agg = rank_table; ratio_agg = ratio_agg_case
    rank_agg$case = paste('Case', i, sep = '')
    ratio_agg$case = paste('Case', i, sep = '')
  }
}

# [rank plot]
long_rank = melt(rank_agg)
png(paste('/projectnb/dmfgrp/dmf_revision/figure/Reply/Rank', np_name, '.png', sep = ''), 
    units="in", width=18, height=10, res=300)
ggplot(long_rank) +geom_boxplot(aes(x=variable, y= value, colour = label))+
  xlab('DMF family') + ylab(TeX('$\\hat{q}$')) + ggtitle(paste('n=', n, ', p = ', p, sep = '')) + 
  geom_hline(yintercept = q_star,  linetype='dashed')+
  theme_light() + theme(plot.title = element_text(hjust = 0.5),
                        legend.position="bottom", 
                        axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
                        text=element_text(size=12)) +
  ylim(c(0,15))+
  annotate("text", x = 2, y = q_star, label = "True Rank", vjust = -0.5) + facet_wrap(~case)
dev.off()


# [ratio plot]
long_ratio = melt(ratio_agg)
long_ratio$ratio_ = as.numeric(long_ratio$ratio_)
ggplot(long_ratio) +geom_boxplot(aes(x=family, y= ratio_))+
  xlab('DMF family') + ylab(TeX('$\\hat{q}$')) + ggtitle(np_name) + 
  theme_light() + theme(plot.title = element_text(hjust = 0.5),
                        legend.position="bottom", 
                        axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
                        text=element_text(size=12)) +facet_wrap(~case, scales = "free")