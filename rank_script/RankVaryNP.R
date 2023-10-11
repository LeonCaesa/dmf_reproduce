setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source("utils.R")
library(MASS)


argv <- commandArgs(TRUE)
if (length(argv) > 0){
  repeats <- as.numeric( argv[1] )
  case <- as.character( argv[2] )
  n <- as.numeric( argv[3] )
  p <- as.numeric( argv[4] )
  q_star <- as.numeric( argv[5] )
  q_max  <- as.numeric( argv[6] )
}


#Rscript RankVaryNP.R 500 Case2 500 50 6 30


#
repeats = 30
# case = 'Case3'
# n <- 500
# p <- 50
# #q_star <- 6
# q_star <- 15
# #q_max = p-5
# q_max = p-5



phi = 2
family_list = c()
family_list[[1]] = negative.binomial(phi)
family_list[[2]] = gaussian()
family_list[[3]] = poisson()
family_list[[4]] = poisson('sqrt')
family_list[[5]] = Gamma('log')
family_list[[6]] = binomial('probit')
family_list[[7]] = binomial()

act_agg_act = matrix(0, nrow = length(family_list), ncol = repeats)
act_agg_onatski = matrix(0, nrow = length(family_list), ncol =repeats)
ratio_names = c('ratio_', 'family', 'case', 'repeat')
ratio_agg_case = data.frame(matrix( nrow = 0, ncol = length(ratio_names))); colnames(ratio_agg_case) = ratio_names



save_dir = '/projectnb/dmfgrp/dmf_revision/rank_result'
save_dir = paste(save_dir, paste(n,p, sep = '_'), sep = '/')
dir.create(save_dir)
save_name = paste('Rank', case, '_addcorrect.RData', sep = '')
save_ratio = paste('Ratio', case, '_addcorrect.RData', sep = '')



if (!file.exists(paste(save_dir, save_name, sep = '/'))){

  paste("Number of Simulations:",repeats)
  paste("Case Num:",case)
  paste("n:",n)
  paste("p:",p)
  paste("q_star:",q_star)
  paste("q_max:",q_max)

  start = Sys.time()
for (family_index in 1:length(family_list)){
    glm_family = family_list[[family_index]]

    print(glm_family)

    temp_act = rep(0, repeats)
    temp_onatski = rep(0,repeats)
    #temp_onatski_correct = rep(0,repeats)
    for (iter in 1:repeats){
      #print(iter)
      temp_rank = rank_simu(p, q_star, n, phi, glm_family, q_max = q_max, case = case)
      temp_act[iter] = temp_rank$q_act
      temp_onatski[iter] = temp_rank$q_onaski
      #temp_onatski_correct [iter] = temp_rank$q_onaski_correct
      ratio_agg_case[nrow(ratio_agg_case)+1, ] = c(temp_rank$ratio_, case, iter)
    }

    act_agg_act[family_index,] = temp_act
    act_agg_onatski[family_index,] = temp_onatski


    print(mean(temp_act))
    print(mean(temp_onatski))
}
end = Sys.time()
# print(end-start)

# [manipulating table]

rank_table = data.frame(rbind(t(act_agg_act), t(act_agg_onatski)))
rank_table$label = factor(c(rep('ACT', repeats), rep('Onatski', repeats)))
colnames(rank_table) = c('Negbin(log)', 'Gaussian(identity)', 'Possion(log)', 'Poisson(sqrt)', 'Gamma(log)', 'Binomial(probit)', 'Binomial(logit)', 'label')

# [saving]
save(rank_table, file = paste(save_dir, save_name, sep = '/'))
save(ratio_agg_case,  file = paste(save_dir, save_ratio, sep = '/'))

}
