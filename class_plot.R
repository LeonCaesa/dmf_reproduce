setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source('utils.R')
if(!exists("foo", mode="function")) source('class_utils.R')
if(!exists("foo", mode="function")) source('plot.R')
library(reshape2)
library(plotly)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(nnet) # multinom
library(rpart) # tree
library(caret) #knn3 with full prob
library(HandTill2001)


# [for leukemia]
# load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Leukemia/dmf_phi=1.9_onrank2and4.RData")
# save(Y, file = '/projectnb/dmfgrp/dmf_revision/data/leukemia.RData')
# save_dir = '/projectnb2/dmfgrp/dmf_revision/class_result/leukemia'
# y = recode(label, '1'='ALL' , '2' ='AML')
# save(y, file = '/projectnb/dmfgrp/dmf_revision/data/leukemia_label.RData')
# save_df = data.frame(cbind(center_identify(fit_nmf2)$V, y))
# save(save_df,  file = paste(save_dir, 'nmf2.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_nmf)$V, y))
# save(save_df,  file = paste(save_dir, 'nmf4.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_negbinom)$V, y))
# save(save_df,  file = paste(save_dir, 'negbinom.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_poisson)$V, y))
# save(save_df,  file = paste(save_dir, 'poisson.RData', sep = '/'))


# [for nlp]
# load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/NLP/4dnewsgroup_phi006_onrank3.RData")
# save(Y, file = '/projectnb/dmfgrp/dmf_revision/data/nlp.RData')
# save_dir = '/projectnb2/dmfgrp/dmf_revision/class_result/nlp'
# y = recode(label, '0'='rec.autos' , '1' ='sci.electronics', '2' = 'rec.sport.baseball', '3' = 'comp.sys.mac.hardware' )
# save(y, file = '/projectnb/dmfgrp/dmf_revision/data/nlp_label.RData')
# save_df = data.frame(cbind(center_identify(fit_negbinom)$L, y))
# save(save_df,  file = paste(save_dir, 'negbinom.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_poisson)$L, y))
# save(save_df,  file = paste(save_dir, 'poisson.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_nmf)$L, y))
# save(save_df,  file = paste(save_dir, 'nmf.RData', sep = '/'))


# [for cbcl]
# load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Computer Vision/onrank/CBCL_total.RData")
# save(Y, file = '/projectnb/dmfgrp/dmf_revision/data/cbcl.RData')
# save_dir = '/projectnb2/dmfgrp/dmf_revision/class_result/cbcl'
# y = label
# save(y, file = '/projectnb/dmfgrp/dmf_revision/data/cbcl_label.RData')
# save_df = data.frame(cbind(center_identify(fit_negbinom)$L, y))
# save(save_df,  file = paste(save_dir, 'negbinom.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_poisson)$L, y))
# save(save_df,  file = paste(save_dir, 'poisson.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_nmf_format)$L, y))
# save(save_df,  file = paste(save_dir, 'nmf.RData', sep = '/'))
# save_df = data.frame(cbind(center_identify(fit_binom)$L, y))
# save(save_df,  file = paste(save_dir, 'binom.RData', sep = '/'))


# [for network]
# load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Network/Onrank/polblogs_onrank3_1222.RData")
# save(Y, file = '/projectnb/dmfgrp/dmf_revision/data/network.RData')
# save_dir = '/projectnb2/dmfgrp/dmf_revision/class_result/network'
# y = str_label
# save(y, file = '/projectnb/dmfgrp/dmf_revision/data/network_label.RData')
# save_df = data.frame(cbind(center_identify(fit_binom)$V, y))
# save(save_df,  file = paste(save_dir, 'binom.RData', sep = '/'))
# save_df = data.frame(cbind(eigen_results$vectors[,(ncol(eigen_results$vectors)-2):ncol(eigen_results$vectors)]),y)
# save(save_df,  file = paste(save_dir, 'spectrum.RData', sep = '/'))




class_dir = '/projectnb2/dmfgrp/dmf_revision/class_result'
data_names = c("cbcl", "leukemia" , "network", "nlp")
n_repeats = 20
# q_hat = 20
q_hat = 20


auc_df <- data.frame( auc_tree= numeric(0), auc_multi = numeric(0), auc_knn = numeric(0),
                      dataset = character(0), model = character(0), repeats = numeric(0))



for (data_name in data_names){

    model_names = list.files(paste(class_dir, data_name, sep = '/'))

    for (model_name in model_names){
          if (endsWith(model_name, paste('q', q_hat, '.RData', sep = ''))){
            load( paste(class_dir, data_name, model_name, sep = '/'))

            for (j in 1:n_repeats){
                SplitFlag = TrainTest_Flag(save_df, seed_ = j, 0.5)
                splitdata = TrainTest_Split(save_df,splitflag = SplitFlag)
                q_hat = dim(splitdata$train)[2]-1
                splitdata$train[, 1:q_hat] = mutate_all(splitdata$train[,1: 1:q_hat], function(x) as.numeric(as.character(x)))
                splitdata$test[, 1:q_hat] = mutate_all(splitdata$test[,1: 1:q_hat], function(x) as.numeric(as.character(x)))

                response = splitdata$test$y
                tree_fit = Tree_tuned(splitdata$train)
                multi_fit = multinom(y~., splitdata$train, trace = FALSE)
                knn_fit = knn3(y~., splitdata$train, k =9)

                tree_pred = predict(tree_fit, newdata = splitdata$test)
                multi_pred = predict(multi_fit, newdata = splitdata$test, type= 'prob')
                if (is.null(dim(multi_pred)[2])){multi_pred = as.matrix(cbind( 1- multi_pred, multi_pred)); colnames(multi_pred) = colnames(tree_pred)}
                knn_pred = predict(knn_fit, newdata = splitdata$test)


                htauc_tree = auc(multcap(response = factor(response),predicted = tree_pred))
                htauc_multi = auc(multcap(response = factor(response), predicted = multi_pred))
                htauc_knn = auc(multcap(response = factor(response),predicted = knn_pred))
                ht_scores = c(htauc_tree, htauc_multi, htauc_knn)

                auc_df[nrow(auc_df)+1, ] <- c(round(ht_scores,4), data_name,  strsplit(model_name,split='.', fixed=TRUE)[[1]][1],   j )
            }# end of n_repeats
          } # end of if else statement
    }# end of factor models
    print(c(data_name, model_name))
}# end of dataset
# save(auc_df, file = '/projectnb/dmfgrp/dmf_revision/class_result/auc_q20.RData')
# save(auc_df, paste("/projectnb/dmfgrp/dmf_revision/class_result/auc_", 'q', q_hat, '.RData', sep = ''))
# load(paste("/projectnb/dmfgrp/dmf_revision/class_result/auc_", 'q', q_hat, '.RData', sep = ''))

colnames(auc_df) = c("tree" ,"multi", "knn",   "dataset" ,  "model" , "repeats"  )
auc_long <- auc_df %>%
  gather(class_model, value, tree:knn)

mapfrom_dict = c( paste("binomial_q", q_hat, sep = ''), "Negative Binomial(0", "Negative Binomial(1", "Negative Binomial(10)_q20",
   paste("nmf_q", q_hat, sep = ''), paste("poisson_q", q_hat, sep = ''), paste("spectrum_q", q_hat, sep = ''))
mapto_dict =c ('DMF(Binom)', 'DMF(Negbinom)', 'DMF(Negbinom)', 'DMF(Negbinom)', 'NMF',"DMF(Poisson)","Spectrum")

level_key = c(mapfrom_dict = mapto_dict); names(level_key) = mapfrom_dict
data_key = c('CBCL', 'Leukemia', 'Network', 'NLP'); names(data_key) = c('cbcl', 'leukemia', 'network', 'nlp')

auc_long$model = recode(auc_long$model, !!!level_key)
auc_long$dataset = recode(auc_long$dataset, !!!data_key)


png('/projectnb/dmfgrp/dmf_revision/figure/Final/auc_q20.png', units="in", width=9, height=6, res=300)
ggplot(auc_long) + geom_boxplot(aes(x = as.factor(class_model),
                                    y = as.numeric(value),
                                    colour = as.factor(model))) +
  ylab('Multi-AUC') + xlab('Dataset') + theme_light(base_size = 14) + theme(legend.position="bottom") +
  facet_wrap(~dataset, scales = "free") + guides(colour=guide_legend(title="Factorization Model"))
dev.off()


