setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source('utils.R')

library(reshape2)
library(readr)
library(tidyr)
library(ggplot2)
library(grid)
library(NMF)
library(dmf)


get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

data_raw <- read_csv("leukemia_big.csv", col_names = FALSE )[2:5001, 2:39]

Y = as.matrix(data_raw)
#phi = mean(Y)^2/(sd(Y)^2)
#phi = 0.0001157483
phi = 1.928
#sum(Y>40000)/prod(dim(Y))

#sort(Y, decreasing = TRUE)

n = dim(Y)[1]
d = dim(Y)[2]

glm_family = negative.binomial(phi)
glm_family2 = poisson('log')
glm_family3 = poisson('identity')
glm_weights = matrix(rep(1,n*d), nrow=n , ncol = d)

#glm_weights[Y>40000] = 0


# rank determination
Rank_negbin = ACT_rank(Y, glm_family, glm_weights)
Rank_poisson = ACT_rank(Y, glm_family2, glm_weights)
Rank_nmf = ACT_rank(Y, glm_family3, glm_weights)


q_max = d-5
Rank_negbin = dmf::onatski_rank(Y, glm_family, q_max, weights = glm_weights, fit_full = TRUE)
Rank_poisson = dmf::onatski_rank(Y, glm_family2, q_max, weights = glm_weights, fit_full = TRUE)

fit_nmf_temp = nmf(Y, d, nrun=5)
fit_nmf= c()
fit_nmf$V = t(fit_nmf_temp2@fit@H)
fit_nmf$L = fit_nmf_temp2@fit@W

#Rank_nmf = dmf::onatski_rank(Y, glm_family3, q_max, glm_weights)

#dmf(Y,glm_family3, 33)
#dmf(Y, glm_family3, 10)

eigen_values1= eigen(cov(tcrossprod(Rank_negbin$L,Rank_negbin$V)))$value
eigen_values2= eigen(cov(tcrossprod(fit_poisson$L,fit_poisson$V)))$value
eigen_values3= eigen(cov(tcrossprod(fit_nmf$L,fit_nmf$V)))$value

#plot(-diff(eigen_values))

plot_eigen = data.frame(negbinom = -diff(eigen_values1)[1:30], poisson= -diff(eigen_values2)[1:30],
                        nmf = -diff(eigen_values3)[1:30])

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=negbinom)) +
  xlab('q') +ylab('eigen diff')+
  ggtitle('Negbinom Eigen Gap') + geom_vline(xintercept = 2)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=poisson)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('Poisson Eigen Gap') + geom_vline(xintercept = 1)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=nmf)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('NMF Eigen Gap') + geom_vline(xintercept = 4) +
  theme(plot.title = element_text(hjust = 0.5))


#fit_negbinom$mu_hat = glm_family$linkinv(tcrossprod(fit_negbinom$L, fit_negbinom$V))

fit_negbinom = dmf(Y, glm_family, 2, weights = glm_weights)
fit_poisson = dmf(Y, glm_family2, Rank_poisson$q_hat, weights= glm_weights)
fit_nmf_temp2 = nmf(Y, 4, nrun=5)
fit_nmf= c()
#fit_nmf$mu_hat = tcrossprod(fit_nmf_temp@fit@W, t(fit_nmf_temp@fit@H))
fit_nmf$mu_hat = fitted(fit_nmf_temp2)
fit_nmf$V = t(fit_nmf_temp2@fit@H)
fit_nmf$L = fit_nmf_temp2@fit@W


fit_nmf_rank2 = nmf(Y, 2, nrun=5)
fit_nmf2 = c()
#fit_nmf$mu_hat = tcrossprod(fit_nmf_temp@fit@W, t(fit_nmf_temp@fit@H))
fit_nmf2$mu_hat = fitted(fit_nmf_rank2)
fit_nmf2$V = t(fit_nmf_rank2@fit@H)
fit_nmf2$L = fit_nmf_rank2@fit@W
# family determiation
G = prod(dim(Y))/500


# 
# norm_vec_negbinom = GOF_test(G, glm_family, fit_negbinom , Y, glm_weights, return_fit = TRUE, phi_estimate = TRUE)
# norm_vec_poisson = GOF_test(G, glm_family2, fit_poisson, Y, glm_weights, return_fit = TRUE)
# norm_vec_nmf = GOF_test(G, glm_family3, fit_nmf , Y, glm_weights, return_fit = TRUE)

norm_vec_negbinom = family_test(Y, fit_negbinom, G, glm_weights)
norm_vec_poisson = family_test(Y, fit_poisson, G, glm_weights)
norm_vec_nmf = GOF_test(G, glm_family3, fit_nmf , Y, glm_weights, return_fit = TRUE)
#norm_vec_nmf = family_test(Y, fit_nmf,G)
norm_vec_nmf2 = GOF_test(G, glm_family3, fit_nmf2 , Y, glm_weights, return_fit = TRUE)


norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom, norm_vec_nmf, norm_vec_nmf2))
colnames(norm_df) = c('Poisson','Neg Binom', 'NMF(rank 4)', 'NMF(rank 2)')

norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom, norm_vec_nmf))
colnames(norm_df) = c('Poisson','Neg Binom', 'NMF')

norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom))
colnames(norm_df) = c('Poisson','Neg Binom')

library(reshape2)
long_norm_df = melt(norm_df)


g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) + 
  geom_density(alpha = 0.5)+
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  stat_function(fun = dnorm,  linetype = "dashed")+
  labs(y="Denstiy Value", x = "Y") + xlim(-5,5)  + ylab('KDE Value') + xlab('Pearson Residual')

g1 = ggplot(long_norm_df) +
  stat_qq(aes(sample = value , color = variable)) +
  scale_color_brewer(palette="Dark2")+ theme_classic() +
  geom_abline(slope =1, intercept =0, linetype = "dashed", aes(colour = 'Standard Normal'))+
  labs(y="Quantile(Test)", x = "Quantile(Norm)") + 
  theme(legend.position="bottom") + guides(colour=guide_legend(title="Factorization Models")) + ylim(-50,50)

legend <- get_legend(g1)
grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                         g2 + theme(legend.position="none"),
                         nrow=1),
             legend, nrow=2,heights=c(10, 1))




1-pchisq(crossprod(norm_vec_poisson)[1], G-1)
1-pchisq(crossprod(norm_vec_nmf)[1], G-1)
1-pchisq(crossprod(norm_vec_nmf2)[1], G-1)
1-pchisq(crossprod(norm_vec_negbinom)[1], G-1)



pc_plot<-function(label, fit_result, plot_title){
  label_string = label
  label = as.numeric(factor(label))
  Rank_temp = dim(fit_result$V)[2]
  principle_result = data.frame(fit_result$V[,1:Rank_temp])
  principle_result = principle_result[order(label),]
  
  
  sep_list = rep(0, length(unique(label))-1)
  sep_list[1] = sum(label==1)
  for (j in 2:length(unique(label))){
    sep_list[j] = sep_list[j-1] + sum(label==j)
  }
  
  if(dim(fit_result$V)[2] ==1){
    principle_result = data.frame(cbind( index =1:length(principle_result),  value = principle_result, variable = label))
    return(ggplot(principle_result) + 
             ggtitle(plot_title) +
             geom_point( aes(x= index, y = value, colour = variable)) +  
             geom_vline(xintercept= sep_list ) +
             annotate("text", x = sep_list, y =  mean(principle_result, na.rm = TRUE), angle = -90, label = unique(label_string), 
                      vjust = 1.2) + 
             xlab('Cell Index') + ylab('Principle Score')+ theme(plot.margin=unit(c(1,0,1,0),"cm"), plot.title = element_text(hjust = 0.5),legend.position = "none")
           
    )
    
  }
  
  long_pc = melt(principle_result)
  long_pc$index = rep(1:dim(principle_result)[1], dim(principle_result)[2])
  return(ggplot(long_pc) + 
           ggtitle(plot_title) +
           geom_point( aes(x= index, y = value, colour = variable)) +  
           geom_vline(xintercept= sep_list ) +
           annotate("text", x = sep_list, y = mean(principle_result, na.rm = TRUE), angle = -90, label = unique(label_string), 
                    vjust = 1.2)+
           xlab('Cell Index') + ylab('Principle Score')+ theme(plot.margin=unit(c(1,0,1,0),"cm"), plot.title = element_text(hjust = 0.5), legend.position = "none")
  )
}


fit_nmf$W_fit = t(fit_nmf_temp@fit@H)
fit_nmf2@W_fit = t(fit_nmf_rank2@fit@H)
data("esGolub")
label = as.numeric(factor(esGolub@phenoData@data[[2]]))
#label = c(rep(1, 19), rep(2,8), rep(3,11))
g1 = pc_plot(label, dmf_center(fit_negbinom), 'Negbinom')
g2 = pc_plot(label, dmf_center(fit_poisson), 'Poisson') 
g3 = pc_plot(label, dmf_center(fit_nmf), 'NMF(rank4)')
g4 = pc_plot(label, dmf_center(fit_nmf2), 'NMF(rank2)')
grid.arrange(g1,g2,g3, g4, nrow =2)


plot_negbin = data.frame(dmf_center(fit_negbinom)$V)
colnames(plot_negbin) = c('PC1', 'PC2')

plot_poisson= data.frame(dmf_center(fit_poisson)$V)
colnames(plot_poisson) = c('PC1')

plot_nmf= data.frame(dmf_center(fit_nmf)$V)
colnames(plot_nmf) = c('PC1', 'PC2','PC3','PC4')

plot_nmf2= data.frame(dmf_center(fit_nmf2)$V)
colnames(plot_nmf2) = c('PC1', 'PC2')


plot_negbin$label = factor(label)
plot_poisson$label = factor(label)
plot_nmf$label = factor(label)
plot_nmf2$label = factor(label)


ggplot(plot_negbin) + geom_point(aes(x= PC1, y= PC2, colour =label)) + xlab('Cell Number')
ggplot(plot_poisson) + geom_point(aes(x= 1:38, y= PC1), color = 'red') + 
  xlab('Cell Number') + geom_vline(xintercept = c(27,38)) + 
  theme(plot.margin=unit(c(1,0,1,0),"cm"), plot.title = element_text(hjust = 0.5),legend.position = "none")+
  ggtitle('Poisson')
ggplot(plot_negbin) + geom_point(aes(x= PC1, y= PC2, colour =label))
ggplot(plot_nmf) + geom_point(aes(x= PC1, y= PC2, colour =label))
ggplot(plot_nmf2) + geom_point(aes(x= PC1, y= PC2, colour =label))

pc_plot(label, dmf_center(fit_negbinom), 'Negbinom') + theme(plot.margin=unit(c(1,0,1,0),"cm"))
g2 = pc_plot(label, dmf_center(fit_poisson), 'Poisson') +  theme(plot.margin=unit(c(1,0,1,0),"cm"))
g3 = pc_plot(label, dmf_center(fit_nmf), 'NMF(rank4)')+  theme(plot.margin=unit(c(1,0,1,0),"cm"))
g4 = pc_plot(label, dmf_center(fit_nmf2), 'NMF(rank2)')+  theme(plot.margin=unit(c(1,0,1,0),"cm"))

ggplot(plot_poisson) + geom_point(aes(x= 1:38, y= PC1), color = 'red') + xlab('Cell Number') + geom_vline(xintercept = )
# 
plot_df = data.frame(fit_negbinom$V)
plot_df$X4 = label
ggplot(plot_df) + geom_point(aes(x=X2, y= X3, colour = X4))

ggplot(plot_df) + geom_point(aes(x=X1, y= X2, colour = X4))

