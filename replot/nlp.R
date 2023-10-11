setwd('/Users/caesa/Desktop/GLM_Matrix')
if(!exists("foo", mode="function")) source("GLM_Matrix_Oct.R")
if(!exists("foo", mode="function")) source("Rank_Boost.R")
if(!exists("foo", mode="function")) source("JOR_Stat.R")
library(reshape2)
library(NMF)
library(dmf)
Y = as.matrix(read.csv('news20.csv'))
#Y = as.matrix(read.csv('news20tfidf.csv'))
#Y = as.matrix(read.csv('news20_aesc.csv'))

label= Y[, 1001]
Y = Y[, 1:1000] 
label = label[rowSums(Y)!=0]
Y= Y[rowSums(Y)!=0, ]


#phi = mean(Y)^2/(sd(Y)^2)
phi = mean(Y)^2/(sd(Y)^2 - mean(Y))
#phi = 0.06
#phi = 0.5 # for newsgroup


n = dim(Y)[1]
d = dim(Y)[2]

glm_family = negative.binomial(phi, 'log')
glm_family2 = poisson('log')
glm_family3 = poisson('identity')
glm_weights = matrix(rep(1,n*d), nrow=n , ncol = d)




# rank determination
#Rank_negbin = ACT_rank(Y, glm_family, glm_weights)
#Rank_poisson = ACT_rank(Y, glm_family2, glm_weights)
#Rank_nmf = ACT_rank(Y, glm_family3, glm_weights)
q_max = 200
Rank_negbin = onatski_rank(Y, glm_family, q_max, fit_full = TRUE)
Rank_poisson = onatski_rank(Y, glm_family2, q_max, fit_full = TRUE)


eigen_values1= eigen(cov(tcrossprod(Rank_negbin$L,Rank_negbin$V)))$value
eigen_values2= eigen(cov(tcrossprod(Rank_poisson$L,Rank_poisson$V)))$value
eigen_values3= eigen(cov(tcrossprod(fit_nmf$L,fit_nmf$V)))$value
plot(-diff(eigen_values1)[1:50])
plot(-diff(eigen_values2)[1:50])

plot_eigen = data.frame(negbinom = -diff(eigen_values1)[1:30], poisson= -diff(eigen_values2)[1:30],
                        nmf =  -diff(eigen_values3)[1:30]
)

#ggplot(cbind(melt(plot_eigen ),  index = rep(1:50,2))) + geom_point(aes(x= index, y=value,colour = variable  ))
ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=negbinom)) +
  xlab('q') +ylab('eigen diff')+
  ggtitle('Negbinom Eigen Gap') + geom_vline(xintercept = 5)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=poisson)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('Poisson Eigen Gap') + geom_vline(xintercept = 5)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=nmf)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('NMF Eigen Gap') + geom_vline(xintercept = 3)+
  theme(plot.title = element_text(hjust = 0.5))


# DMF fit
iter_matrix =  10
#fit_negbinom = IRLS_Matrix(Y, Rank_negbin, glm_family = glm_family, glm_weights = glm_weights, iter_matrix)
#fit_negbinom = IRLS_Matrix(Y, 2, glm_family = glm_family, glm_weights = glm_weights, iter_matrix)
fit_negbinom = dmf(Y, family = glm_family, 5)
#fit_poisson = IRLS_Matrix(Y, Rank_poisson, glm_family = glm_family2, glm_weights = glm_weights, iter_matrix)
#fit_poisson = IRLS_Matrix(Y, 2, glm_family = glm_family2, glm_weights = glm_weights, iter_matrix)
fit_poisson = dmf(Y, family = glm_family2, 5)

#fit_nmf_temp = nmf(Y, Rank_nmf, nrun=1, seed = 'nndsvd')
fit_nmf_temp = nmf(Y, 200, nrun=1, seed = 'nndsvd')
fit_nmf= c()
#fit_nmf$mu_hat = fitted(fit_nmf_temp)
fit_nmf$V = t(fit_nmf_temp@fit@H)
fit_nmf$L = fit_nmf_temp@fit@W
fit_nmf$family = glm_family3

# 
# fit_negbinom_adjusted = fit_negbinom
# fit_negbinom_adjusted$mu_hat = fit_negbinom_adjusted$mu_hat-2



# family determiation
#G = prod(dim(Y))/10000
#G = 150
G = 100
#G = 20 # for Four different groups with phi = 0.06 with rank 3



#norm_vec_negbinom = GOF_test(G, glm_family, fit_negbinom , Y, glm_weights, return_fit = TRUE, phi_estimate = TRUE)
#norm_vec_poisson = GOF_test(G, glm_family2, fit_poisson, Y, glm_weights, return_fit = TRUE, phi_estimate = TRUE)
#norm_vec_nmf = GOF_test(G, glm_family3, fit_nmf , Y, glm_weights, return_fit = TRUE)


# norm_vec_negbinom = family_test(Y, fit_negbinom, G, chisq_stat = FALSE)
# norm_vec_poisson = family_test(Y, fit_poisson, G, chisq_stat = FALSE)
# norm_vec_nmf = family_test(Y, fit_nmf, G, chisq_stat = FALSE)

norm_vec_negbinom = family_test(Y, fit_negbinom, G)
norm_vec_poisson = family_test(Y, fit_poisson, G)
norm_vec_nmf = family_test(Y, fit_nmf, G)

norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom, norm_vec_nmf))
colnames(norm_df) = c('Poisson','Neg Binom', 'NMF')


norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom))
colnames(norm_df) = c('Poisson','Neg Binom')

#norm_df = data.frame(cbind(norm_vec_negbinom, norm_vec_nmf))
#colnames(norm_df) = c('Neg Binom', 'NMF')


long_norm_df = melt(norm_df)

g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) + 
  geom_density(alpha = 0.5)+
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  stat_function(fun = dnorm, aes(colour = 'Standard Normal'), linetype = "dashed")+
  labs(y="Denstiy Value", x = "Y") +xlim(-10,10)


g1 = ggplot(long_norm_df) +
  stat_qq(aes(sample = value , color = variable)) +
  scale_color_brewer(palette="Dark2")+ theme_classic() +
  geom_abline(slope =1, intercept =0)+
  labs(y="Quantile(Test)", x = "Quantile(Norm)")


grid.arrange(g1,g2, ncol =2)

1-pchisq(crossprod(norm_vec_poisson)[1], G-1)
1-pchisq(crossprod(norm_vec_nmf)[1], G-1)
1-pchisq(crossprod(norm_vec_negbinom)[1], G-1)



pc_plot<-function(label, fit_result, plot_title){
  label_string = label
  label = as.numeric(factor(label))
  Rank_temp = dim(fit_result$L)[2]
  principle_result = data.frame(fit_result$L[,1:Rank_temp])
  
  principle_result = principle_result[order(label),]
  
  sep_list = rep(0, length(unique(label))-1)
  sep_list[1] = sum(label==1)
  for (j in 2:length(unique(label))){
    sep_list[j] = sep_list[j-1] + sum(label==j)
  }
  
  long_pc = melt(principle_result)
  long_pc$index = rep(1:dim(principle_result)[1], dim(principle_result)[2])
  return(ggplot(long_pc) + 
           ggtitle(plot_title) +
           geom_point( aes(x= index, y = value, colour = variable)) +  
           geom_vline(xintercept= sep_list ) +
           annotate("text", x = sep_list, y = 0.4, angle = -90, label = unique(label_string), 
                    vjust = 1.2)
  )
}

fit_negbinom$mu_hat =  glm_family$linkinv(crossprod(t(fit_negbinom$L), t(fit_negbinom$V)))

fit_nmf$W_fit = t(fit_nmf_temp@fit@H)
fit_nmf$Z_fit = fit_nmf_temp@fit@W
#label = as.numeric(factor(esGolub@phenoData@data[[2]]))
g1 = pc_plot(label, dmf_center(fit_negbinom), 'Negbinom')
g2 = pc_plot(label, dmf_center(fit_poisson), 'Poisson')
g3 = pc_plot(label, dmf_center(fit_nmf_format), 'NMF')

grid.arrange(g1,g2,g3, nrow =1)


str_label = recode(label, '0'='rec.autos' , '1' ='sci.electronics', '2' = 'rec.sport.baseball', '3' = 'comp.sys.mac.hardware' )
pc_result = data.frame( dmf_center(fit_negbinom)$L)
pc_result$label = as.factor(str_label)

pc_result2 = data.frame( dmf_center(fit_poisson)$L)
pc_result2$label = as.factor(str_label)

pc_result3 = data.frame( dmf_center(fit_nmf)$L)
pc_result3$label = as.factor(str_label)


g3 = ggplot(pc_result) + 
  geom_point(aes(x= X1, y =X3, colour = label)) +
  theme_classic()+ 
  ggtitle('Negbin(log)') #+ xlim(0, 10000) + ylim(-1000, 2000)

library(plotly)
#filter(pc_result, label %in% c(0,2,3))
plot_ly(pc_result, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3, alpha = 0.5)) %>%
  layout(scene = list(xaxis = list(title = 'PC1',
                                   range = c(-1000, 1000)),
                      yaxis = list(title = 'PC2',
                                   range = c(-1000, 1000)),
                      zaxis = list(title = 'PC3',
                                   range = c(-150, 180)))
         #paper_bgcolor = 'rgb(243, 243, 243)',
         #plot_bgcolor = 'rgb(243, 243, 243)'
  )%>%layout(
    #title = 'Negbin PC',
    legend = list(orientation = "h",   # show entries horizontally
                  xanchor = "center",  # use center of legend as anchor
                  x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))


plot_ly(pc_result2, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3, alpha = 0.5)) %>%
  layout(scene = list(xaxis = list(title = 'PC1',
                                   range = c(-100, 1000)),
                      yaxis = list(title = 'PC2',
                                   range = c(-500, 400)),
                      zaxis = list(title = 'PC3',
                                   range = c(-150, 180)))
         #paper_bgcolor = 'rgb(243, 243, 243)',
         #plot_bgcolor = 'rgb(243, 243, 243)'
  )%>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))


plot_ly(pc_result3, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3)) %>%layout(scene = list(xaxis = list(range = c(-0.1, 0.2)),
                                                        yaxis = list(range = c(-0.2, 0.3)),
                                                        zaxis = list(range = c(-0.3, 0)))
                                           #paper_bgcolor = 'rgb(243, 243, 243)',
                                           #plot_bgcolor = 'rgb(243, 243, 243)'
        ) %>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))



# cats = [ 'rec.autos',
#          'rec.motorcycles',
#          'rec.sport.baseball',
#          'rec.sport.hockey']
# 
cats = c( 'rec.autos',
          'sci.electronics',
          'rec.sport.baseball',
          'comp.sys.mac.hardware')




g4 = ggplot(pc_result2) + geom_point(aes(x=V1, y =V2, colour = label, shape = factor(label))) +
  theme_classic() + ggtitle('NMF(identity)')

g5 = ggplot(pc_result3) + 
  geom_point(aes(x=V1, y =V2, colour = label, shape = factor(label))) +
  theme_classic()+ 
  ggtitle('QuasiPoisson(log)') #+ xlim(-500, 500)


grid.arrange(g3,g4,g5, nrow =1)







fit_negbinom = dmf(Y, glm_family, rank= dim(Y)[2])
centered_l = dmf_center(fit_negbinom)
left_temp = tcrossprod(centered_l$V , diag(sqrt(diag(crossprod(centered_l$L)))))
corr_temp = sweep(left_temp, 1, sqrt(diag(tcrossprod(left_temp))), '/')
eigen_temp2 = eigen(tcrossprod(corr_temp), symmetric = TRUE)
plot(-diff(eigen_temp2$values), col= 'red')
sum(eigen_temp2$values>1)