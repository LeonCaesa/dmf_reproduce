setwd('/Users/caesa/Desktop/GLM_Matrix')
if(!exists("foo", mode="function")) source("GLM_Matrix_Oct.R")
if(!exists("foo", mode="function")) source("Rank_Boost.R")
if(!exists("foo", mode="function")) source("JOR_Stat.R")
library(reshape2)
library(readr)
library(tidyr)
library(ggplot2)
library(networkdata)
library(igraph)
library(network)
library(SNAData)



#data(karate)
#Y = as.matrix(as_adjacency_matrix(karate))
# data(surfersb)
# Y = as.matrix(as_adjacency_matrix(surfersb, attr = 'weight'))

# 
# data(dolphins_2)
# Y = as.matrix(as_adjacency_matrix(dolphins_2))

# eg = graph_from_data_frame(read.table('email-Eu-core.txt'), directed = TRUE) 
# Y = as.matrix(as_adjacency_matrix(eg))

#Y = as.matrix( as_adjacency_matrix(read.graph('polbooks.gml', format = 'gml'), attr = 'label'))
#g = read.graph('polbooks.gml', format = 'gml')
#Y = as.matrix(as_adjacency_matrix(g))

g = read.graph('polblogs.gml', format = 'gml')
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g_sub = igraph::induced_subgraph(g, vert_ids)
Y_sub = as.matrix(as_adjacency_matrix(as.undirected(g_sub, mode= 'mutual')))

Y = as.matrix(as_adjacency_matrix(g))
Y[Y==2] = 1


#Y = as.matrix(read.csv('revije.csv', sep = ' ',header = FALSE))

# Bank Agent Data
#Y = as.matrix(read.csv('da6.csv', sep =',')[1:55, 2:47])


# data("messages")
# Y = as(messages, "matrix")




n = dim(Y)[1]
d = dim(Y)[2]

glm_family = binomial('logit')
glm_family2 = binomial('cloglog')
glm_family3 = binomial('probit')

#glm_weights = matrix(rep(1, n*d), nrow =n, ncol =d)
#glm_weights = matrix(rep(rowSums(Y), d), nrow =n, ncol =d)
#glm_weights = t(matrix(rep(colSums(Y), n), nrow = n, ncol =d))


# onrank without weight
q_max = 200
Rank_logit = onatski_rank(x =Y, family = glm_family, q_max)
Rank_cloglog= onatski_rank(x =Y, family = glm_family2, q_max)
Rank_probit = onatski_rank(x =Y, family = glm_family3, q_max)


# onrank with weight
Rank_logit = onatski_rank(x =Y/glm_weights, family = glm_family, q_max = d-5, glm_weights)
Rank_cloglog= onatski_rank(x =Y/glm_weights, family = glm_family2, q_max = d-5, glm_weights)
Rank_probit = onatski_rank(x =Y/glm_weights, family = glm_family3, q_max = d-5, glm_weights)


eigen_values1= eigen(cov(tcrossprod(Rank_logit$L,Rank_logit$V)))$value
eigen_values2= eigen(cov(tcrossprod(Rank_cloglog$L,Rank_cloglog$V)))$value
eigen_values3= eigen(cov(tcrossprod(Rank_probit$L,Rank_probit$V)))$value
plot(-diff(eigen_values)[1:30])



plot_eigen = data.frame(logit = -diff(eigen_values1)[1:30], cloglog= -diff(eigen_values2)[1:30],
                        probit =  -diff(eigen_values3)[1:30]
)

#ggplot(cbind(melt(plot_eigen ),  index = rep(1:50,2))) + geom_point(aes(x= index, y=value,colour = variable  ))
ggplot(plot_eigen) + geom_point(aes(x= 1:30, y= logit)) +
  xlab('q') +ylab('eigen diff')+
  ggtitle('Logit Eigen Gap') + geom_vline(xintercept = 3)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=cloglog)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('Cloglog Eigen Gap') + geom_vline(xintercept = 2)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=probit)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('Polblogs Eigen Gap') + geom_vline(xintercept = 3)+
  theme(plot.title = element_text(hjust = 0.5))




plot_eigen = data.frame(club = -diff(eigen_values)[1:30])
#plot_eigen = data.frame(club = -diff(eigen_values1)[1:30], polblogs= -diff(eigen_values2)[1:30])

#ggplot(cbind(melt(plot_eigen ),  index = rep(1:50,2))) + geom_point(aes(x= index, y=value,colour = variable  ))
ggplot(plot_eigen) + geom_point(aes(x= 1:30, y =club)) +
  xlab('q') +ylab('eigen diff')+
  ggtitle('Karate Club Eigen Gap') + geom_vline(xintercept = 2)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=polblogs)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('Polblogs Eigen Gap') + geom_vline(xintercept = 2)+
  theme(plot.title = element_text(hjust = 0.5))




# rank determination
Rank_logit = ACT_rank(Y, glm_family, glm_weights)
Rank_cloglog = ACT_rank(Y, glm_family2, glm_weights)
Rank_probit = ACT_rank(Y, glm_family3, glm_weights)



# dmf fit
fit_logit = dmf(Y, glm_family, 3)
fit_cloglog = dmf(Y, glm_family2, 3)
fit_probit = dmf(Y, glm_family3, 3)

fit_logit = dmf(Y/glm_weights, glm_family, 3, glm_weights)
fit_cloglog = dmf(Y/glm_weights, glm_family2, 3, glm_weights)
fit_probit = dmf(Y/glm_weights, glm_family3, 3, glm_weights)


# family determination
#G = prod(dim(Y))/(dim(Y)[1])
G = 16
#G = 80
#G = 10 # for bank
#G = 60
#GLM fit
# norm_vec_logit = GOF_test(G, glm_family, fit_logit, Y, glm_weights, return_fit = TRUE)
# norm_vec_cloglog = GOF_test(G, glm_family2, fit_cloglog, Y, glm_weights, return_fit = TRUE)
# norm_vec_probit = GOF_test(G, glm_family3, fit_probit, Y, glm_weights, return_fit = TRUE)

# dmf fit without weights
norm_vec_logit = family_test(Y, fit_logit, G)
norm_vec_cloglog = family_test(Y, fit_cloglog, G)
norm_vec_probit = family_test(Y, fit_probit, G)


#dmf fit with weights
norm_vec_logit = family_test(Y, fit_logit, G, glm_weights)
norm_vec_cloglog = family_test(Y, fit_cloglog, G, glm_weights)
norm_vec_probit = family_test(Y, fit_probit, G, glm_weights)




norm_df = data.frame(cbind(norm_vec_logit, norm_vec_cloglog, norm_vec_probit))
colnames(norm_df) = c('logit','cloglog', 'probit')

norm_df = data.frame(cbind(norm_vec_logit, norm_vec_probit))
colnames(norm_df) = c('logit', 'probit')

norm_df = data.frame(cbind(norm_vec_cloglog, norm_vec_probit))
colnames(norm_df) = c('cloglog', 'probit')

norm_df = data.frame(cbind(norm_vec_cloglog))
colnames(norm_df) = c('cloglog')


long_norm_df = melt(norm_df)

g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) + 
  geom_density(alpha = 0.5)+
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  stat_function(fun = dnorm, aes(colour = 'Standard Normal'), linetype = "dashed")+
  labs(y="Denstiy Value", x = "Y") +xlim(-5,5)


g1 = ggplot(long_norm_df) +
  stat_qq(aes(sample = value , color = variable)) +
  scale_color_brewer(palette="Dark2")+ theme_classic() +
  geom_abline(slope =1, intercept =0)+
  labs(y="Quantile(Test)", x = "Quantile(Norm)")


grid.arrange(g1,g2, ncol =2)



pc_plot<-function(label, fit_result, plot_title){
  label_string = label
  label = as.numeric(factor(label))
  #Rank_temp = dim(fit_result$W_fit)[2]
  #principle_result = data.frame(fit_result$W_fit[,1:Rank_temp])
  Rank_temp = dim(fit_result$V)[2]
  principle_result = data.frame(fit_result$V[,1:Rank_temp])
  
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


#pol book?
# label = factor(vertex.attributes(g)$value)

# Bank agent
# label = c(rep('Accountancy',5), rep('Advertising I', 6), rep('Advertising II',5),
#           rep('Banking/Finance I',3), rep('Banking/Finance II',4), rep('Banking/Finance III',5), rep('Bank/Finance IV', 2),
#           rep('Law I',5), rep('Law II', 8), rep('Law III', 3))


label = rep(1, dim(Y)[1])
indexs = c(1, 13, 37, 16, 35, 12, 11, 10, 8, 32, 28, 25)
indexs2 = c(43,34,36, 38,7, 9, 29, 31,30 ,42,41, 39, 40)
label[indexs] =2
label[indexs2] =3

# zachary club
label = rep(1, dim(Y)[1])
indexs = c(3, 4, 14, 2, 1, 8, 22, 20, 18, 13, 12, 7, 17, 6, 5, 11)

label[indexs] =0
g1= pc_plot(label, fit_logit, 'logit')
g2= pc_plot(label, fit_cloglog, 'cloglog')
g3= pc_plot(label, fit_probit, 'probit')
grid.arrange(g1,g2,g3, nrow =1)

plot_df = data.frame(dmf_center(fit_cloglog)$V)
plot_df$label = factor(label)
ggplot(plot_df) + geom_point(aes(x =X1, y= X2, colour = label )) + xlab('PC1') + ylab('PC2') +
  ggtitle('DMF Logit')+ theme_set(theme_grey()) +theme(plot.title = element_text(hjust = 0.5))



#+ 
# geom_text(aes(x = X1, y= X2, label = 1:34), position=position_jitter(width=0.08,height=0.08))


#dolphin 2
group1 = c('Quasi', 'Mus', 'Notch', 'Number1',
           'NM23', 'Jet', 'Knit', 'Beescratch', 'DN63',
           'Upbang', 'SN89', 'SN90', 'DN21', 'Wave', 'Feather',
           'Web', 'Zig', 'Ripplefluke', 'DN16', 'TR82', 'Gallatin')
group2 = setdiff(colnames(Y), group1) 
label = rep(0, 62)
label[is.na(match(colnames(Y), group1))]=1


principle_result = data.frame(fit_cloglog$V)
principle_result = data.frame(value = principle_result[order(label),])
principle_result$index = 1:dim(principle_result)[1]
sep_list = rep(0, length(unique(label))-1)
sep_list[1] = sum(label==1)
for (j in 2:length(unique(label))){
  sep_list[j] = sep_list[j-1] + sum(label==j)
}

ggplot(principle_result) + 
  geom_point( aes( x=index, y = value)) +  
  geom_vline(xintercept= sep_list )

# bank agent
#pc_plot(label, fit_probit, 'probit')
# pc_plot(label, fit_cloglog, 'cloglog')
# pc_plot(label, fit_logit, 'logit')

# for polbooks
plot_df = data.frame(fit_logit$V)
plot_df = data.frame(dmf_center(fit_probit)$V)
plot_df = data.frame(dmf_center(fit_probit)$V)
plot_df$label = label
ggplot(plot_df) + geom_point(aes(x =X2, y= X3, colour = label ))
#ggplot(plot_df) + geom_point(aes(x =X1, y= X2, colour = label ))
plot_ly(plot_df, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label, marker = list(size = 10, alpha = 0.5))

# for polblogs
plot_df = data.frame(dmf_center(fit_probit)$V)
colnames(plot_df) <- c('PC1', 'PC2', 'PC3')
str_label = recode(factor(vertex.attributes(g)$value), '0'='liberal' , '1' ='conservative')
plot_df$label = str_label
ggplot(plot_df) + geom_point(aes(x =PC1, y= PC2, colour = label )) +
  xlim(-0.09,0.09) +ylim(-0.15,0.15)+ theme(plot.title = element_text(hjust = 0.5))

library(plotly)
#filter(pc_result, label %in% c(0,2,3))
plot_ly(plot_df, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3, alpha = 0.5)) %>%
  layout(scene = list(xaxis = list(title = 'PC1',
                                   range = c(-0.09, 0.09)),
                      yaxis = list(title = 'PC2',
                                   range = c(-0.15, 0.15)),
                      zaxis = list(title = 'PC3',
                                   range = c(-0.09, 0.09))),
         paper_bgcolor = 'rgb(243, 243, 243)',
         plot_bgcolor = 'rgb(243, 243, 243)')%>%layout(
           title = 'Negbin PC',
           legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))


a = read.table('http://moreno.ss.uci.edu/beach.dat', 
               header=FALSE, skip=6)

#as.sociomatrix(beach$bb,attr="edgevalue")

L = diag(rowSums(Y))- Y
eigen_results = eigen(L)
eigen_values = eigen_results$values
plot(1:34, -diff(eigen_values)[1:34])
#plot(eigen_results$vectors[,2], eigen_results$vectors[,3])

plot_df = data.frame(cbind(eigen_results$vectors[,33], eigen_results$vectors[,34]))
plot_df$label = factor(label)
ggplot(plot_df) + geom_point(aes(x =X1, y= X2, colour = label)) +xlab('PC1') +
  ylab('PC2')+ ggtitle('Spectral Clustering') + theme_set(theme_grey()) +theme(plot.title = element_text(hjust = 0.5))




#+ 
#  geom_text(aes(x = X1, y= X2, label = 1:34), position=position_jitter(width=0.08,height=0.08))



g =  read.graph('polblogs copy.gml', format = 'gml')
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g_sub = igraph::induced_subgraph(g, vert_ids)
Y_sub = as.matrix(as_adjacency_matrix(as.undirected(g_sub, mode= 'collapse')))
Y = Y_sub
Y[Y==2]=1


D = diag(rowSums(Y))
L = D-Y

# LRW
# Dneg1 = solve(D)
# Lrw= tcrossprod(Dneg1, L)
# eigen_results = eigen(Lrw)

# Lsym
D2 = diag(1/sqrt(rowSums(Y)))
#D2 = solve(sqrtm(D))
L_sym = tcrossprod(crossprod(D2,L), D2)
eigen_results = eigen(L_sym)



eigen_values = eigen_results$values
plot(eigen_values)


plot_df = data.frame(PC1 = eigen_results$vectors[,1220], PC2 =eigen_results$vectors[,1221] ,PC3= eigen_results$vectors[,1222])
#plot_df = data.frame(PC1 = eigen_results$vectors[1,], PC2 =eigen_results$vectors[2,] ,PC3= eigen_results$vectors[3,])
str_label = factor(vertex.attributes(g_sub)$value)
plot_df$label = str_label

ggplot(plot_df) + geom_point(aes(x= PC1, y= PC3, colour = label))+
  xlim(-0.12, 0.12) + ylim(-0.15, 0.15)+theme_light()+theme( panel.grid.minor = element_blank())



library(plotly)
#filter(pc_result, label %in% c(0,2,3))
plot_ly(plot_df, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", opacity=1, mode="markers", color=~label,
        marker = list(size = 2, alpha = 1)) %>%
  layout(scene = list(xaxis = list(title = 'PC1',
                                   range = c(-0.09, 0.09)),
                      yaxis = list(title = 'PC2',
                                   range = c(-0.05, 0.05)),
                      zaxis = list(title = 'PC3',
                                   range = c(-0.09, 0.09)))#,
         #paper_bgcolor = 'rgb(243, 243, 243)',
         #plot_bgcolor = 'rgb(243, 243, 243)'
  )%>%layout(
    #title = 'Negbin PC',
    legend = list(orientation = "h",   # show entries horizontally
                  xanchor = "center",  # use center of legend as anchor
                  x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))


library(dmf)
fit_binom = dmf(Y, binomial('logit'),3)

plot_df = data.frame(dmf_center(fit_binom)$V)
colnames(plot_df) <- c('PC1', 'PC2', 'PC3')
str_label = recode(factor(vertex.attributes(g_sub)$value), '0'='liberal' , '1' ='conservative')
plot_df$label = str_label
ggplot(plot_df) + geom_point(aes(x =PC1, y= PC2, colour = label )) +
  xlim(-0.09,0.09) +ylim(-0.15,0.15)+ theme(plot.title = element_text(hjust = 0.5))


#filter(pc_result, label %in% c(0,2,3))
plot_ly(plot_df, x=~PC1, y=~PC2, z=~PC3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 2, alpha = 0.5)) %>%
  layout(scene = list(xaxis = list(title = 'PC1',
                                   range = c(-0.02, 0)),
                      yaxis = list(title = 'PC2',
                                   range = c(-0.02, 0.01)),
                      zaxis = list(title = 'PC3',
                                   range = c(-0.05, 0.05)))
         #paper_bgcolor = 'rgb(243, 243, 243)',
         #plot_bgcolor = 'rgb(243, 243, 243)'
  )%>%layout(
    #title = 'Negbin PC',
    legend = list(orientation = "h",   # show entries horizontally
                  xanchor = "center",  # use center of legend as anchor
                  x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
