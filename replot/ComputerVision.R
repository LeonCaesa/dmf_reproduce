setwd('/Users/caesa/Desktop/GLM_Matrix')
if(!exists("foo", mode="function")) source("GLM_Matrix_Oct.R")
if(!exists("foo", mode="function")) source("Rank_Boost.R")
if(!exists("foo", mode="function")) source("JOR_Stat.R")
library(reshape2)
library(readr)
library(tidyr)
library(ggplot2)
library(NMF)
library(dmf)
library(R.matlab)

# cbml mit of 2 person vs non-person
# Y <- t(read_csv("face.csv", col_names = FALSE )[2:362, 2:2430])
# Y <- readMat('CBCL_Face.mat')$fea
# label <- readMat('CBCL_Face.mat')$label


# yale face of 15 people, focusing on 1 3 5 11 7 or 1, 5, 11, 7 or 1,3, 7, 10, 11
Y<- t(readMat('Yale_32x32.mat')$fea)
label = readMat('Yale_32x32.mat')$gnd
Y= Y[,label %in% c(1, 3, 7, 10, 11)]
label = label[label %in% c(1,3, 7, 10, 11)]



# cifar 10
# Y= readMat('cifar10_test.mat')$data[1:2000, 1:(32*32)]
# label = readMat('cifar10_test.mat')$labels[1:2000]
# Y= readMat('cifar10_test.mat')$data
# label = readMat('cifar10_test.mat')$labels
# Y = Y[label %in% c(1,5), 1:(32*32)]
# label = label[label %in% c(1,5)]

#mnist 10
Y <- read_csv("mnist_train.csv",col_names = FALSE)
label = Y$X1[1:2000]
Y <- as.matrix(Y[1:2000, 2:785])

#phi = mean(Y)^2/(sd(Y)^2 - mean(Y))
#phi = 88.88
#phi = 0.0411728
#phi = 3.600904 #yale face
#phi = 2.131676 #cbcl
#phi = 3.845798 #cbcl log
#phi = 5.141406 #cifar10
#phi = 0.0001209285 #mnist

# large phi experiments
phi = 50
n = dim(Y)[1]
d = dim(Y)[2]

glm_family = negative.binomial(phi, 'log')
glm_family2 = poisson('log')
glm_family3 = poisson('identity')
glm_family4 = gaussian()
#glm_family5 = Gamma('log')
#glm_family5 = quasibinomial('logit')
glm_family5 = binomial('probit')
glm_weights = matrix(rep(1,n*d), nrow=n , ncol = d)
glm_weights = matrix(255, nrow=n , ncol = d)




# rank determination
q_max = 200
Rank_negbin = onatski_rank(x = Y, family = glm_family, q_max = 100, fit_full = TRUE)
Rank_poisson = onatski_rank(x = Y, family = glm_family2, q_max = 100,fit_full = TRUE)
#Rank_nmf = onatski_rank(x = Y, family = glm_family3, q_max = q_max)
#Rank_nmf = ACT_rank(Y, glm_family3,  matrix(rep(1,n*d), nrow=n , ncol = d))
Rank_nmf = nmf(Y+0.001, q_max, nrun=1, seed = 'nndsvd')
Rank_gaussian = onatski_rank(x = Y, family = glm_family4, q_max = q_max)
#Rank_gamma = onatski_rank(x = Y, family = glm_family5, q_max = q_max)
Rank_logit = onatski_rank(x = Y/glm_weights, family = glm_family5, q_max = q_max, weights = glm_weights)
Rank_probit = onatski_rank(x = Y/glm_weights, family = binomial('probit'), q_max = q_max, weights = glm_weights)
Rank_cloglog = onatski_rank(x = Y/glm_weights, family = binomial('cloglog'), q_max = q_max, weights = glm_weights)


#Rank_probit = onatski_rank(x = Y/glm_weights, family = quasibinomial('probit'), q_max = q_max, weights = glm_weights)


eigen_values1= eigen(cov(tcrossprod(Rank_negbin$L,Rank_negbin$V)))$value 
# CBCL rank4
# yale rank3
# cifar10 rank3/ 4
# cifar4 rank 14
# mnist rank 12
eigen_values2= eigen(cov(tcrossprod(Rank_poisson$L,Rank_poisson$V)))$value
# CBCL rank 3
# yale rank3
# cifar10 rank3/ 7/ 
# cifar4 rank 9
# mnist rank 12
eigen_values3 = eigen(cov(tcrossprod(Rank_nmf@fit@W, t(Rank_nmf@fit@H))))$value
eigen_values4 = eigen(cov(tcrossprod(Rank_logit$L,Rank_logit$V)))$value
eigen_values5 = eigen(cov(tcrossprod(Rank_probit$L,Rank_probit$V)))$value
eigen_values6 = eigen(cov(tcrossprod(Rank_cloglog$L,Rank_cloglog$V)))$value

#plot(log(-diff(eigen_values)[1:50]))
#plot(-diff(eigen_values4)[2:50])

plot_eigen = data.frame(negbinom = -diff(eigen_values1)[1:30], poisson= -diff(eigen_values2)[1:30],
                        nmf = -diff(eigen_values3)[1:30], binom = -diff(eigen_values6)[1:30])

#ggplot(cbind(melt(plot_eigen ),  index = rep(1:50,2))) + geom_point(aes(x= index, y=value,colour = variable  ))
ggplot(plot_eigen) + geom_point(aes(x= 1:30, y =negbinom)) +
  xlab('q') +ylab('eigen diff')+
  ggtitle('Negbinom Eigen Gap') + geom_vline(xintercept = 3)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y=poisson)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('Poisson Eigen Gap') + geom_vline(xintercept = 3)+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(plot_eigen) + geom_point(aes(x= 1:30, y= nmf)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('NMF Eigen Gap') + geom_vline(xintercept = 4)+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(plot_eigen) + geom_point(aes(x= 1:30, y= binom)) + xlab('q') +
  ylab('eigen diff')+
  ggtitle('Binomial Eigen Gap') + geom_vline(xintercept = 3)+
  theme(plot.title = element_text(hjust = 0.5))


# GLM fit
# iter_matrix =  50
# fit_negbinom = IRLS_Matrix(Y, Rank_negbin, glm_family = glm_family, glm_weights = glm_weights, iter_matrix)
# fit_poisson = IRLS_Matrix(Y, Rank_poisson, glm_family = glm_family2, glm_weights = glm_weights, iter_matrix)
# fit_nmf = nmf(Y, Rank_nmf, nrun=1, seed = 'nndsvd')
# fit_nmf_format = c()
# fit_nmf_format$mu_hat = fitted(fit_nmf)
# fit_gaussian = IRLS_Matrix(Y, Rank_gaussian, glm_family = glm_family4, glm_weights = glm_weights, iter_matrix)
# fit_gamma = IRLS_Matrix(Y+0.01, Rank_gamma, glm_family = glm_family5, glm_weights = glm_weights, iter_matrix)

# dmf fit
fit_negbinom = dmf(Y, glm_family, 4)
fit_poisson = dmf(Y, glm_family2, 3)
fit_nmf = nmf(Y+0.001, 4, nrun=1, seed = 'nndsvd')
fit_nmf_format = c()
fit_nmf_format$V = t(fit_nmf@fit@H)
fit_nmf_format$L = fit_nmf@fit@W
fit_nmf_format$family = glm_family3
fit_gaussian = dmf(Y, glm_family4, 2)
#fit_binom = dmf(Y/glm_weights, glm_family5, 2, weights= glm_weights)
fit_binom = dmf(Y/glm_weights, binomial('probit'), 3, weights= glm_weights)
#fit_gamma = dmf(Y, glm_family5, Rank_gamma)
# fit_gamma$L = fit_gamma$Z_fit
# fit_gamma$V = fit_gamma$W_fit
# fit_gamma$family = glm_family5


# Pearson DF
pearson_df = data.frame(cbind(
  pearson_binom = c((Y-fit_negbinom$mu_hat)^2/glm_family$variance(fit_negbinom$mu_hat)/ phi),
  pearson_poisson= c((Y-fit_poisson$mu_hat)^2/glm_family2$variance(fit_poisson$mu_hat)),
  pearson_nmf = c((Y-fit_nmf_format$mu_hat)^2/glm_family3$variance(fit_nmf_format$mu_hat)),
  pearson_gaussian = c((Y-fit_gaussian$mu_hat)^2/glm_family4$variance(fit_gaussian$mu_hat)/ mean((Y-fit_gaussian$mu_hat)^2)),
  pearson_gamma = c((Y-fit_negbinom$mu_hat)^2/glm_family5$variance(fit_gamma$mu_hat)* mean((Y-fit_gamma$mu_hat)^2))
))
colnames(pearson_df) = c('Neg Binom(log)', 'Poisson(log)', 'NMF', 'Gaussian(identity)', 'Gamma(log)')

# long_pearson = melt(pearson_df)
# long_pearson$Y = rep(c(Y), dim(pearson_df)[2])
# ggplot(long_pearson) + geom_point(aes(x= Y, y= value, colour = variable), alpha =0.1)+
#   scale_color_brewer(palette="Dark2")+ theme_classic() +
#   labs(y= 'Pearson_Reisdual', x = "Y")

ggplot(long_pearson, aes(variable, log(value))) + geom_boxplot() +
  scale_color_brewer(palette="Dark2")+ theme_classic() +
  labs(y= 'Pearson_Reisdual', x = "Differnet Family")

# family determiation
G = prod(dim(Y))/dim(Y)[1]
G = prod(dim(Y))/dim(Y)[2]
#G = 200
# norm_vec_negbinom = GOF_test(G, glm_family, fit_negbinom , Y, glm_weights, return_fit = TRUE, phi_estimate = TRUE)
# norm_vec_poisson = GOF_test(G, glm_family2, fit_poisson, Y, glm_weights, return_fit = TRUE)
# norm_vec_nmf = GOF_test(G, glm_family3, fit_nmf_format , Y, glm_weights, return_fit = TRUE)
# norm_vec_gaussian = GOF_test(G, glm_family4, fit_gaussian , Y, glm_weights, return_fit = TRUE, phi_estimate = TRUE)
# norm_vec_gamma = GOF_test(G, glm_family5, fit_gamma , Y, glm_weights, return_fit = TRUE, phi_estimate = TRUE)



norm_vec_negbinom = family_test(Y, fit_negbinom, G)
norm_vec_poisson = family_test(Y, fit_poisson, G)
norm_vec_nmf =family_test(Y, fit_nmf_format, G)
norm_vec_gaussian = family_test(Y, fit_gaussian, G)
norm_vec_gamma = family_test(Y, fit_gamma, G, estimate = TRUE)
norm_vec_binom = family_test(Y, fit_binom, G, weights= glm_weights)


norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom))
colnames(norm_df) = c('Poisson','Neg Binom')

norm_df = data.frame(cbind(norm_vec_nmf, norm_vec_negbinom))
colnames(norm_df) = c('NMF','Neg Binom')

norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom, norm_vec_gaussian))
colnames(norm_df) = c('Poisson','Neg Binom', 'Gaussian')


norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom, norm_vec_nmf, norm_vec_gaussian))
colnames(norm_df) = c('Poisson(log)','Neg Binom(log)', 'NMF', 'Gaussian(identity)')

norm_df = data.frame(cbind(norm_vec_poisson, norm_vec_negbinom, norm_vec_nmf, norm_vec_binom))
colnames(norm_df) = c('Poisson(log)','Neg Binom(log)', 'NMF', 'Binomial(probit)')


#long_norm_df = melt(norm_df[, c('Poisson', 'Neg Binom')])
long_norm_df = melt(norm_df)

g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) + 
  geom_density(alpha = 0.5)+
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  stat_function(fun = dnorm, aes(colour = 'Standard Normal'), linetype = "dashed")+
  labs(y="Denstiy Value", x = "Y") +xlim(-5, 5)


g1 = ggplot(long_norm_df) +
  stat_qq(aes(sample = value , color = variable)) +
  scale_color_brewer(palette="Dark2")+ theme_classic() +
  geom_abline(slope =1, intercept =0)+
  labs(y="Quantile(Test)", x = "Quantile(Norm)")


grid.arrange(g1, g2, ncol =2)




plot_fit <- function(mu_hat, rnum_pixel, cnum_pixel, num_pic, plot_title){
  pixels_gathered <- mu_hat %>%
    mutate(instance = row_number()) %>%
    gather(pixel, value, -instance) %>%
    tidyr::extract(pixel, "pixel", "(\\d+)", convert = TRUE) %>%
    mutate(pixel = pixel - 1,
           x = pixel %% rnum_pixel,
           y = cnum_pixel - pixel %/% rnum_pixel)
  theme_set(theme_light())
  
  pixels_gathered %>%
    filter(instance <= num_pic) %>%
    ggplot(aes(x, y)) + geom_raster(aes(fill= value))+
    facet_wrap(~ instance)+scale_fill_gradient(low="black",high="white")+
    ggtitle(plot_title)
  
}

# Restored picture plotting
num_faces = 4
# rnum_pixel = 19
# cnum_pixel = 19
rnum_pixel = 32
cnum_pixel = 32

# for yale face
fit_poisson$mu_hat = glm_family2$linkinv(t(crossprod(t(fit_poisson$L), t(fit_poisson$V))))
fit_negbinom$mu_hat = glm_family$linkinv(t(crossprod(t(fit_negbinom$L), t(fit_negbinom$V))))
fit_nmf_format$mu_hat = t(crossprod(t(fit_nmf_format$L), t(fit_nmf_format$V)))
fit_gaussian$mu_hat = glm_family4$linkinv(t(crossprod(t(fit_gaussian$L), t(fit_gaussian$V))))
fit_gamma$mu_hat = crossprod(t(fit_gamma$L), t(fit_gamma$V))
fit_binom$mu_hat = glm_family5$linkinv(t(crossprod(t(fit_binom$L), t(fit_binom$V))))

fit_poisson$mu_hat = glm_family2$linkinv(crossprod(t(fit_poisson$L), t(fit_poisson$V)))
fit_negbinom$mu_hat = glm_family$linkinv(crossprod(t(fit_negbinom$L), t(fit_negbinom$V)))
#fit_nmf_format$mu_hat = t(crossprod(t(fit_nmf_format$L), t(fit_nmf_format$V)))
#fit_gaussian$mu_hat = glm_family4$linkinv(t(crossprod(t(fit_gaussian$L), t(fit_gaussian$V))))



#GLM fit
p1 = plot_fit(as_tibble(fit_poisson$mu_hat), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[1])
p2 = plot_fit(as_tibble(fit_negbinom$mu_hat), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[2])
p3 = plot_fit(as_tibble(fit_nmf_format$mu_hat), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[3])
p4 = plot_fit(as_tibble(fit_gaussian$mu_hat), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[4])
p5 = plot_fit(as_tibble(fit_gamma$mu_hat), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[5])
#p6 = plot_fit(as_tibble(t(Y)), rnum_pixel, cnum_pixel, num_faces, 'Original')
p6 = plot_fit(as_tibble(Y), rnum_pixel, cnum_pixel, num_faces, 'Original')
grid.arrange(p1, p2, p3, 
             p4, p6,
             nrow = 2)

grid.arrange(p1, p2, p6, ncol =3)

# basis picture plotting
rnum_pixel = 19
cnum_pixel = 19
# rnum_pixel = 32
# cnum_pixel = 32
# GLM fit
# p1 = plot_fit(as_tibble(t(fit_negbinom$W_fit)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[1])
# p2 = plot_fit(as_tibble(t(fit_poisson$W_fit)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[2])
# p3 = plot_fit(as_tibble(fit_nmf@fit@H), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[3])
# p4 = plot_fit(as_tibble(t(fit_gaussian$W_fit)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[4])
# p5 = plot_fit(as_tibble(t(fit_gamma$W_fit)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[5])

# dmf fit
p1 = plot_fit(as_tibble(t(fit_negbinom$V)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[1])
p2 = plot_fit(as_tibble(t(fit_poisson$V)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[2])
p3 = plot_fit(as_tibble(fit_nmf@fit@H), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[3])
p4 = plot_fit(as_tibble(t(fit_gaussian$V)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[4])
p5 = plot_fit(as_tibble(t(fit_gamma$V)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[5])


# for yale face
p1 = plot_fit(as_tibble(t(fit_negbinom$L)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[2])
p2 = plot_fit(as_tibble(t(fit_poisson$L)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[1])
p3 = plot_fit(as_tibble(t(fit_nmf@fit@W)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[3])
p4 = plot_fit(as_tibble(t(fit_gaussian$L)), rnum_pixel, cnum_pixel, num_faces, colnames(norm_df)[4])


grid.arrange(p1, p2, p3, p4,
             ncol = 2)

grid.arrange(p1, p2, ncol =2)




# yale face with V
# plot_negbin_df = data.frame(dmf_center(fit_negbinom)$V)
# plot_nmf_df = data.frame(dmf_center(fit_nmf_format)$V)
# plot_poisson_df = data.frame(dmf_center(fit_poisson)$V)
# plot_binom_df = data.frame(dmf_center(fit_binom)$V)

#plot_df = data.frame(dmf_center(fit_poisson)$V)

# cifar with L
plot_negbin_df = data.frame(dmf_center(fit_negbinom)$L)
plot_nmf_df = data.frame(dmf_center(fit_nmf_format)$L)
plot_poisson_df = data.frame(dmf_center(fit_poisson)$L)
plot_binom_df = data.frame(dmf_center(fit_binom)$L)

#label[6977]=2
plot_negbin_df$label = as.factor(label)
plot_nmf_df$label = as.factor(label)
plot_poisson_df$label = as.factor(label)
plot_binom_df$label = as.factor(label)



g1 = ggplot(plot_negbin_df) + geom_point(aes(x=X2, y= X3, colour = label),size= 1) + ggtitle('Fit Neg Binom')




plot_nmf_df = data.frame(dmf_center(fit_nmf_format)$V)
#plot_nmf_df = data.frame(dmf_center(fit_nmf_format)$L)
#plot_nmf_df$label = as.factor(label)
#plot_df= plot_df[label == c(1,2,3,4,5),]
g2= ggplot(plot_nmf_df) + geom_point(aes(x=X1, y= X2, colour = label)) +  ggtitle('Fit NMF')

plot_poisson_df = data.frame(dmf_center(fit_poisson)$L)
plot_poisson_df$label = as.factor(label)

plot_gaussian_df = data.frame(dmf_center(fit_gaussian)$L)
plot_gaussian_df$label = as.factor(label)


g3 = ggplot(plot_poisson_df) + geom_point(aes(x=X1, y= X2, colour = label)) +  ggtitle('Fit Poisson')
g4 = ggplot(plot_gaussian_df) + geom_point(aes(x=X1, y= X2, colour = label)) +  ggtitle('Fit Gaussian')
grid.arrange(g1, g2, g3, g4, ncol =2)

grid.arrange(g1, g2, ncol =2)

# 3D plot
library(plotly)


plot_3d <- function(fit_dmf_df, plot_title){
  
  plot_ly(fit_dmf_df, x=~X1, y=~X2, z=~X3, opacity =1, type="scatter3d",
          mode="markers", color=~label, marker = list(size = 3))%>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3'),
                        paper_bgcolor = 'rgb(243, 243, 243)',
                        plot_bgcolor = 'rgb(243, 243, 243)'))%>%layout(
                          title = plot_title,
                          legend = list(orientation = "h",   # show entries horizontally
                                        xanchor = "center",  # use center of legend as anchor
                                        x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))
}

# Yale face
# plot_3d(filter(plot_negbin_df, label %in% c(1:15)), 'Negative Binomial(2.6)' )
# plot_3d(filter(plot_binom_df, label %in% c(1:15)), 'Binomial(cloglog)' )
# plot_3d(filter(plot_nmf_df, label %in% c(1:15)), 'NMF' )
# plot_3d(filter(plot_poisson_df, label %in% c(1,3, 7, 10, 11)), 'Poisson(log)' )


plot_ly(plot_negbin_df, x=~X1, y=~X2, z=~X3, opacity =0.5, type="scatter3d",
        mode="markers", color=~label, marker = list(size = 3))%>%
  layout(scene = list(xaxis = list(title = 'PC1',  range = c(-30, 0)),
                      yaxis = list(title = 'PC2', range = c(-20, 20)),
                      zaxis = list(title = 'PC3', range = c(-15,15))))%>%
  layout(legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))

plot_ly(plot_nmf_df, x=~X1, y=~X2, z=~X3, opacity =0.5, type="scatter3d",
        mode="markers", color=~label, marker = list(size = 3))%>%
  layout(scene = list(xaxis = list(title = 'PC1', range=c(-4,6)),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3'))
  )%>%layout(
    title = 'NMF',
    legend = list(orientation = "h",   # show entries horizontally
                  xanchor = "center",  # use center of legend as anchor
                  x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))

plot_ly(plot_binom_df, x=~X1, y=~X2, z=~X3, opacity =0.5, type="scatter3d",
        mode="markers", color=~label, marker = list(size = 3))%>%
  layout(scene = list(xaxis = list(title = 'PC1', range = c(-35,0)),
                      yaxis = list(title = 'PC2',range = c(-15,15)),
                      zaxis = list(title = 'PC3', range = c(-10, 10))))%>%
  layout(
    # title = 'Binomial(cloglog)',
    legend = list(orientation = "h",   # show entries horizontally
                  xanchor = "center",  # use center of legend as anchor
                  x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))




plot_ly(plot_poisson_df, x=~X1, y=~X2, z=~X3, opacity =0.5, type="scatter3d",
        mode="markers", color=~label,  marker = list(size = 3))%>%
  layout(scene = list(xaxis = list(title = 'PC1', range = c(-20,20)),
                      yaxis = list(title = 'PC2', range = c(-10,10)),
                      zaxis = list(title = 'PC3', range = c(-10,10))),
         #paper_bgcolor = 'rgb(243, 243, 243)',
         #plot_bgcolor = 'rgb(243, 243, 243)',
         #title = 'Poisson(log)' ,
         legend = list(orientation = "h",   # show entries horizontally
                       xanchor = "center",  # use center of legend as anchor
                       x = 0.5), margin = list(t = 0, l = 0, r= 0, b =0))






# plot_ly(plot_binom_df, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers",
#         color=~label, marker = list(size = 5))


# cifar with L

pc_plot<-function(label, fit_result, plot_title){
  label_string = label
  label = as.numeric(factor(label))
  Rank_temp = dim(fit_result$L)[2]
  principle_result = data.frame(fit_result$L[,1:Rank_temp])
  #principle_result = data.frame(fit_result$L[,1:2])
  
  # Rank_temp = dim(fit_result$V)[2]
  # principle_result = data.frame(fit_result$V[,1:Rank_temp])
  
  
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

pc_plot(label, dmf_center(fit_negbinom), 'negbinom')
pc_plot(label, dmf_center(fit_poisson), 'poisson')
pc_plot(label, dmf_center(fit_nmf_format), 'poisson')



#plot(dmf_center(fit_negbinom)$L[order(label),1])
plot(fit_poisson$L[order(label),1])
abline(v=1000)
abline(v=2000)
abline(v=3000)


#'''
#------------------------------------------------------------------------------------
# for binomial recoginition
#'''

n = dim(Y)[1]
d = dim(Y)[2]

glm_family = binomial('logit')
glm_family2 = quasibinomial('logit')
glm_family3 = binomial('probit')
glm_family4 = poisson('identity')




glm_weights = matrix(rep(255, n*d), nrow =n, ncol =d)

# dmf fit
fit_logit = dmf(Y/glm_weights, glm_family, 3, glm_weights)
fit_cloglog = dmf(Y/glm_weights, glm_family2, 3, glm_weights)
fit_probit = dmf(Y/glm_weights, glm_family3, 3, glm_weights)
fit_nmf = nmf(Y, 3, nrun=1, seed = 'nndsvd')
fit_nmf_format = c()
fit_nmf_format$V = t(fit_nmf@fit@H)
fit_nmf_format$L = fit_nmf@fit@W
fit_nmf_format$family = glm_family4

#dmf fit with weights
G = 100
norm_vec_logit = family_test(Y/glm_weights, fit_logit, G, glm_weights)
norm_vec_cloglog = family_test(Y/glm_weights, fit_cloglog, glm_weights)
norm_vec_probit = family_test(Y/glm_weights, fit_probit, G, glm_weights)
norm_vec_nmf =family_test(Y, fit_nmf_format, G)

norm_df = data.frame(cbind(norm_vec_logit, norm_vec_probit, norm_vec_nmf))
colnames(norm_df) = c('logit','probit', 'NMF')
long_norm_df = melt(norm_df)

g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) + 
  geom_density(alpha = 0.5)+
  scale_fill_brewer(palette="Dark2") +
  theme_classic() +
  stat_function(fun = dnorm, aes(colour = 'Standard Normal'), linetype = "dashed")+
  labs(y="Denstiy Value", x = "Y") +xlim(-5, 5)


g1 = ggplot(long_norm_df) +
  stat_qq(aes(sample = value , color = variable)) +
  scale_color_brewer(palette="Dark2")+ theme_classic() +
  geom_abline(slope =1, intercept =0)+
  labs(y="Quantile(Test)", x = "Quantile(Norm)")


grid.arrange(g1, g2, ncol =2)


# plot_binom_df = data.frame(dmf_center(fit_cloglog)$V)
# plot_nmf_df = data.frame(dmf_center(fit_nmf_format)$V)

plot_binom_df = data.frame(dmf_center(fit_cloglog)$L)
plot_nmf_df = data.frame(dmf_center(fit_nmf_format)$L)

plot_binom_df$label = as.factor(label)
plot_nmf_df$label = as.factor(label)


g2= ggplot(plot_nmf_df) + geom_point(aes(x=X1, y= X2, colour = label)) +  ggtitle('Fit NMF')


g1 = ggplot(plot_binom_df) + geom_point(aes(x=X1, y= X2, colour = label),size= 1) + ggtitle('Fit Binom')

grid.arrange(g1,g2, ncol =2)

library(plotly)

plot_ly(plot_binom_df, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3))
plot_ly(plot_nmf_df, x=~X1, y=~X2, z=~X3, type="scatter3d", mode="markers", color=~label,
        marker = list(size = 3))

plot(fit_logit$L[order(label),1])
plot(fit_logit$V[order(label),1])
