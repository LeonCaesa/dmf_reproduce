setwd('/projectnb2/dmfgrp/dmf_revision/')
if(!exists("foo", mode="function")) source('utils.R')
if(!exists("foo", mode="function")) source('plot.R')
library(reshape2)
library(plotly)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(igraph)
library(orca)
library(scales)
plot_family = TRUE
plot_rank = FALSE
plot_familytotal = FALSE
plot_eigen = FALSE
plot_vis = ' '
# save_path = "/projectnb/dmfgrp/dmf_revision/figure/Replot/FamilyTest"
save_path = "/projectnb/dmfgrp/dmf_revision/figure/Replot/VisPlot"
#save_path = "/projectnb/dmfgrp/dmf_revision/figure/Replot/RankTest"
# save_name = 'GOF_agg1.png'


# [for leukemia]
# load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Leukemia/dmf_phi=1.9_onrank2and4.RData")
# #labels_name = c("Poisson", TeX("NegBinom($\\hat{phi}  = 1.9$)"), "NMF(rank 4)", "NMF(rank 2)" )
# labels_name = c("Poisson", TeX("NegBinom($phi  = 1.9$)"), "NMF(rank 4)", "NMF(rank 2)" )
# ylims = c(-50, 50)
# xlims = c(-10, 10)
# #save_name = 'Leukemia_gof.png'
# save_name = 'Leukemia_vis.png'

# [for nlp]
load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/NLP/4dnewsgroup_phi006_onrank3.RData")
labels_name = unname(TeX(c("Poisson", "NegBinom($\\hat{phi} = .006$)", "NMF")))
ylims = c(-10, 10)
xlims = c(-10, 10)
str_label = recode(label, '0'='rec.autos' , '1' ='sci.electronics', '2' = 'rec.sport.baseball', '3' = 'comp.sys.mac.hardware' )
save_name = 'nlp_gof.png'


# [for cbcl]
load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Computer Vision/onrank/CBCL_total.RData")
labels_name = c("Poisson", TeX("NegBinom($\\hat{phi} = .462$)"), "NMF", "Binomial")
ylims = c(-20, 20)
xlims = c(-10, 10)
save_name = 'cbcl_gof.png'




# [factorization plot]
if (plot_vis == 'nlp'){
    dmf_indentified = dmf_identify(tcrossprod(dmf_center(fit_negbinom)$L, dmf_center(fit_negbinom)$V), 3)
    pc_result = data.frame( dmf_indentified$L)
    dmf_indentified2 = dmf_identify(tcrossprod(dmf_center(fit_poisson)$L, dmf_center(fit_poisson)$V), 3)
    pc_result2 = data.frame(dmf_indentified2$L)
    pc_result3 = data.frame(dmf_center(fit_nmf)$L)

    str_label = recode(str_label, comp.sys.mac.hardware = 'comp', rec.sport.baseball= 'baseball',
           sci.electronics = 'electronics', rec.autos = 'autos')
    pc_result$label = as.factor(str_label)
    pc_result2$label = as.factor(str_label)
    pc_result3$label = as.factor(str_label)


    scene2 = list(camera = list(eye = list(x=2, y=0.1, z=0.1)))

    scene3 = list(camera = list(eye = list(x=2, y=0.1, z=0.1)))

    #600 x 400
    Sys.setenv("plotly_username" = "thelightofking")
    Sys.setenv("plotly_api_key" = "abpAH3gahfJwyuhR3wbu")
    plotly_IMAGE(plot_3dUpdate(pc_result, c(-1000, 1000), c(-1000, 1000), c(-150, 180), font_size = 16),
                 width = 400, height = 400,
                 format = "png", out_file = paste(save_path,  "nlp_negbin.png", sep = '/'), scale = 5)
    plotly_IMAGE(plot_3dUpdate(pc_result2, c(-100, 200),  c(-200, 200), c(-150, 180), font_size = 16) %>% layout(scene = scene2),
                 width = 400, height = 400,
                 format = "png", out_file = paste(save_path,  "nlp_poisson.png", sep = '/'), scale = 5)
    plotly_IMAGE(plot_3dUpdate(pc_result3,  c(-0.1, 0.1), c(-0.1, 0.3), c(-0.15, 0), font_size = 16)%>% layout(scene = scene3),
                 width = 400, height = 400,
                 format = "png", out_file = paste(save_path,  "nlp_nmf.png", sep = '/'), scale = 5)

    library("plot3D")
    output_nlp =rbind( cbind(pc_result, model = 'Negbinom' ),
                       cbind(pc_result2, model = 'Poisson' ),
                       cbind(pc_result3, model = 'NMF' ))
    write.csv(output_nlp, 'nlp_plot.csv')

    scatter3D(pc_result$X1, pc_result$X2, pc_result$X3,phi = 0, bty = "g", pch = 18,
              col.var = as.integer(pc_result$label),
              pch = 18, ticktype = "detailed",
              colkey = list(at = c(2, 3, 4, 5), side = 1,
                            addlines = TRUE, length = 0.5, width = 0.5,
                            labels = c("setosa", "versicolor", "virginica", '1')) )

}else if (plot_vis == 'cbcl'){

    dmf_negbin_indentified = dmf_identify(tcrossprod(dmf_center(fit_negbinom)$L, dmf_center(fit_negbinom)$V), 3)
    dmf_poisson_indentified = dmf_identify(tcrossprod(dmf_center(fit_poisson)$L, dmf_center(fit_poisson)$V), 3)
    dmf_binom_indentified = dmf_identify(tcrossprod(dmf_center(fit_binom)$L, dmf_center(fit_binom)$V), 3)

    plot_negbin_df = data.frame(dmf_negbin_indentified$L)
    plot_nmf_df = data.frame(dmf_center(fit_nmf_format)$L)
    plot_poisson_df = data.frame(dmf_poisson_indentified$L)
    plot_binom_df = data.frame(dmf_binom_indentified$L)

    plot_negbin_df$label = as.factor(label)
    plot_nmf_df$label = as.factor(label)
    plot_poisson_df$label = as.factor(label)
    plot_binom_df$label = as.factor(label)

    Sys.setenv("plotly_username" = "thelightofking")
    Sys.setenv("plotly_api_key" = "abpAH3gahfJwyuhR3wbu")
    scene2 = list(camera = list(eye = list(x=-0, y=0, z=0)))
    plotly_IMAGE(plot_3dUpdate(plot_negbin_df, c(-30, 30), c(-20, 20), c(-15, 15), opacity = 0.5, font_size = 24),
                 width = 400, height = 350,
                 format = "png", out_file = paste(save_path,  "cbcl_negbin.png", sep = '/'), scale = 2)
    plotly_IMAGE(plot_3dUpdate(plot_nmf_df, c(-6, 4),  c(-6, 4), c(-6, 4), opacity = 0.5, font_size = 24),
                 width = 400, height = 350,
                 format = "png", out_file = paste(save_path,  "cbcl_nmf.png", sep = '/'), scale = 2)
    plotly_IMAGE(plot_3dUpdate(plot_binom_df,   c(-35,35),  c(-15,15),  c(-10, 10), opacity = 0.5, font_size = 24),
                 width = 400, height = 350,
                 format = "png", out_file = paste(save_path,  "cbcl_binom.png", sep = '/'), scale = 2)
    plotly_IMAGE(plot_3dUpdate(plot_poisson_df,   c(-20,20),  c(-10,10),  c(-10, 10), opacity = 0.5, font_size = 24),
                 width = 400, height = 350,
                 format = "png", out_file = paste(save_path,  "cbcl_poisson.png", sep = '/'), scale = 2)

}else if(plot_vis == 'network'){

  # [polblog]
    # [for dmf]
    load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Network/Onrank/polblogs_onrank3_1222.RData")
    #scene2 = list(camera = list(eye = list(x=2, y=0.1, z=0.1)))

    scene2 = list(camera = list(eye = list(x=2, y=1, z=0.1)))
    plot_df = data.frame(dmf_center(fit_binom)$V)
    plot_df$label = str_label
    plotly_IMAGE(plot_3dUpdate(plot_df, c(-0.02, 0.02), c(-0.02, 0.02), c(-0.02, 0.02), font_size = 24) %>% layout(scene = scene2),
                 width = 480,  height = 500,
                 format = "png", out_file = paste(save_path,  "polblog_pc1.png", sep = '/'), scale = 5)

    output_network =cbind(plot_df, model = 'DMF_Binom' )


    # [for spectrum]
    plot_df = data.frame(eigen_results$vectors[,(ncol(eigen_results$vectors)-2):ncol(eigen_results$vectors)])
    plot_df$label = str_label

    Sys.setenv("plotly_username" = "thelightofking")
    Sys.setenv("plotly_api_key" = "abpAH3gahfJwyuhR3wbu")
    plotly_IMAGE(plot_3dUpdate(plot_df, c(-0.05, 0.05), c(-0.05, 0.05), c(-0.05, 0), font_size = 24) %>% layout(scene = scene2),
                 width = 480, height = 500,
                 format = "png", out_file = paste(save_path,  "polblog_pc2.png", sep = '/'), scale = 5)

    output_network =rbind(output_network, cbind(plot_df, model = 'Spectrum' ))
    write.csv(output_network, 'network_plot.csv')

    colnames(plot_df) = c('PC1', 'PC2', 'PC3', 'label')


    save_name = 'polblog_pc3.png'
    png(paste(save_path, save_name, sep = '/'), units="in", width=3, height=3.2, res=300)
    ggplot(plot_df) + geom_point(aes(x= PC1, y= PC3, colour = label), alpha = 0.5)+
      xlim(-0.12, 0.12) + ylim(-0.1, 0) + theme_light(base_size = 12) +
      scale_colour_manual(values = c("#1f77b4", "#4AC6B7")) +
      theme( panel.grid.minor = element_blank()) + theme(
        legend.title = element_text(size=12),  legend.text = element_text(size=12),
        plot.margin=unit(x=c(10,0,0,0),units="mm"),
        plot.title = element_text(hjust = 0.5), text=element_text(size=12), legend.position="bottom")
    dev.off()


    # [for combined]
    load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Network/Onrank/polblogs_onrank3_1222.RData")
    scene2 = list(camera = list(eye = list(x=2, y=1, z=0.1)))
    plot_df = data.frame(dmf_center(fit_binom)$V)
    plot_df$label = str_label
    fig <- plot_ly(
      type = 'scatter3d',
      mode = 'markers'
    )
    fig <- fig %>%
      add_markers(
        data = plot_df, x=~X1, y=~X2, z=~X3, color=~label, opacity =1, marker = list(size = 4, alpha = 0.5)) %>%
      layout(font=t,
             scene = list(xaxis = list(title = 'PC1',
                                       range = c(-0.02, 0.02),
                                       titlefont = list(size = 24-8)
             ),
             yaxis = list(title = 'PC2',
                          range = c(-0.02, 0.02),
                          titlefont = list(size = 24-8)),
             zaxis = list(title = 'PC3',
                          range = c(-0.02, 0.02),
                          titlefont = list(size = 24-8))),
             margin = list(
               l = 0,r = 0,b = 0,t = 0,pad = 1),
             legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = 0.5,
                           font = list(size = 24)
             ))

    plot_df = data.frame(eigen_results$vectors[,(ncol(eigen_results$vectors)-2):ncol(eigen_results$vectors)])
    plot_df$label = str_label
    #p2 <- plot_3dUpdate(plot_df, c(-0.05, 0.05), c(-0.05, 0.05), c(-0.05, 0), font_size = 24) %>% layout(scene = scene2)

    fig <- fig %>%
      add_markers(
        data = plot_df, x=~X1, y=~X2, z=~X3, color=~label, opacity =1, marker = list(size = 4, alpha = 0.5)) %>%
      layout(font=t,
             scene = list(xaxis = list(title = 'PC1',
                                       range = c(-0.05, 0.05),
                                       titlefont = list(size = 24-8)
             ),
             yaxis = list(title = 'PC2',
                          range = c(-0.05, 0.05),
                          titlefont = list(size = 24-8)),
             zaxis = list(title = 'PC3',
                          range = c(-0.05, 0),
                          titlefont = list(size = 24-8))),
             margin = list(
               l = 0,r = 0,b = 0,t = 0,pad = 1),
             legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = 0.5,
                           font = list(size = 24)
             ))


    subplot(style(p1, showlegend = F), p1, nrows =2)

  # [karate club]
    # [for dmf]
    load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Network/Onrank/Zachary_club_onrank.RData")
    save_name = 'network_club.png'
    #dmf_indentified = dmf_identify(tcrossprod(dmf_center(fit_logit)$L , dmf_center(fit_logit)$V), 2)
    label[label == 0] = 'class 0'
    label[label == 1] = 'class 1'
    dmf_indentified = dmf_center(fit_logit)
    plot_df = data.frame(dmf_indentified$V)
    plot_df$label = factor(label)
    g1 = ggplot(plot_df) + geom_point(aes(x =X1, y= X2, colour = label )) + xlab('PC1') + ylab('PC2') +
      ggtitle('DMF Logit')+ theme_light(base_size = 12) +theme(plot.title = element_text(hjust = 0.5), text=element_text(size=8), legend.position="bottom")

    # [for spectrum]
    L = diag(rowSums(Y))- Y
    eigen_results = eigen(L)
    eigen_values = eigen_results$values
    #plot(1:34, -diff(eigen_values)[1:34])

    plot_df = data.frame(eigen_results$vectors[,(ncol(eigen_results$vectors)-1):ncol(eigen_results$vectors)])
    plot_df$label = factor(label)
    g2 = ggplot(plot_df) + geom_point(aes(x =X1, y= X2, colour = label)) +xlab('PC1') +
      ylab('PC2')+ ggtitle('Spectral Clustering') + theme_light() +
      theme(plot.title = element_text(hjust = 0.5),  text=element_text(size=8),
            legend.title = element_text(size=12),  legend.text = element_text(size=12)
            ,legend.position="bottom")


    legend <- get_legend(g2)

    png(paste(save_path, save_name, sep = '/'), units="in", width=6, height=3, res=300)
    grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                             g2 + theme(legend.position="none"),
                             nrow=1),
                 legend, nrow=2,heights=c(10, 1))
    dev.off()

}else{
label[label==1] = 'ALL'
label[label==2] = 'AML'

g1 = pc_plotUpdate(label, dmf_center(fit_poisson), labels_name[1])
g2 = pc_plotUpdate(label, dmf_center(fit_negbinom), labels_name[2])
g3 = pc_plotUpdate(label, dmf_center(fit_nmf), labels_name[3])
g4 = pc_plotUpdate(label, dmf_center(fit_nmf2), labels_name[4])
legend <- get_legend(g3)
png(paste(save_path, save_name, sep = '/'), units="in", width=8, height=6, res=300)
grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                         g2 + theme(legend.position="none"),
                         g3 + theme(legend.position="none"),
                         g4 + theme(legend.position="none"),
                         nrow=2),
  legend, nrow=2, heights=c(10, 1))
dev.off()
}


# [family plot]
if (plot_family == TRUE){
    long_norm_df = melt(norm_df)

    g2 = ggplot(long_norm_df, aes(x= value, fill = variable)) +
      geom_density(alpha = 0.3, labels = labels_name)+
      scale_fill_brewer(palette="Dark2", labels = labels_name) +
      theme_light( base_size = 10) +
      stat_function(fun = dnorm,  linetype = "dashed")+
       ylab('KDE | Pearson Residual') + xlab('Pearson Residual') +
      theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position="bottom") +
    guides(colour=guide_legend(title="Factorization Models")) +
      xlim(xlims[1],xlims[2])


    g1 = ggplot(long_norm_df) +
      stat_qq(aes(sample = value , color = variable)) +
      scale_color_brewer(palette="Dark2", labels = labels_name ) +
      theme_light( base_size = 10) +
      geom_abline(slope =1, intercept =0, linetype = "dashed", aes(colour = 'Standard Normal'))+
      labs(y="Quantile(Pearson Residuals)", x = "Quantile(Standard Normal)") +
      theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position="bottom") + guides(colour=guide_legend(title="Factorization Models")) +
      ylim(ylims[1], ylims[2])


    legend <- get_legend(g1)
    plot_save =grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                             g2 + theme(legend.position="none"),
                             nrow=1),
                 legend, nrow=2,heights=c(10, 1))

    png(paste(save_path, save_name, sep = '/'), units="in", width=6, height=3, res=300)
    grid.arrange(arrangeGrob(g1 + theme(legend.position="none"),
                             g2 + theme(legend.position="none"),
                             nrow=1),
                 legend, nrow=2,heights=c(10, 1))
    dev.off()
}



if (plot_rank == TRUE){
  agg_rank = data.frame(matrix(ncol = 10))
  colnames(agg_rank) = c('Negbin(log)', 'Gaussian(identity)', 'Possion(log)', 'Poisson(sqrt)', 'Gamma(log)', 'Binomial(probit)', 'Binomial(logit)', 'label', 'Case', 'q_star')
    for(i in 1:6){
      load_name = paste("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/Rank_experiment/rank_simu", i, ".RData", sep = '')
      load(load_name)
      rank_table = data.frame(rbind(t(act_agg_act), t(act_agg_onatski)))
      rank_table$label = factor(c(rep('ACT', repeats), rep('Onatski', repeats)))
      rank_table$case = paste('Case', i)
      rank_table$q_star = paste(q_star)
      colnames(rank_table) = c('Negbin(log)', 'Gaussian(identity)', 'Possion(log)', 'Poisson(sqrt)', 'Gamma(log)', 'Binomial(probit)', 'Binomial(logit)', 'label', 'Case', 'q_star')
      agg_rank = rbind(agg_rank, rank_table)
    }
    agg_rank = agg_rank[-1,]

    long_rank = melt(agg_rank, id.vars = c('Case', 'label', 'q_star'))
    long_rank[, c('q_star', 'value')] =apply(long_rank[, c('q_star', 'value')], 2, as.integer)
    png(paste(save_path, 'rank_simu.png', sep = '/'), units="in", width=10, height=8, res=300)
    ggplot(long_rank, aes(x = variable , y = value, colour = label)) + geom_boxplot() +
      xlab('DMF family') + ylab(TeX('$\\hat{q}$')) + theme_light() +
      geom_hline(aes(yintercept = q_star, linetype='')) +
      scale_y_continuous(breaks= pretty_breaks()) +
      scale_linetype_manual(name = TeX('$q^*$'), values = c(4)) +
      # annotate("text", aes( x = 2, y = q_star), label = "True Rank", vjust = 0) +
      facet_wrap(~Case, ncol = 3, scale = 'free_y') +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position="bottom",
            axis.title.y = element_text(hjust = 0.5, vjust = 0.5, angle = 0),
            axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
            text=element_text(size=13))
    dev.off()

}


if(plot_familytotal == TRUE){



  load("/projectnb/dmfgrp/dmf_revision/replot/ReproduceData/GOF_Experiments/Agg1.RData")
  base_size = 14
  g1 = remove_geom(g1, 'GeomAbline') + geom_abline(slope =1, intercept =0, linetype = "dashed", aes(colour = 'Standard Normal'))+
    labs(y="Quantile(Pearson Residuals)", x = "Quantile(Standard Normal)")  + theme_light(base_size = base_size)+ theme(text = element_text(size = base_size))
  g2 = g2 + stat_function(fun = dnorm,  linetype = "dashed")+ theme_light(base_size = base_size) +
    ylab('KDE | Pearson Residual') + xlab('Pearson Residual') + xlim(-20,20) +
    guides(fill=guide_legend(title="Factorization Models")) +
    theme(legend.position = c(0.75, 0.75)) + theme(text = element_text(size = base_size))

  g3 = g3 + theme_light(base_size = base_size) + ylab ('Pearson Residual') + xlab('Data Index')+
    theme(text = element_text(size = base_size)) +  scale_color_manual(values = c("#D95F02","#1B9E77"))


  g4 = remove_geom(g4, 'GeomAbline') + theme_light(base_size = base_size)+ geom_abline(slope =1, intercept =0, linetype = "dashed", aes(colour = 'Standard Normal'))+
    labs(y="Quantile(Pearson Residuals)", x = "Quantile(Standard Normal)") + theme_light() + theme(text = element_text(size = 12))
  g5 = g5 + stat_function(fun = dnorm,  linetype = "dashed")+ theme_light(base_size = base_size) +
    ylab('KDE | Pearson Residual') + xlab('Pearson Residual') + xlim(-20,20)+
    guides(fill=guide_legend(title="Factorization Models")) + theme(legend.position = c(0.75, 0.75)) + theme(text = element_text(size = base_size))
  g6 = g6 + theme_light(base_size = base_size) + ylab ('Pearson Residual') + xlab('Data Index') +
    theme(text = element_text(size = base_size)) +  scale_color_manual(values = c("#D95F02","#1B9E77"))

  g7 = remove_geom(g7, 'GeomAbline') +  theme_light(base_size = base_size) + geom_abline(slope =1, intercept =0, linetype = "dashed", aes(colour = 'Standard Normal'))+
    labs(y="Quantile(Pearson Residuals)", x = "Quantile(Standard Normal)") + theme_light()+ theme(text = element_text(size = 12))
  g8 = g8 + stat_function(fun = dnorm,  linetype = "dashed")+ theme_light(base_size = base_size) +
    ylab('KDE | Pearson Residual') + xlab('Pearson Residual') + xlim(-20,20)+
    guides(fill=guide_legend(title="Factorization Models")) + theme(legend.position = c(0.75, 0.75))+ theme(text = element_text(size = 12))

  g9 = g9 + theme_light(base_size = base_size) + ylab ('Pearson Residual') + xlab('Data Index') +
    theme(text = element_text(size = base_size)) +  scale_color_manual(values = c("#D95F02","#1B9E77"))


  png(paste(save_path, save_name, sep = '/'), units="in", width=10, height=10, res=300)
  grid.arrange(arrangeGrob(g1, addSmallLegend(g2, pointSize = 1, textSize = 8, spaceLegend = 0.5),
                           g3,  ncol =3, top = grid::textGrob("Gamma(log) vs Gaussian(log)", gp=grid::gpar(fontsize=base_size+1))),
               arrangeGrob(g4, addSmallLegend(g5, pointSize = 1, textSize = 8, spaceLegend = 0.5),
                           g6, ncol =3, top = grid::textGrob("Negbinom(log) vs Poisson(log)", gp=grid::gpar(fontsize=base_size+1))),
               arrangeGrob(g7, addSmallLegend(g8, pointSize = 1, textSize = 8, spaceLegend = 0.5),
                           g9, ncol=3, top = grid::textGrob("Binom(clog) vs Binom(logit)", gp=grid::gpar(fontsize=base_size+1))),
               nrow = 3)
  dev.off()

}

if (plot_eigen){


}


#Sys.setenv(PATH=paste0("/share/pkg.7/orca/5.0.1/install/bin/orca/;", Sys.getenv("PATH")))
#Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/anaconda3/bin", sep = .Platform$path.sep))

#api_create(plot_3dUpdate(pc_result, c(-1000, 1000), c(-1000, 1000), c(-150, 180)), filename = "r-docs-midwest-boxplots")


#orca(plot_3dUpdate(pc_result, c(-1000, 1000), c(-1000, 1000), c(-150, 180)),  "output.png")

# png(paste(save_path, 'nlp_negbin.png', sep = '/'), units="in", width=6, height=3, res=300)

# export_plotly2SVG(plot_3dUpdate(pc_result, c(-1000, 1000), c(-1000, 1000), c(-150, 180)),
#                               filename = 'nlp_negbin.png',
#                               parent_path = save_path,
#                               width = 800,
#                               height = 600,
#                               remove_title = FALSE,
#                               font_family = "Arial",
#                               incl_PDF_copy = FALSE,
#                               incl_PNG_copy = TRUE,
#                               png_scaling_factor = 1.8,
#                               autocrop_png = TRUE)

