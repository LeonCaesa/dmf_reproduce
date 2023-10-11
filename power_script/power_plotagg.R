library(tidyverse)

result_dir = '/projectnb/dmfgrp/dmf_revision/power_script/power_result'
file_names = list.files(result_dir)
for (file_idx in 1:length(file_names)){
  #  if (grepl('stdchange', file_names[file_idx], fixed = TRUE)){
  #if (grepl('stdlargeG', file_names[file_idx], fixed = TRUE)){
  #if (grepl('RefinedCvg', file_names[file_idx], fixed = TRUE)){ # first result
  #if (grepl('FixGaussian', file_names[file_idx], fixed = TRUE)){
  # if (grepl('GaussianVariance', file_names[file_idx], fixed = TRUE)){
  # if (grepl('GammaExp_', file_names[file_idx], fixed = TRUE)){
  # if (grepl('GammaFix_', file_names[file_idx], fixed = TRUE)){
  if (grepl('GammaInit_', file_names[file_idx], fixed = TRUE)){
      load(paste(result_dir, file_names[file_idx], sep = '/'))
      # summary_table$add_name = str_extract(file_names[file_idx], 'GaussianVariance_.')
      # summary_table$add_name = str_extract(file_names[file_idx], 'GammaExp_.')
      # summary_table$add_name = str_extract(file_names[file_idx], 'GammaFix_.')
      summary_table$add_name = str_extract(file_names[file_idx], 'GammaInit_.')
    }else{next}
    if (exists("agg_table")== FALSE){agg_table = summary_table}else{
      colnames(summary_table) = colnames(agg_table)
      agg_table = rbind(agg_table, summary_table)
    }
}
agg_table$pvalue = as.numeric(agg_table$pvalue)
agg_table$n = as.numeric(agg_table$n)
agg_table$p = as.numeric(agg_table$p)

#group_by(n,p,std_eta, mean_eta)



# power_table = filter(agg_table, (n !=600) | (p !=60), (n !=200) | (p !=20),
#                      std_eta<= 0.4) %>% group_by(add_name, true_family, comp_family, n,p, std_eta, mean_eta) %>% summarize(power = mean(pvalue<=0.05, na.rm = TRUE))
power_table = filter(agg_table, std_eta<= 0.4) %>% group_by(add_name, true_family, comp_family, n,p, std_eta, mean_eta) %>% summarize(power = mean(pvalue<=0.05, na.rm = TRUE))

#filter(agg_table, std_eta<= 0.4) %>% group_by(add_name, true_family, comp_family, n,p, std_eta, mean_eta, G) %>% summarize(count_ = n())

power_table = power_table %>% group_by(true_family, comp_family, n,p, std_eta, mean_eta) %>% mutate(mean_power = mean(power, na.rm = TRUE))
unique(power_table$n)

power_table$mean_power = round(power_table$mean_power,2)
ggplot(filter(power_table, n == 200, std_eta<=0.9), aes(x = true_family, y = comp_family, fill = mean_power)) +
  geom_tile(color = 'white', alpha = 0.5) +scale_fill_gradient(low = 'white', high = 'steelblue') +
  geom_text(aes(label= mean_power)) + facet_wrap(~n + p + std_eta,  labeller = label_wrap_gen(multi_line=FALSE), ncol = 4)


ggplot(filter(power_table, true_family == 'Gamma')) +
  geom_boxplot(aes(x = std_eta ,y = power, colour= comp_family)) +
  facet_wrap(~ n + p, labeller = label_wrap_gen(multi_line=FALSE))


# std = 0.4/0.5; n = 200 + 600
ggplot(filter(power_table, true_family == 'Negbinom')) +
  geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family)) +
  facet_wrap(~ std_eta , labeller = label_wrap_gen(multi_line=FALSE))


# std = 0.4; n = 200 + 600
ggplot(filter(power_table, true_family == 'Poisson')) +
  geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family)) +
  facet_wrap(~ std_eta , labeller = label_wrap_gen(multi_line=FALSE))

# std = 0.4/0.5; n = 200 + 600
ggplot(filter(power_table, true_family == 'Binom')) +
  geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family)) +
  facet_wrap(~ std_eta , labeller = label_wrap_gen(multi_line=FALSE))

# std = 0.2; n = 200 + 600
ggplot(filter(power_table, true_family == 'Gamma')) +
  geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family)) +
  facet_wrap(~ std_eta , labeller = label_wrap_gen(multi_line=FALSE))




# [publication plot]
png(filename = '/projectnb/dmfgrp/dmf_revision/power_script/power_figure/extreme_power.png', units="in",
    width=8, height=4, res=300)
text.on.legend ="n, p, sigma"
g1 = ggplot(filter(power_table, n == 200, p ==40, std_eta==0.1), aes(x = true_family, y = comp_family, fill = mean_power)) +
  geom_tile(color = 'white', alpha = 0.5) +scale_fill_gradient(low = 'white', high = 'steelblue') +
  ylab('Comparing Family') + xlab ('True Family') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position="bottom") +
  guides(fill=guide_legend(title="Power")) +
  geom_text(size = 5, aes(label= mean_power)) + facet_wrap(~n + p + std_eta,
              labeller= label_bquote(.(text.on.legend)-.(paste(n, p, std_eta, sep = ', '))))

g2 = ggplot(filter(power_table, n == 600, p ==420, std_eta==0.4), aes(x = true_family, y = comp_family, fill = mean_power)) +
  geom_tile(color = 'white', alpha = 0.5) +scale_fill_gradient(low = 'white', high = 'steelblue') +
  ylab('Comparing Family') + xlab ('True Family') + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12), legend.position="bottom") +
  guides(fill=guide_legend(title="Power")) +
  geom_text(size = 5, aes(label= mean_power)) + facet_wrap(~n + p + std_eta,
     labeller= label_bquote(.(text.on.legend)-.(paste(n, p, std_eta, sep = ', '))))

ggpubr::ggarrange(g1,g2, ncol = 2, common.legend = TRUE, legend="bottom")
dev.off()



# [size effect]
library(latex2exp)
png(filename ='/projectnb/dmfgrp/dmf_revision/power_script/power_figure/size_power.png',
    units="in", width=9, height=6, res=300)
g1 = ggplot(filter(power_table, true_family == 'Poisson', std_eta == '0.4')) +
  ggtitle(TeX ('Data generated with Poisson,  $\\sigma = 0.4$')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family"))+ xlab("np") + ylab('Power') +
geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family))

g2 = ggplot(filter(power_table, true_family == 'Negbinom', std_eta == '0.4')) +
  ggtitle(TeX ('Data generated with Negbinom($\\phi$ = 5), $\\sigma = 0.4$')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family")) + xlab("np") + ylab('Power') +
  geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family))


g3 = ggplot(filter(power_table, true_family == 'Gamma', std_eta == '0.4')) +
  ggtitle(TeX ('Data generated with Gamma($\\phi$ = 1), $\\sigma = 0.4$')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family"))+ xlab("np") + ylab('Power') +
  geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family))

g4 = ggplot(filter(power_table, true_family == 'Binom', std_eta == '0.4')) +
  ggtitle(TeX ('Data generated with Binom, $\\sigma = 0.4$')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family")) + xlab("np") + ylab('Power') +
  geom_boxplot(aes(x = as.factor(n * p) ,y = power, colour= comp_family))

ggpubr::ggarrange(g1, g2, g3, g4, ncol = 2, nrow =2, common.legend = TRUE, legend="bottom")
dev.off()


# [sigma effect]
check_n = 200; check_p =80
png(filename ='/projectnb/dmfgrp/dmf_revision/power_script/power_figure/sigma_power.png',
    units="in", width=9, height=6, res=300)
g1 = ggplot(filter(power_table, true_family == 'Poisson', n == check_n, p == check_p)) +
  ggtitle(TeX ('Data generated with Poisson $(n=200, p =80)$')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family")) + xlab(TeX("$\\sigma$")) + ylab('Power') +
  geom_boxplot(aes(x = std_eta ,y = power, colour= comp_family))


g2 = ggplot(filter(power_table, true_family == 'Negbinom', n == check_n, p ==check_p)) +
  ggtitle(TeX ('Data generated with Negbinom($\\phi$ = 5, n=200, p =80)')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family")) + xlab(TeX("$\\sigma$")) + ylab('Power') +
  geom_boxplot(aes(x = std_eta ,y = power, colour= comp_family))

g3 = ggplot(filter(power_table, true_family == 'Gamma', n == check_n, p == check_p)) +
  ggtitle(TeX ('Data generated with Gamma($\\phi$ = 1, n=200, p =80)')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family")) + xlab(TeX("$\\sigma$")) + ylab('Power') +
  geom_boxplot(aes(x = std_eta ,y = power, colour= comp_family))

g4 = ggplot(filter(power_table, true_family == 'Binom', n == check_n, p == check_p)) +
  ggtitle(TeX ('Data generated with Binom $(n=200, p =80)$')) + theme_bw( base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.text.x= element_text(hjust = 0.5, vjust = 0.5, angle = 45),
        legend.position="bottom") +
  guides(colour=guide_legend(title="Comparing Family")) + xlab(TeX("$\\sigma$")) + ylab('Power') +
  geom_boxplot(aes(x = std_eta ,y = power, colour= comp_family))
ggpubr::ggarrange(g1, g2, g3, g4, ncol = 2, nrow =2, common.legend = TRUE, legend="bottom")

dev.off()



# [chisq pvalue check]
library(tidyverse)
  # [G effect]
#load("/projectnb/dmfgrp/dmf_revision/power_script/power_result/n200p80std03G25_chisqcheck.RData")
#load("/projectnb/dmfgrp/dmf_revision/power_script/power_result/n200p80std03G100_chisqcheck.RData")

  # [Size effect]
load("/projectnb/dmfgrp/dmf_revision/power_script/power_result/n600p420std04_chisqcheck.RData")
#load("/projectnb/dmfgrp/dmf_revision/power_script/power_result/n200p40std01_chisqcheck.RData")


power_plot <- function(summary_table, check_family = c("Gamma", "Binom", "Negbinom", "Poisson")){
  summary_table$chistat = as.numeric(summary_table$chistat)
  summary_table$pvalue = 1- mapply(pchisq, summary_table$chista, df = as.numeric(summary_table$G) - 1)
  summary_table$pvalue = round( as.numeric(summary_table$pvalue), 2)
  check_table = filter(summary_table, true_family %in% c(check_family))
  return(check_table)
}

check_table = power_plot(summary_table, check_family = c('Binom'))
#check_table = filter(check_table, chistat <= 400)
qqplot(check_table$chistat, rchisq(length(check_table$chistat), as.numeric(check_table$G)) -1)
abline(0,1)


check_table$power = check_table$pvalue<=0.05
check_table %>% group_by(true_family, comp_family) %>% summarise(p_ = mean(power))

