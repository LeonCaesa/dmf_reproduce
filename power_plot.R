library(tidyverse)
#load_name = '/projectnb/dmfgrp/dmf_revision/family_results/power_desired.RData' # ideal offset, ideal power
load_name = '/projectnb/dmfgrp/dmf_revision/family_results/power_result5.RData' # offset 5, power 0.5
#load_name = '/projectnb/dmfgrp/dmf_revision/family_results/power_result6.RData' # offset 0, power 1
# mean magitude of \eta and std of \eta
load(load_name)

# Generate with F_2, g_2
# H_0: X_ij sim F_1(wrong) with g_1(wrong)
# H_alpha: X_ij does not following F1 and g1
# test F1 and g1, we should reject

summary_table[,3:4] = apply(summary_table[,3:4], 2, function(x) round(as.numeric(x), 4))

long_summary = pivot_longer(summary_table, cols = c('P_true', 'P_esti'), names_to = "model", values_to = 'P_value')
long_summary = unique(long_summary)


# power: probability of correctly rejecting the null hypothesis
ggplot(long_summary) + geom_boxplot(aes(x = Estimate, y = P_value, colour = model)) + facet_wrap(~True)


power_table = summary_table %>% group_by(True, Estimate) %>% summarize(power = mean(P_esti<=0.05))

ggplot(power_table, aes(x = True, y = Estimate, fill = power)) +
  geom_tile(color = 'white', alpha = 0.5) +scale_fill_gradient(low = 'white', high = 'steelblue') +
  geom_text(aes(label=power))

#filter(long_summary, True %in% c('Poisson(log)'), Estimate %in% c('Poisson(sqrt)'))
