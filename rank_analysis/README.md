Given sample size `n`, data dimension `p`, and latent 
dimension $q^*$, `rank_bash.R` generate matrix data following 
 different exponential family distributions `$f$`
('Gamma', 'Binom', 'Negbinom', 'Poisson') and different links `$g$` (sqrt, log, logit etc).

$$ X_{ij} \sim f(\mu_{ij} = g(\eta_{ij})) $$

Given data matrix X, we then conduct rank estimation according to 
maximum-eigen gap (2010 Onatski) and adjust-correlation 
threshold (2020 Fan). The estimation is conducted under different
magnitude of $\eta$ and repeated multiple times. The experiment results 
are then used to produce Fig 6.

Example run:
``` 
 Rscript ./rank_bash.R ${nsim} ${case} ${n} ${p} ${q_star} ${q_max}
```

