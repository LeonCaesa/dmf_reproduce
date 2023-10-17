Given sample size `n`, data dimension `p`, and latent 
dimension $q^*$, `power_bash.R` generate matrix data following 
 different exponential family distributions `f`
('Gamma', 'Binom', 'Negbinom', 'Poisson' in DMF paper).
```math
X_{ij} \sim f(\mu_{ij} = g(\eta_{ij}))
```
The magnitude and variation of the data `X` is controlled by parameters 
$\sigma_{\eta}$ and $\mu_{\eta}$ 
   
Given data matrix X, we then conduct family `f` and link `g` selection 
according to a generalized Hosmer-Lemeshow test.  


`family_residual.R` compare the pearson residual statistics using
DMF under correct and incorrect `f` and `g`. It is the main script used 
to generate Fig 1-2.

`power_bash.R` compute the power of the family/link selection 
via repeated simulation. The experiment results are then used to produce Fig 3-5.
Example run script 
```shell
  Rscript ./power_bash.R ${n} ${p} ${q_star} ${std_eta} ${mean_eta} ${G} ${nrepeats} ${exp_name}
```
 


