Given sample size n, data dimension p, and latent 
dimension q, `power_bash.R` generate matrix data following 
four distributions 'Gamma', 'Binom', 'Negbinom', 'Poisson'.

The generation is detailed as:

- with matrix  $x \sim N$
- This sentence uses $\` and \`$ delimiters to show math inline:  $`\sqrt{3x-1}+(1+x)^2`$
- test
$$\left( \sum_{k=1}^n a_k b_k \right)^2 \leq \left( \sum_{k=1}^n a_k^2 \right) \left( \sum_{k=1}^n b_k^2 \right)$$

  
```math
 $x \sim N$
```

Runing script 
```
```shell
  Rscript ./power_bash.R ${n} ${p} ${q_star} ${std_eta} ${mean_eta} ${G} ${nrepeats} ${add_name}
 ```