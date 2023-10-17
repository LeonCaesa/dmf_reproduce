This folder contains the code to apply DMF on various benchmark datasets:
- `CBCL`: LEE, D. D. and CHOI, S. H. (1999). Learning the parts of objects by nonnegative
matrix factorization. Nature 401.
- `20NewsGroup`: LANG, K. (1995). Newsweeder: Learning to filter netnews. 331–339.
- `Karate club`: GIRVAN, M. and NEWMAN, M. E. J. (2002). Community structure in social and bio-
logical networks. Proceedings of the National Academy of Sciences 99 7821–7826.
- `Polblog`: GIRVAN, M. and NEWMAN, M. E. J. (2002). Community structure in social and bio-
logical networks. Proceedings of the National Academy of Sciences 99 7821–7826.
- `Leukemia`: GOLUB, T. R., SLONIM, D. K., TAMAYO, P., HUARD, C., GAASENBEEK, M.,
MERISOV, J. P., COLLER, H., LOH, M. L., DOWNING, J. R., CALIGIURI, M. A., BLOOMFIELD, C. D. and LANDER, E. S. (1999). Molecular classification of can- cer: class discovery and class prediction by gene expression monitoring. Science 286 531–537.

For the convenience of reproduction, we provided the `.RData` of those
 four benchmark datasets in the parent `data` directory.


The rank estimation and family test results are straightforward 
application of the `family_test` and `onatski_rank` function from 
`dmf`. For the 



