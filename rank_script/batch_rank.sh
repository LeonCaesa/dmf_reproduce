#!/bin/bash -l

#$-P dmfgrp
#$-l h_rt=48:00:00
#$-j y

ratio=$SGE_TASK_ID
#ratio=3

module load R/4.0.2
#nsim, case, n, p, q_star, q_max
nsim=20
n=500
p=$((50*$ratio))


start=`date +%s`
for case in "Case1" "Case2" "Case3" "Case4"
do
  q_star=6
  q_max=$((p-5))
  Rscript /projectnb/dmfgrp/dmf_revision/rank_script/RankVaryNP.R ${nsim} ${case} ${n} ${p} ${q_star} ${q_max}
done

for case in "Case5" "Case6"
do
  q_star=15
  q_max=$((p-5))
  Rscript /projectnb/dmfgrp/dmf_revision/rank_script/RankVaryNP.R ${nsim} ${case} ${n} ${p} ${q_star} ${q_max}
done
end=`date +%s`

runtime=$((end-start))

echo ${runtime}


