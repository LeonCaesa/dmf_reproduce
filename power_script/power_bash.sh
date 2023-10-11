#!/bin/bash -l

#$-P dmfgrp
#$-l h_rt=12:00:00
#$-j y


module load R/4.0.2
n=$1
p=$2
q_star=5
std_eta=$3
mean_eta=$4
G=$5
nrepeats=30
add_name=$6



start=`date +%s`
echo "n = $n | p = $p | q = $q_star | sd = $std_eta | mu = $mean_eta | G = $G | nsimu = $nrepeats | Exp = $add_name |"
Rscript /projectnb/dmfgrp/dmf_revision/power_script/power_bash.R ${n} ${p} ${q_star} ${std_eta} ${mean_eta} ${G} ${nrepeats} ${add_name}
end=`date +%s`
runtime=$((end-start))
echo " Run time of the script is ${runtime}"

