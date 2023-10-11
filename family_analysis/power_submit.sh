for exp_idx in 6 7 8 9 10
do
  add_name="GammaInit_${exp_idx}"
  for n in 200 600
  do #do n
      for ratio in 2 3 4 5 6 7
      do # do p
      let p="n/10*ratio"
      let G="n*p/50"
         for std_eta in 0.1 0.2 0.3 0.4
         do # do std_eta
              for mean_eta in 0
              do # do mean_eta
                  qsub power_bash.sh $n $p $std_eta $mean_eta $G $add_name
              done # mean eta
          done #std_eta
      done # end of p
  done # end of n
done
