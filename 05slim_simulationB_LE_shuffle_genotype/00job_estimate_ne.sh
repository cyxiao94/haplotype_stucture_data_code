mkdir -p 04estimate_Ne


#estimate for neutral population with popsize 300 
for run in {1..10}
do 
  echo "Rscript 00scripts/ne_estimation.R --popsize 300 --input_file Run${run}_popsize300_n400_10loci_s0_freq_matrix.txt.gz"
done | parallel -j 10

#estimate for neutral population with popsize 1250
##no popsize parameter is needed, by defalut it is 1250
for run in {1..10}
do 
  echo "Rscript 00scripts/ne_estimation.R --input_file Run${run}_popsize1250_n400_10loci_s0_freq_matrix.txt.gz"
done | parallel -j 10


#estimate for population under selection
for run in {1..10}
do 
  for num_target in 10 100 1000 5000 10000 20000
  do
    for strength in 10 50 100
    do
    echo "Rscript 00scripts/ne_estimation.R --input_file Run${run}_popsize1250_n400_${num_target}loci_s${strength}_freq_matrix.txt.gz"
    done
  done 
done | parallel -j 10
