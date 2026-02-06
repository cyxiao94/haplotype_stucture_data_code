mkdir -p 04estimate_Ne





#estimate for population under selection
for run in {1..10}
do 
  for sample_size in 400 800 1200 2500
  do
    echo "Rscript 00scripts/ne_estimation.R --input_file Run${run}_popsize1250_n${sample_size}_100loci_s10_freq_matrix.txt.gz"
  done 
done | parallel -j 20
