#reference
#slim -d "input='vcffile/n400_t100_s10_p1.vcf'" selection.slim > n400_t100_s10/p1_rep1.log
#slim -d "input='../01vcffile_add_sc/n2_p1_t1000_sc0.1.vcf'" -d "sample_size=2" selection.slim >n2_t1000/p1_rep1.log


#sleep 96h



#function to run the slim simulation
run_slim_simulation(){
  num_sample=$1
  num_target=$2
  strength=$3
  runs=$4
  popsize=$5

  MAX_REP_JOBS=15

  seed_file="slim_simulation_seed/Run${runs}_popsize${popsize}_n${num_sample}_${num_target}loci_s${strength}_seeds.txt"

  for pop in p1 p2 p3
  do
    for rep in {1..5}
    do
      seed=$(shuf -i 1-2147483647 -n 1)
      echo -e "Run${runs}\tPopsize${popsize}\tn${num_sample}\t${num_target}loci\ts${strength}\t${pop}\tRep${rep}\t${seed}" >> "${seed_file}" 
      slim -s ${seed} \
         -d "vcffile='/home/vetlinux04/Changyi/04Dsim_SA_flucVScons/06slim_2025/01simulation_10runs_chr2L/01generate_founder_population/06vcffile_shuffled_genotype_LE/run${runs}_n${num_sample}_${num_target}loci_s10_${pop}.vcf'" \
         -d "popsize=${popsize}" \
         -d "strength=${strength}" \
         -d "num_target=${num_target}" \
         selection.slim > freq_matrix/run${runs}_popsize${popsize}_n${num_sample}_${num_target}loci_s${strength}_${pop}_rep${rep}.log &
    done
    while (( $(jobs -r | wc -l) >= MAX_REP_JOBS )); do
      sleep 1
    done
  done

  wait

  python get_matrix.py --runs ${runs} --popsize ${popsize} --num_sample ${num_sample} --num_target ${num_target} --strength ${strength} &&
  rm freq_matrix/run${runs}_popsize${popsize}_n${num_sample}_${num_target}loci_s${strength}_*.log

}


#start the simulations
time_start=$(date +%s)
time=$(date "+%Y-%m-%d %H:%M:%S")
echo "Start at $time"
MAX_PARAM_JOBS=3

mkdir -p slim_simulation_seed
mkdir -p freq_matrix

#neutral simulation
num_sample=400
num_target=10
strength=0
for runs in {1..10}
do
  for popsize in 1250 300
  do
    run_slim_simulation ${num_sample} ${num_target} ${strength} ${runs} ${popsize} &
    echo "[num_target=${num_target}, n=${num_sample}, s=${strength}, run=${runs}, popsize=${popsize}] launched"
  
    while (( $(jobs -r | wc -l) >= MAX_PARAM_JOBS )); do
      sleep 3
    done
  done
done


#simulation with selection
popsize=1250
for runs in {1..10}
do
  for num_sample in 400
  do
    for num_target in 10 100 1000 5000 10000 20000
    do
      for strength in 10 50 100
      do
        #mkdir -p "n${num_sample}_${num_target}loci" #create folder to store the information
        run_slim_simulation ${num_sample} ${num_target} ${strength} ${runs} ${popsize} &
        echo "[num_target=${num_target}, n=${num_sample}, s=${strength}, run=${runs}, popsize=${popsize}] launched"
        while (( $(jobs -r | wc -l) >= MAX_PARAM_JOBS )); do
          sleep 3
        done
      done
    done
  done
done 

wait


#combine seed information
echo -e "Runs\tNum_sample\tNum_target\tStrength\tSupergroup\tReplicate\tSeed" > seed_simulation
cat slim_simulation_seed/* >> seed_simulation

time_end=$(date +%s)
echo "simulations take in total $(( time_end - time_start)) seconds"

time=$(date "+%Y-%m-%d %H:%M:%S")
echo "Finish at $time"

