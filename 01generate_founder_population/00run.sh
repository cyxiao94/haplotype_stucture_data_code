
#1. generate index for snps and individuals
##10 selection target
##25,50,100,400,800 individuals
##i corresponding to the indepdent runs
mkdir -p 02sample_index
mkdir -p 03vcffile_LD
mkdir -p 04vcffile_LE

for num_target in 10 100 1000 5000 10000 20000
do
  for run in {1..10}
  do 
  echo "python 00script/01sample_snp_indv.py ${num_target} ${run} "
  done
done | parallel -j 20


#2.create directory to save the temporary file
for run in {1..10}
do
  for num_sample in 400
  do
    for num_target in 10 100 1000 5000 10000 20000
    do
      for selection_strength in 10
      do
        for pop in p1 p2 p3
        do    
          if [ -d 03vcffile_LD/run${run}_n${num_sample}_${num_target}loci_s${selection_strength}_${pop} ]
          then 
          echo "directory run${run}_n${num_sample}_${num_target}loci_s${selection_strength}_${pop} exists, clean the folder"
          rm -r 03vcffile_LD/run${run}_n${num_sample}_${num_target}loci_s${selection_strength}_${pop}
          fi
          mkdir 03vcffile_LD/run${run}_n${num_sample}_${num_target}loci_s${selection_strength}_${pop}
        done        
      done
    done
  done
done

#3.tranform haplo file format into vcf format for the simulation
##transform chr_pos (SLiM v4 doesn;t support multiple chromosome)
##sample indivudals from 1st step index
##designate effect size for selection target from 1st step index
for run in {1..10}
do
  for num_target in 10 100 1000 5000 10000 20000
  do
    echo "Generate founder population for Run${run}_${num_target}loci"
    for file in `ls 01hap_file`
    do
      echo "python 00script/02hap2vcf.py ${num_target} ${file} ${run}"
    done | parallel -j 30
    cat 00script/header_n400 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p1/* > 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p1.vcf && rm -r 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p1
    cat 00script/header_n400 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p2/* > 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p2.vcf && rm -r 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p2
    cat 00script/header_n400 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p3/* > 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p3.vcf && rm -r 03vcffile_LD/run${run}_n400_${num_target}loci_s10_p3
    #cat 00script/header_n400 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p1/* > 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p1.vcf && rm -r 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p1
    #cat 00script/header_n400 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p2/* > 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p2.vcf && rm -r 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p2
    #cat 00script/header_n400 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p3/* > 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p3.vcf && rm -r 04vcffile_LE/run${run}_n400_${num_target}loci_s10_p3
  done
done

#4.generate vcffile for LE simulation
#re-position to speed up the simulation
for run in {1..10}
do
  for num_target in 10 100 1000 5000 10000 20000
  do
    for pop in p1 p2 p3
    do
    echo "python 00script/03generate_LE_vcf.py 03vcffile_LD/run${run}_n400_${num_target}loci_s10_${pop}.vcf"
    done
  done
done | parallel -j 50