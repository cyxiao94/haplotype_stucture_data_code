mkdir -p 05vcffile_shuffled_genotype_LD
mkdir -p 06vcffile_shuffled_genotype_LE 

#shuffle genotype
for file in `ls 03vcffile_LD/*p1.vcf`
do
    filename=$(basename "$file")
    echo "python 00script/04shuffle_genotype.py $filename"
done | parallel -j 50


#4.generate vcffile for LE simulation
#re-position to speed up the simulation
for file in `ls 05vcffile_shuffled_genotype_LD/*.vcf`
do
  filename=$(basename "$file")
  echo "python 00script/05generate_LE_vcf_shufflegenotype.py 05vcffile_shuffled_genotype_LD/$filename"
done | parallel -j 50