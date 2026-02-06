import sys
import random

seed = random.randint(1, 10**9)
rng = random.Random(seed)

#finput="01vcf_withS/n400_t10000_s100_p1/xaa"
f_name=sys.argv[1]
finput="03vcffile_LD/" + f_name

#finput="01vcf_withS/n400_"+target+"_"+ss+"_p1/"+index
out_p1="05vcffile_shuffled_genotype_LD/" + f_name
out_p2=out_p1.replace('p1','p2')
out_p3=out_p1.replace('p1','p3')


with open(finput,'r') as f, \
     open(out_p1,'w') as f_p1, \
     open(out_p2,'w') as f_p2, \
     open(out_p3,'w') as f_p3:
         
    for line in f:
        if line.startswith('#'):
            f_p1.write(line)
            f_p2.write(line)
            f_p3.write(line)
        else:
            item=line.strip().split('\t')
            fixed = item[:9]
            geno=item[9:]
            
            geno1 = geno.copy(); rng.shuffle(geno1)
            geno2 = geno.copy(); rng.shuffle(geno2)
            geno3 = geno.copy(); rng.shuffle(geno3)
            
            f_p1.write('\t'.join(fixed + geno1) + '\n')
            f_p2.write('\t'.join(fixed + geno2) + '\n')
            f_p3.write('\t'.join(fixed + geno3) + '\n')

print(f"{f_name} with seed number {seed} finished! with seed number ")
