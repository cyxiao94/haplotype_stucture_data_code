#this script generate index for randonly sampled selection targets and sampled individuals

# %%
import sys
import numpy as np
import random
import os
#%load_ext memory_profiler

# %%
seed = np.random.randint(0, 2**32 - 1)
random.seed(seed)
num_target=int(sys.argv[1])
#num_sample=int(sys.arv[2])
output_prefix=sys.argv[2]
#os.chdir("/home/vetlinux04/Changyi/07_simulation_num_ILs/18_rerun_simulation")

# %%
#sample selection targets
num_snps=801298  ##total number of segregating sites
snps=[]
with open("00script/snps") as f:
    for line in f:
        snps.append(line.strip())
        
sampled_target=random.sample(snps, num_target)
np.savetxt(f"02sample_index/run{output_prefix}_{num_target}loci_snp_index.txt", np.array(sampled_target), fmt="%s",header="#SNP index",comments=f"#Seed {seed}\n#index showed the (transformed) postion of the SNP\n")


# %%
#sample individuals
fo=open(f"02sample_index/run{output_prefix}_{num_target}loci_sample_index.txt","w")
fo.write(f"#Seed {seed}\n")
fo.write(f"#Run {output_prefix}\n")

for num_sample in [400]:
    sample_indv=random.sample(range(4,8004), num_sample*3)
    sampled_col_p1=list(range(0,4))+sample_indv[0:num_sample]
    fo.write(f"n{num_sample}_p1\t"+",".join(map(str,sampled_col_p1))+'\n')
    sampled_col_p2=list(range(0,4))+sample_indv[num_sample:(2*num_sample)]
    fo.write(f"n{num_sample}_p2\t"+",".join(map(str,sampled_col_p2))+'\n')
    sampled_col_p3=list(range(0,4))+sample_indv[(2*num_sample):(3*num_sample)]
    fo.write(f"n{num_sample}_p3\t"+",".join(map(str,sampled_col_p3))+'\n')
    #np.savetxt(f'02sample_index/run{output_prefix}_{num_target}loci_n{num_sample}_sample_index.txt', sampled_col, delimiter="\t",fmt="%d", header=f"Run {output_prefix}", comments=f"#Seed {seed}\n" )
    # sampled_col_p2=tuple(list(range(0,4))+sampled_col[num_sample:(2*num_sample)])
    # sampled_col_p3=tuple(list(range(0,4))+sampled_col[(2*num_sample):(3*num_sample)])
fo.close()
#sample_arr=np.column_stack((sampled_col_p1,sampled_col_p2,sampled_col_p3))


