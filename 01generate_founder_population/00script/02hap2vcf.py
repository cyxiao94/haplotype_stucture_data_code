# %%
import sys
import numpy as np
import random
from datetime import datetime
import os

#num_sample=sys.argv[1]
num_target=sys.argv[1]
input_file=sys.argv[2]
runs=sys.argv[3]
#print(datetime.now().strftime("%Y-%m-%d %H:%M:%S"), f" Process {input_file}", flush=True)
#os.chdir("/home/vetlinux04/Changyi/07_simulation_num_ILs/18_rerun_simulation")


# %%
#load sample individuals
#load sampled index of individuals
def convert_numbers(x):
    return list(map(int, x.split(',')))

f_sample=f"02sample_index/run{runs}_{num_target}loci_sample_index.txt"

np_sample = np.loadtxt(f_sample, delimiter="\t", comments="#", dtype=object, converters={1:convert_numbers})

dict_sample = {row[0]: row[1] for row in np_sample}

#load snp index
f_snp = f"02sample_index/run{runs}_{num_target}loci_snp_index.txt"
np_snp = np.genfromtxt(f_snp, comments="#", dtype=str, autostrip=True)


# %%
#read haplo file

#fill in the information that required for vcf file
##print("adding SNP quality information", flush=True)
##information that are shared across samples
chr_2L=23000001
#chr_2R=21100001
#chr_3L=24500001
#chr_3R=27900001

chr_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=0, dtype=str)
pos_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=1, dtype=int)
num_snp=chr_arr.shape[0]

#generate new chr and pos
new_chr=np.array(["1"]*num_snp)
new_pos = pos_arr.copy()
#new_pos[chr_arr=="2L"] = pos_arr[chr_arr=="2L"]
#new_pos[chr_arr=="2R"] = pos_arr[chr_arr=="2R"] + chr_2L   
#new_pos[chr_arr=="3L"] = pos_arr[chr_arr=="3L"] + chr_2L + chr_2R
#new_pos[chr_arr=="3R"] = pos_arr[chr_arr=="3R"] + chr_2L + chr_2R + chr_3L

new_pos_LD_str = new_pos.astype(str)

#new_pos_LE =  np.arange(1, num_snp + 1)
#new_pos_LE_str =  new_pos_LE.astype(str)
 

id=new_pos_LD_str.copy()
ref=np.array(["A"]*num_snp)    
alt=np.array(["C"]*num_snp)
qual=np.array(["1000"]*num_snp)  
filter=np.array(["PASS"]*num_snp)

#temp_info=[]
info=np.array(["MID=1;S=0;DOM=0.5;PO=0;TO=1;MT=1;DP=1000"]*num_snp, dtype=object)
#if the new snp position is within the sampled selection target, update the info for that snp
mask = np.isin(new_pos_LD_str, np_snp)
#s10
effect_size_s10=10/int(num_target)
info[mask] = f"MID=1;S={effect_size_s10};DOM=0.5;PO=0;TO=1;MT=1;DP=1000"


format=np.array(["GT"]*num_snp)

merged_arr_LD=np.vstack((new_chr, new_pos_LD_str, id,ref,alt,qual,filter,info,format)).T
#merged_arr_LE=np.vstack((new_chr, new_pos_LE_str, id,ref,alt,qual,filter,info,format)).T

# %%

#n400_p1
##LD
vcf_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=tuple(dict_sample['n400_p1']), dtype=str)
vcf_arr=np.hstack([merged_arr_LD, vcf_arr[:,4:]])
#replace genotype into vcf format
vcf_arr_sub=np.where(vcf_arr == "AA", "0|0", np.where(vcf_arr=="CC","1|1", vcf_arr))
#save re-formated file into vcf file
np.savetxt(f"03vcffile_LD/run{runs}_n400_{num_target}loci_s10_p1/{input_file}.vcf", vcf_arr_sub,fmt="%s", delimiter="\t")

##LE
#vcf_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=tuple(dict_sample['n400_p1']), dtype=str)
#vcf_arr=np.hstack([merged_arr_LE, vcf_arr[:,4:]])
#replace genotype into vcf format
#vcf_arr_sub=np.where(vcf_arr == "AA", "0|0", np.where(vcf_arr=="CC","1|1", vcf_arr))
#save re-formated file into vcf file
#np.savetxt(f"04vcffile_LE/run{runs}_n400_{num_target}loci_s10_p1/{input_file}.vcf", vcf_arr_sub,fmt="%s", delimiter="\t")


#n400_p2
##LD
vcf_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=tuple(dict_sample['n400_p2']), dtype=str)
vcf_arr=np.hstack([merged_arr_LD, vcf_arr[:,4:]])
#replace genotype into vcf format
vcf_arr_sub=np.where(vcf_arr == "AA", "0|0", np.where(vcf_arr=="CC","1|1", vcf_arr))
#save re-formated file into vcf file
np.savetxt(f"03vcffile_LD/run{runs}_n400_{num_target}loci_s10_p2/{input_file}.vcf", vcf_arr_sub,fmt="%s", delimiter="\t")

##LE
#vcf_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=tuple(dict_sample['n400_p2']), dtype=str)
#vcf_arr=np.hstack([merged_arr_LE, vcf_arr[:,4:]])
#replace genotype into vcf format
#vcf_arr_sub=np.where(vcf_arr == "AA", "0|0", np.where(vcf_arr=="CC","1|1", vcf_arr))
#save re-formated file into vcf file
#np.savetxt(f"04vcffile_LE/run{runs}_n400_{num_target}loci_s10_p2/{input_file}.vcf", vcf_arr_sub,fmt="%s", delimiter="\t")

#n400_p3
##LD
vcf_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=tuple(dict_sample['n400_p3']), dtype=str)
vcf_arr=np.hstack([merged_arr_LD, vcf_arr[:,4:]])
#replace genotype into vcf format
vcf_arr_sub=np.where(vcf_arr == "AA", "0|0", np.where(vcf_arr=="CC","1|1", vcf_arr))
#save re-formated file into vcf file
np.savetxt(f"03vcffile_LD/run{runs}_n400_{num_target}loci_s10_p3/{input_file}.vcf", vcf_arr_sub,fmt="%s", delimiter="\t")

##LE
#vcf_arr=np.loadtxt(f'01hap_file/{input_file}', usecols=tuple(dict_sample['n400_p3']), dtype=str)
#vcf_arr=np.hstack([merged_arr_LE, vcf_arr[:,4:]])
#replace genotype into vcf format
#vcf_arr_sub=np.where(vcf_arr == "AA", "0|0", np.where(vcf_arr=="CC","1|1", vcf_arr))
#save re-formated file into vcf file
#np.savetxt(f"04vcffile_LE/run{runs}_n400_{num_target}loci_s10_p3/{input_file}.vcf", vcf_arr_sub,fmt="%s", delimiter="\t")