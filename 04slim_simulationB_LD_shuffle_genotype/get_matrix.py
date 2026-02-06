import sys
import time
import argparse
import gzip

parser = argparse.ArgumentParser(description="parameters")
parser.add_argument("--num_target", type=str, required=True, help="Number of selection target")
#parser.add_argument("--fitness_mean", type=str, required=True, help="Number of selection target")
parser.add_argument("--strength", type=str, required=True, help="Total selection strength")
parser.add_argument("--popsize", type=str, required=True, help="Census population size")
parser.add_argument("--num_sample", type=str, required=True, help="Number of foudner haplotypes")
parser.add_argument("--runs", type=str, required=True, help="Simulation runs")
args = parser.parse_args()
#dir=sys.argv[1]

#populatin settings
supergroup=["p1","p2","p3"]
reps=range(1,6)
generations=[0,10,20,30,40,50,60,70,80]

#print("Merging results from 10 replicates")
#load variant information
#print("##load snp information")
start=time.time()
snps=[]
snp_traj={}
with open('snps', 'r') as f:
    for line in f:
        item=line.strip()
        pos=str(int(item)-1)
        snps.append(pos)
        snp_traj[pos]=["0.0"]*len(generations)*len(supergroup)*len(reps)
#print("took "+str(time.time()-start)+"s")


#print("##load snp trajectories")
f_info=open(f'freq_matrix/Run{args.runs}_popsize{args.popsize}_n{args.num_sample}_{args.num_target}loci_s{args.strength}_info.txt','w')
f_info.write(f"Run\tPopsize\tNumber_sample\tNumber_target\tSupergroup\tReplicate\tseed\n")
#load trajectory
for pop in [1,2,3]:
    for rep in reps:
        start=time.time()
        input=f"freq_matrix/run{args.runs}_popsize{args.popsize}_n{args.num_sample}_{args.num_target}loci_s{args.strength}_p{pop}_rep{rep}.log"  #e.g. 10loci_n25/Run1_rep1.log
        with open(input,'r') as f:
            #read seed
            f.readline() #skip comment
            seed_line = f.readline()
            f_info.write(f"{args.runs}\t{args.popsize}\tn{args.num_sample}\t{args.num_target}loci\ts{args.strength}\tp{pop}\trep{rep}\t{seed_line.strip()}\n")

            for line in f:
                item=line.strip().split('\t')
                if len(item)==3:
                    gen=item[0]
                    pos=item[1]
                    freq=item[2]
                    index=int(len(generations)*len(reps)*(pop-1)+len(generations)*(rep-1)+int(gen)/10) 
                    snp_traj[pos][index]=freq
                    
        end=time.time()
    #print("load "+input+" took "+str(end-start)+"s")
f_info.close()
#print("##write matrix")
start=time.time()
sample_list=[]

for pop in supergroup:
    for rep in reps:
        for gen in generations:
            sample_list.append(f"{pop}_Rep{rep}_F{gen}")

#write matrix
output_file=f'freq_matrix/Run{args.runs}_popsize{args.popsize}_n{args.num_sample}_{args.num_target}loci_s{args.strength}_freq_matrix.txt.gz'
#fo=gzip.open(f'{args.num_target}loci_n{args.num_sample}/Fitness{args.fitness_mean}_run{args.runs}_freq_matrix.txt.gz','wt')
with gzip.open(output_file,'wt') as fo:
    #write header
    fo.write('pos\t'+'\t'.join(sample_list)+'\n')
    #write SNP rows
    for snp in snps:
        if any(x != "0.0" for x in snp_traj[snp]):
            fo.write(snp+'\t'+'\t'.join(snp_traj[snp])+'\n')
#fo.close()

#print("took "+str(time.time()-start)+"s")
#print(f"Merging Finished! Took {(time.time()-start)/60:.2f} minutes!!")
