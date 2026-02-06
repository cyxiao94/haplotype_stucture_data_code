import sys

fi_name=sys.argv[1]
fo_name=fi_name.replace("03vcffile_LD","04vcffile_LE")

fo=open(fo_name,"w")

with open(fi_name,"r") as fi, open(fo_name,"w") as fo:
    i=1
    for line in fi:
        if line.startswith("#"):
            fo.write(line)
        else:
            item=line.rstrip().split('\t')
            chrom=item[0];pos=item[1]
            new_pos=i
            fo.write(f"{chrom}\t{new_pos}\t"+"\t".join(item[2:])+"\n")
            i = i + 1
   