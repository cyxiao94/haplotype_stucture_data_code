
for file in `ls /home/vetlinux05/Changyi/04Dsim_SA_flucVScons/06slim_2025/01simulation_10runs/01generate_founder_population/03vcffile_LD/`
do
echo "python extract_chr2L.py $file"
done | parallel -j 40
