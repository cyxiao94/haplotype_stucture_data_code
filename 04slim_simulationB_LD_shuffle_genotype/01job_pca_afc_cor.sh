mkdir -p 02pca/pca_result
mkdir -p 02pca/pca_plot
mkdir -p 03afc_cor/results


for file in `ls freq_matrix/*freq*`; do
    base=$(basename $file _freq_matrix.txt.gz)
    #echo "Processing $base"
    echo "Rscript 00scripts/perform_PCA_afc_cor.R $base"
done | parallel -j 20
