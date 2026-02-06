suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(corrplot))


#1. load AF file----
args<-commandArgs(trailingOnly = T)
para <- args[1]
#para <- "Run10_popsize1250_n400_100loci_s50"
print(para)
inputf <- paste0('./freq_matrix/',para,"_freq_matrix.txt.gz")
#print(inputf)
af<-fread(inputf, header = T)
af_mat<-as.matrix(af[,2:ncol(af)])
rownames(af_mat)<-af$pos
#rm(f)

#2. get sample information----
samples <- colnames(af_mat)
splits<-unlist(strsplit(samples,'_'))
generations <- as.factor(splits[grep("F",splits)])
replicates <- as.factor(splits[grep("rep",splits)])
groups <- as.factor(splits[grep("^p",splits)])
rep_gen <- as.factor(paste(replicates, generations,sep="_"))

#remove invariant site
var <- apply(af_mat, 1, function(x) var(na.omit(x)))
ind<-which(var==0)
if (length(ind) != 0){
  af_mat<-af_mat[-ind,]
}
#3. PCA----
##3.1 perform PCA----
dat <- 2*asin(sqrt(af_mat))
pcadata <- as.data.frame(t(dat))
pca <- prcomp(pcadata, retx=TRUE, center=TRUE, scale. = TRUE)
pca_score <- pca$x[, 1:20] #keep information for the first 20 PCs
pca_var <- get_eigenvalue(pca)[1:20, 2] 
pc1_loadings <- pca$rotation[, 1] #pc1 loadings
rm(af_mat)
rm(var, ind, dat, pcadata)
saveRDS(pca_score, file = paste0('02pca/pca_result/',para,".pca.res.rds"))
saveRDS(pca_var, file = paste0('02pca/pca_result/',para,".pca.var.rds"))
saveRDS(pc1_loadings, file = paste0('02pca/pca_result/',para,".pc1_loadings.rds"))

##3.2 plot PCA result----
col <- colorRampPalette(c("#E7726F","#FBEC92","#79BC81"))(8)
#gen_pop
p1 <- ggplot(as.data.frame(pca_score), aes(x = PC1, y = PC2, col = generations, shape = groups))+
  geom_point()+
  scale_color_manual(values = c("black",col))+
  labs(x=paste("PC1(",round(pca_var[1],2),"%)",sep=""),
       y=paste("PC2(",round(pca_var[2],2),"%)",sep=""))+
  labs(title=para)
ggsave(paste0('02pca/pca_plot/',para,'.pdf'), plot = p1, width = 7, height = 7)


#4. calculate AFC correlation----
##4.1 extract AFC
prefix <- sub("_F[0-9]+$", "", samples)
gen <- sub(".*_F", "", samples)

#identify F0 columns
f0_cols <- samples[gen == "0"]
names(f0_cols) <- prefix[gen == "0"]  # name by pop+rep

dt_afc <- copy(af)
for (i in seq_along(samples)){
  if (gen[i] != "0"){
    dt_afc[[samples[i]]] <- af[[samples[i]]]-af[[f0_cols[prefix[i]]]]
  }
}
dt_afc <- dt_afc[, !f0_cols, with = FALSE]

##4.2 calculate pearson correlation_all_SNPS
afc_cols <- names(dt_afc)[-1]
meta <- data.table(
  samples = afc_cols,
  pop = sub("_Rep.*","", afc_cols),
  rep = sub(".*_Rep([0-9]+)_F.*","\\1",afc_cols),
  gen = sub(".*_F","", afc_cols)
)


afc_cor <- meta[, {
  mat <- as.matrix(dt_afc[, samples, with = FALSE])
  cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  
  #extract upper triangle
  pairs <- which(upper.tri(cor_mat), arr.ind = TRUE)
  
  data.table(
    sample1 = colnames(cor_mat)[pairs[,1]],
    sample2 = colnames(cor_mat)[pairs[,2]],
    r = cor_mat[pairs]
  )
}, by = gen]

afc_cor[, c("grp1", "rep1") := tstrsplit(sample1, "_Rep", keep = 1:2)]
afc_cor[, c("grp2", "rep2") := tstrsplit(sample2, "_Rep", keep = 1:2)]

afc_cor[, group_pair := paste(grp1, grp2, sep="_")]
afc_cor[, comparison := ifelse(grp1 == grp2, "within", "between")]
afc_cor[, c("grp1","grp2","rep1","rep2") := NULL]
fwrite(afc_cor, file =  paste0("03afc_cor/results_allSNPs/", para, "_cor.txt.gz"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

##4.3 calculate pearson correlation_all_SNPS

###extract information from SNPs with 5% largest PC1 loadings(absolute value)
top_n <- ceiling(length(pc1_loadings) * 0.05)
top_vars <- names(sort(abs(pc1_loadings), decreasing = TRUE)[1:top_n])
dt_afc_sig <- dt_afc[pos %in% top_vars]

afc_cols_sig <- names(dt_afc)[-1]
meta <- data.table(
  samples = afc_cols_sig,
  pop = sub("_Rep.*","", afc_cols_sig),
  rep = sub(".*_Rep([0-9]+)_F.*","\\1",afc_cols_sig),
  gen = sub(".*_F","", afc_cols_sig)
)


afc_cor_sig <- meta[, {
  mat <- as.matrix(dt_afc_sig[, samples, with = FALSE])
  cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  
  #extract upper triangle
  pairs <- which(upper.tri(cor_mat), arr.ind = TRUE)
  
  data.table(
    sample1 = colnames(cor_mat)[pairs[,1]],
    sample2 = colnames(cor_mat)[pairs[,2]],
    r = cor_mat[pairs]
  )
}, by = gen]

afc_cor_sig[, c("grp1", "rep1") := tstrsplit(sample1, "_Rep", keep = 1:2)]
afc_cor_sig[, c("grp2", "rep2") := tstrsplit(sample2, "_Rep", keep = 1:2)]

afc_cor_sig[, group_pair := paste(grp1, grp2, sep="_")]
afc_cor_sig[, comparison := ifelse(grp1 == grp2, "within", "between")]
afc_cor_sig[, c("grp1","grp2","rep1","rep2") := NULL]
fwrite(afc_cor_sig, file =  paste0("03afc_cor/results_PC1_sigSNPs/", para, "_cor.txt.gz"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
