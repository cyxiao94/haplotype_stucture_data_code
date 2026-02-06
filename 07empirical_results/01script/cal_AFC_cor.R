suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(corrplot))






#1. get sample information----
all_samples <- fread("00data/header", header = FALSE)$V1

#2. load data----
para <- "biallelic.max15.min85.SA.hfluc_hcons.F0-F80.flt"
print(para)

file_af <- paste0('./00data/',para,".af.txt.gz")
file_cov <- paste0('./00data/',para,".cov.txt.gz")

dt_af <- fread(file_af, header = F)
setnames(dt_af, c("CHR","POS","REF","ALT", all_samples))



#3. fluctuating population(2818)----
##3.1 data filtering----
af_mat<-as.matrix(dt_af[, grep("anc|28-18", names(dt_af)), with = FALSE])
samples <- all_samples[grepl("anc|28-18", all_samples)]


##remove invariant site and variants contain missing value
var <- apply(af_mat, 1, function(x) var(na.omit(x)))
ind<-which(var==0 | is.na(var))
if (length(ind) != 0){
  af_mat <- af_mat[-ind,]
  
}

af_mat <- matrix(as.numeric(af_mat), nrow = nrow(af_mat), ncol = ncol(af_mat), dimnames = dimnames(af_mat))
rm(var,ind)


##3.2. calculate AFC correlation----
### extract AFC
prefix_rep <- sub(".*_r", "", samples)
prefix_gen <- sub(".*_F([0-9]+)_.*", "\\1", samples)

#identify F0 columns
f0_cols <- samples[grepl("F0",samples)]


dt_afc <- data.table()
for (i in seq_along(samples)){
  if (prefix_gen[i] != "0"){
    rep_num <- prefix_rep[i]
    f0_col <- f0_cols[f0_cols %in% paste0("Dsim_SA_anc_F0_r", rep_num)]
    dt_afc[,samples[i] := af_mat[,samples[i]]-af_mat[,f0_col]]
  }
}
rm(prefix_rep, prefix_gen, f0_cols, i, rep_num, f0_col)
#dt_afc <- dt_afc[, !f0_cols, with = FALSE]

## calculate pearson correlation_all_SNPS
afc_cols <- names(dt_afc)
meta <- data.table(
  samples = afc_cols,
  #supergroup = sub("_Rep.*","", afc_cols),
  rep = sub(".*_r","",afc_cols),
  gen = sub(".*_F([0-9]+)_.*","\\1", afc_cols)
)

meta[, supergroup := ifelse(rep %in% seq(1,5), "SupergroupI",
                            ifelse(rep %in% seq(6,10), "SupergroupII",
                                   ifelse(rep %in% seq(11,15), "SupergroupIII", NA)))]

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

afc_cor[, rep1 := tstrsplit(sample1, "_", keep = 5)]
afc_cor[, grp1 := ifelse(rep1 %in% paste0("r",1:5), "SupergroupI",
                         ifelse(rep1 %in% paste0("r",6:10), "SupergroupII","SupergroupIII"))]
afc_cor[, rep2 := tstrsplit(sample2, "_", keep = 5)]
afc_cor[, grp2 := ifelse(rep2 %in% paste0("r",1:5), "SupergroupI",
                         ifelse(rep2 %in% paste0("r",6:10), "SupergroupII","SupergroupIII"))]

afc_cor[, group_pair := paste(grp1, grp2, sep="_")]
afc_cor[, comparison := ifelse(grp1 == grp2, "within", "between")]
afc_cor[, c("grp1","grp2","rep1","rep2") := NULL]
fwrite(afc_cor, file =  paste0("03AFC_cor/2818_AFC_cor.txt.gz"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
rm(af_mat, afc_cor, dt_afc, meta,afc_cols,samples)

#4. constant population(23)----
##4.1 data filtering----
af_mat<-as.matrix(dt_af[, grep("anc|23", names(dt_af)), with = FALSE])
samples <- all_samples[grepl("anc|23", all_samples)]

#rep_gen <- as.factor(paste(replicates, generations,sep="_"))
##remove invariant site and variants contain missing value
var <- apply(af_mat, 1, function(x) var(na.omit(x)))
ind<-which(var==0 | is.na(var))
if (length(ind) != 0){
  af_mat <- af_mat[-ind,]
  
}

af_mat <- matrix(as.numeric(af_mat), nrow = nrow(af_mat), ncol = ncol(af_mat), dimnames = dimnames(af_mat))
rm(var,ind)


##4.2. calculate AFC correlation----
### extract AFC
prefix_rep <- sub(".*_r", "", samples)
prefix_gen <- sub(".*_F([0-9]+)_.*", "\\1", samples)

#identify F0 columns
f0_cols <- samples[grepl("F0",samples)]


dt_afc <- data.table()
for (i in seq_along(samples)){
  if (prefix_gen[i] != "0"){
    rep_num <- prefix_rep[i]
    f0_col <- f0_cols[f0_cols %in% paste0("Dsim_SA_anc_F0_r", rep_num)]
    dt_afc[,samples[i] := af_mat[,samples[i]]-af_mat[,f0_col]]
  }
}
rm(prefix_rep, prefix_gen, f0_cols, i, rep_num, f0_col)
#dt_afc <- dt_afc[, !f0_cols, with = FALSE]

## calculate pearson correlation_all_SNPS
afc_cols <- names(dt_afc)
meta <- data.table(
  samples = afc_cols,
  #supergroup = sub("_Rep.*","", afc_cols),
  rep = sub(".*_r","",afc_cols),
  gen = sub(".*_F([0-9]+)_.*","\\1", afc_cols)
)

meta[, supergroup := ifelse(rep %in% seq(1,5), "SupergroupI",
                            ifelse(rep %in% seq(6,10), "SupergroupII",
                                   ifelse(rep %in% seq(11,15), "SupergroupIII", NA)))]

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

afc_cor[, rep1 := tstrsplit(sample1, "_", keep = 5)]
afc_cor[, grp1 := ifelse(rep1 %in% paste0("r",1:5), "SupergroupI",
                         ifelse(rep1 %in% paste0("r",6:10), "SupergroupII","SupergroupIII"))]
afc_cor[, rep2 := tstrsplit(sample2, "_", keep = 5)]
afc_cor[, grp2 := ifelse(rep2 %in% paste0("r",1:5), "SupergroupI",
                         ifelse(rep2 %in% paste0("r",6:10), "SupergroupII","SupergroupIII"))]

afc_cor[, group_pair := paste(grp1, grp2, sep="_")]
afc_cor[, comparison := ifelse(grp1 == grp2, "within", "between")]
afc_cor[, c("grp1","grp2","rep1","rep2") := NULL]
fwrite(afc_cor, file =  paste0("03AFC_cor/23_AFC_cor.txt.gz"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
rm(af_mat, afc_cor, dt_afc, meta)
