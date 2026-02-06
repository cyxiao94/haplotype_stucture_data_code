suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(corrplot))



#load sample information
#1. get sample information----
all_samples <- fread("00data/header", header = FALSE)$V1


#2. load data----
para <- "biallelic.max15.min85.SA.hfluc_hcons.F0-F80.flt"
print(para)

file_af <- paste0('./00data/',para,".af.txt.gz")
file_cov <- paste0('./00data/',para,".cov.txt.gz")

dt_af <- fread(file_af, header = F)
setnames(dt_af, c("CHR","POS","REF","ALT", all_samples))
#dt_cov <- fread(file_cov, header = F)
#setnames(dt_cov, c("CHR","POS",samples))



#3. fluctuating population(2818)----
##3.1 data filtering----
af_mat<-as.matrix(dt_af[, grep("anc|28-18", names(dt_af)), with = FALSE])
samples <- all_samples[grepl("anc|28-18", all_samples)]
treatment <- sapply(strsplit(samples, "_"), `[`, 3)
splits<-unlist(strsplit(samples,'_'))
generations <- as.factor(splits[grep("F",splits)])
replicates <- as.factor(splits[grep("r",splits)])
supergroup <- ifelse(replicates %in% paste0("r", 1:5), "SupergroupI",
                     ifelse(replicates %in% paste0("r", 6:10), "SupergroupII",
                            ifelse(replicates %in% paste0("r", 11:15), "SupergroupIII", NA)))

### remove invariant site and variants contain missing value
var <- apply(af_mat, 1, function(x) var(na.omit(x)))
ind<-which(var==0 | is.na(var))
if (length(ind) != 0){
  af_mat<-af_mat[-ind,]
}

af_mat <- matrix(as.numeric(af_mat), nrow = nrow(af_mat), ncol = ncol(af_mat), dimnames = dimnames(af_mat))

##3.2 perform PCA----
dat <- 2*asin(sqrt(af_mat))
pcadata <- as.data.frame(t(dat))
pca <- prcomp(pcadata, retx = TRUE, center = TRUE, scale. = TRUE)
pca_score <- pca$x[, 1:20] #keep information for the first 20 PCs
pca_var <- get_eigenvalue(pca)[1:20, 2] 
pc1_loadings <- pca$rotation[, 1] #pc1 loadings
saveRDS(pca_score, file = paste0("02pca/2818.pca.res.rds"))
saveRDS(pca_var, file = paste0("02pca/2818.pca.var.rds"))
saveRDS(pc1_loadings, file = paste0("02pca/2818.pc1_loadings.rds"))
rm(af_mat)
rm(var, ind, dat, pcadata)


##3.3 plot PCA result----
col <- colorRampPalette(c("#E7726F","#FBEC92","#79BC81"))(8)
#gen_pop
p1 <- ggplot(as.data.frame(pca_score), aes(x = PC1, y = PC2, col = generations, shape = supergroup))+
  geom_point()+
  scale_color_manual(values = c("black",col))+
  labs(x=paste("PC1(",round(pca_var[1],2),"%)",sep=""),
       y=paste("PC2(",round(pca_var[2],2),"%)",sep=""))+
  labs(title="PCA_2818")
ggsave(paste0("02pca/2818_pca.pdf"), plot = p1, width = 7, height = 7)
rm(pca, pca_score, pca_var, pc1_loadings)
rm(samples, treatment, splits, generations, replicates,supergroup)

#4. constant population(23)----
##4.1 data filtering----
af_mat<-as.matrix(dt_af[, grep("anc|23", names(dt_af)), with = FALSE])
samples <- all_samples[grepl("anc|23", all_samples)]
treatment <- sapply(strsplit(samples, "_"), `[`, 3)
splits<-unlist(strsplit(samples,'_'))
generations <- as.factor(splits[grep("F",splits)])
replicates <- as.factor(splits[grep("r",splits)])
supergroup <- ifelse(replicates %in% paste0("r", 1:5), "SupergroupI",
                     ifelse(replicates %in% paste0("r", 6:10), "SupergroupII",
                            ifelse(replicates %in% paste0("r", 11:15), "SupergroupIII", NA)))

## remove invariant site and variants contain missing value
var <- apply(af_mat, 1, function(x) var(na.omit(x)))
ind<-which(var==0 | is.na(var))
if (length(ind) != 0){
  af_mat<-af_mat[-ind,]
  
}

af_mat <- matrix(as.numeric(af_mat), nrow = nrow(af_mat), ncol = ncol(af_mat), dimnames = dimnames(af_mat))

##4.2 perform PCA----
dat <- 2*asin(sqrt(af_mat))
pcadata <- as.data.frame(t(dat))
pca <- prcomp(pcadata, retx = TRUE, center = TRUE, scale. = TRUE)
pca_score <- pca$x[, 1:20] #keep information for the first 20 PCs
pca_var <- get_eigenvalue(pca)[1:20, 2] 
pc1_loadings <- pca$rotation[, 1] #pc1 loadings
saveRDS(pca_score, file = paste0("02pca/23.pca.res.rds"))
saveRDS(pca_var, file = paste0("02pca/23.pca.var.rds"))
saveRDS(pc1_loadings, file = paste0("02pca/23.pc1_loadings.rds"))
rm(af_mat)
rm(var, ind, dat, pcadata)

##4.3 plot PCA result----
col <- colorRampPalette(c("#E7726F","#FBEC92","#79BC81"))(8)
#gen_pop
p1 <- ggplot(as.data.frame(pca_score), aes(x = PC1, y = PC2, col = generations, shape = supergroup))+
  geom_point()+
  scale_color_manual(values = c("black",col))+
  labs(x=paste("PC1(",round(pca_var[1],2),"%)",sep=""),
       y=paste("PC2(",round(pca_var[2],2),"%)",sep=""))+
  labs(title="PCA_23")
ggsave(paste0("02pca/23_pca.pdf"), plot = p1, width = 7, height = 7)
rm(samples, treatment, splits, generations, replicates,supergroup)
rm(pca, pca_score, pca_var, pc1_loadings)

#5. combine two environment----
##5.1 data filtering----
af_mat<-as.matrix(dt_af[, grep("anc|28-18|23", names(dt_af)), with = FALSE])
samples <- all_samples[grepl("anc|28-18|23", all_samples)]
treatment <- sapply(strsplit(samples, "_"), `[`, 3)
splits<-unlist(strsplit(samples,'_'))
generations <- as.factor(splits[grep("F",splits)])
replicates <- as.factor(splits[grep("r",splits)])
supergroup <- ifelse(replicates %in% paste0("r", 1:5), "SupergroupI",
                     ifelse(replicates %in% paste0("r", 6:10), "SupergroupII",
                            ifelse(replicates %in% paste0("r", 11:15), "SupergroupIII", NA)))
treat_gen <- as.factor(paste(treatment,generations, sep = "_"))
#remove invariant site and variants contain missing value
var <- apply(af_mat, 1, function(x) var(na.omit(x)))
ind<-which(var==0 | is.na(var))
if (length(ind) != 0){
  af_mat<-af_mat[-ind,]
  
}
af_mat <- matrix(as.numeric(af_mat), nrow = nrow(af_mat), ncol = ncol(af_mat), dimnames = dimnames(af_mat))



##5.2 perform PCA----
dat <- 2*asin(sqrt(af_mat))
pcadata <- as.data.frame(t(dat))
pca <- prcomp(pcadata, retx = TRUE, center = TRUE, scale. = TRUE)
pca_score <- pca$x[, 1:20] #keep information for the first 20 PCs
pca_var <- get_eigenvalue(pca)[1:20, 2] 
pc1_loadings <- pca$rotation[, 1] #pc1 loadings
saveRDS(pca_score, file = paste0("02pca/2818_23.pca.res.rds"))
saveRDS(pca_var, file = paste0("02pca/2818_23.pca.var.rds"))
saveRDS(pc1_loadings, file = paste0("02pca/2818_23.pc1_loadings.rds"))
rm(af_mat)
rm(var, ind, dat, pcadata)


##5.3 plot PCA result----
col_cons<-colorRampPalette(colors = c("#203764","#5AB6EA"))
col_fluc<-colorRampPalette(colors = c("#833C0C","#FBC21A"))
col1<-c("black",col_cons(8), col_fluc(8)) #anc, cons, fluc 

p_pca12<-ggplot(pca_score,aes(x=PC1,y=PC2, col=treat_gen, shape=supergroup))+
  geom_point()+
  annotate("text",x=1000,y=-200, label="supergroup I", color="darkred")+
  annotate("text",x=0,y=-500, label="supergroup II", color="darkred")+
  annotate("text",x=0,y=750, label="supergroup III", color="darkred")+
  scale_color_manual(values = col1)+
  guides(shape = "none")+
  xlab(paste("PC1(",round(pca_var[1], 2),"%)",sep=""))+
  ylab(paste("PC2(",round(pca_var[2], 2),"%)",sep=""))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title=element_text(size=15))

#plot(p_pca12)
ggsave("02pca/2818_23_PCA12.pdf", plot=p_pca12, width = 8, height = 6)
rm(samples, treatment, splits, generations, replicates,supergroup,treat_gen)
rm(pca, pca_score, pca_var, pc1_loadings)
