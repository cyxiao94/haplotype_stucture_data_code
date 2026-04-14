library(data.table)
library(poolSeq)



#load sample information
#1. get sample information----
all_samples <- fread("00data/header", header = FALSE)$V1


#2. load data----
para <- "biallelic.max15.min85.SA.hfluc_hcons.F0-F80.flt"
print(para)

file_af <- paste0('./00data/',para,".af.txt.gz")
file_cov <- paste0('./00data/',para,".cov.txt.gz")

dt_af <- fread(file_af, header = F, na.strings = c("NA","."))
setnames(dt_af, c("CHR","POS","REF","ALT", all_samples))
dt_cov <- fread(file_cov, header = F)
setnames(dt_cov, c("CHR","POS",all_samples))


#3. fluctuating population(2818)----
print("calculate Ne for fluctuating population")
freq <- dt_af[, grep("anc|28\\-18", names(dt_af)), with = FALSE]
colnames(freq) <- sapply(strsplit(colnames(freq), "_"), function(x) {
    paste(tail(x, 2), collapse = "_")
  })
freq[, (names(freq)) := lapply(.SD, as.numeric)]
freq[,CHR := dt_af$CHR]

cov <- dt_cov[, grep("anc|28\\-18", names(dt_cov)), with = FALSE]
colnames(cov) <- sapply(strsplit(colnames(cov), "_"), function(x) {
  paste(tail(x, 2), collapse = "_")
})
cov[, (names(cov)) := lapply(.SD, as.numeric)]
cov[, CHR := dt_cov$CHR]

##1. remove 0 variance site and site contain NA
#var <-apply(freq, 1, function(x) var(na.omit(x)))
cols <-setdiff(names(freq), c("CHR"))
freq[, var := apply(.SD, 1, function(x) var(x, na.rm = TRUE)), .SDcols = cols]
ind1 <- freq[var == 0 | is.na(var) , which = TRUE]
ind2 <- which(!complete.cases(freq[, ..cols]))
ind <- unique(c(ind1, ind2))
if (length(ind) !=0){
  freq <- freq[-ind,]
  cov <- cov[-ind,]
  }

## estimate Ne
###autosome
med_Ne_auto <- NULL
ne <- NULL
input_rep <- 15
for (i in 1:input_rep){
  df_af<-as.data.frame(subset(freq, CHR != "X", select = paste0("F",c(0,80),"_r",i)))
  df_cov<-as.data.frame(subset(cov, CHR != "X", select = paste0("F",c(0,80),"_r",i)))
  #get Ne
  colnames(df_af)<-c("F0","F80"); colnames(df_cov)<-c("F0","F80")
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(df_af)[1]
  ne <- NULL
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 1000)
    ne <- rbind(ne, data.frame(trial = j,
                               ne = estimateNe(p0 = df_af[ind,"F0"], pt = df_af[ind,"F80"], 
                                               #cov0 = 1000, covt = 1000,
                                               cov0 = df_cov[ind,"F0"], covt = df_cov[ind,"F80"],
                                               t = 80, ploidy=2, truncAF=0.01, method=c("P.planI"), poolSize=rep(6000, times=2), Ncensus=12500)))
    off <- off+1
  }
  
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  print(paste0("autosome rep",i))
  print(ne$ne)
  Ne_rep <- round(median(ne$ne))
  med_Ne_auto<-c(med_Ne_auto,Ne_rep)
}

###X chromosome
med_Ne_x <- NULL
ne <- NULL
input_rep <- 15
for (i in 1:input_rep){
  df_af<-as.data.frame(subset(freq, CHR == "X", select = paste0("F",c(0,80),"_r",i)))
  df_cov<-as.data.frame(subset(cov, CHR == "X", select = paste0("F",c(0,80),"_r",i)))
  #get Ne
  colnames(df_af)<-c("F0","F80"); colnames(df_cov)<-c("F0","F80")
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(df_af)[1]
  ne <- NULL
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 1000)
    ne <- rbind(ne, data.frame(trial = j,
                               ne = estimateNe(p0 = df_af[ind,"F0"], pt = df_af[ind,"F80"], cov0 = df_cov[ind,"F0"], covt = df_cov[ind,"F80"],
                                               t = 80, ploidy=2, truncAF=0.05, method=c("P.planI"), poolSize=rep(600, times=2), Ncensus=1250)))
    off <- off+1
  }
  
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  print(paste0("X chrom rep",i))
  print(ne$ne)
  Ne_rep <- round(median(ne$ne))
  med_Ne_x<-c(med_Ne_x,Ne_rep)
}

rm(df_af);rm(df_cov);rm(ne);rm(i);rm(ind);rm(j);rm(nb_obs);rm(NB_trials);rm(Ne_rep);rm(off)

#write file
ne_data<-data.frame(Ne_auto=med_Ne_auto,
                    NE_x=med_Ne_x,
                    group=rep(c("SupergroupI","SupergroupII","SupergroupIII"), each=5),
                    rep1=paste("r",rep(1:5,3),sep=""),
                    rep2=paste("r",1:15,sep=""),
                    treatment="2818")
write.table(ne_data,file="./05Ne/2818.ne", quote=F, sep='\t', row.names = FALSE)


#4. constant population(23)----
print("calculate Ne for constant population")
freq <- dt_af[, grep("anc|23", names(dt_af)), with = FALSE]
colnames(freq) <- sapply(strsplit(colnames(freq), "_"), function(x) {
  paste(tail(x, 2), collapse = "_")
})
freq[, (names(freq)) := lapply(.SD, as.numeric)]
freq[,CHR := dt_af$CHR]

cov <- dt_cov[, grep("anc|23", names(dt_cov)), with = FALSE]
colnames(cov) <- sapply(strsplit(colnames(cov), "_"), function(x) {
  paste(tail(x, 2), collapse = "_")
})
cov[, (names(cov)) := lapply(.SD, as.numeric)]
cov[,CHR := dt_cov$CHR]

##. remove 0 variance site and site contain NA
#var <-apply(freq, 1, function(x) var(na.omit(x)))
cols <-setdiff(names(freq), c("CHR"))
freq[, var := apply(.SD, 1, function(x) var(x, na.rm = TRUE)), .SDcols = cols]
ind1 <- freq[var == 0 | is.na(var) , which = TRUE]
ind2 <- which(!complete.cases(freq[, ..cols]))
ind <- unique(c(ind1, ind2))
if (length(ind) !=0){
  freq <- freq[-ind,]
  cov <- cov[-ind,]
  }



## estimate Ne
###autosome
med_Ne_auto <- NULL
ne <- NULL
input_rep <- 15
for (i in 1:input_rep){
  df_af<-as.data.frame(subset(freq, CHR != "X", select = paste0("F",c(0,80),"_r",i)))
  df_cov<-as.data.frame(subset(cov, CHR != "X", select = paste0("F",c(0,80),"_r",i)))
  #get Ne
  colnames(df_af)<-c("F0","F80"); colnames(df_cov)<-c("F0","F80")
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(df_af)[1]
  ne <- NULL
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 1000)
    ne <- rbind(ne, data.frame(trial = j,
                               ne = estimateNe(p0 = df_af[ind,"F0"], pt = df_af[ind,"F80"], cov0 = df_cov[ind,"F0"], covt = df_cov[ind,"F80"],
                                               t = 80, ploidy=2, truncAF=0.05, method=c("P.planI"), poolSize=rep(600, times=2), Ncensus=1250)))
    off <- off+1
  }
  
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  print(paste0("autosome rep",i))
  print(ne$ne)
  Ne_rep <- round(median(ne$ne))
  med_Ne_auto<-c(med_Ne_auto,Ne_rep)
}

###X chromosome
med_Ne_x <- NULL
ne <- NULL
input_rep <- 15
for (i in 1:input_rep){
  df_af<-as.data.frame(subset(freq, CHR == "X", select = paste0("F",c(0,80),"_r",i)))
  df_cov<-as.data.frame(subset(cov, CHR == "X", select = paste0("F",c(0,80),"_r",i)))
  #get Ne
  colnames(df_af)<-c("F0","F80"); colnames(df_cov)<-c("F0","F80")
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(df_af)[1]
  ne <- NULL
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 1000)
    ne <- rbind(ne, data.frame(trial = j,
                               ne = estimateNe(p0 = df_af[ind,"F0"], pt = df_af[ind,"F80"], cov0 = df_cov[ind,"F0"], covt = df_cov[ind,"F80"],
                                               t = 80, ploidy=2, truncAF=0.05, method=c("P.planI"), poolSize=rep(600, times=2), Ncensus=1250)))
    off <- off+1
  }
  
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  print(paste0("X chrom rep",i))
  print(ne$ne)
  Ne_rep <- round(median(ne$ne))
  med_Ne_x<-c(med_Ne_x,Ne_rep)
}
rm(df_af);rm(df_cov);rm(ne);rm(i);rm(ind);rm(j);rm(nb_obs);rm(NB_trials);rm(Ne_rep);rm(off)

#write file
ne_data<-data.frame(Ne_auto=med_Ne_auto,
                    Ne_x=med_Ne_x,
                    group=rep(c("SupergroupI","SupergroupII","SupergroupIII"), each=5),
                    rep1=paste("r",rep(1:5,3),sep=""),
                    rep2=paste("r",1:15,sep=""),
                    treatment="23")
write.table(ne_data,file="./05Ne/23.ne", quote=F, sep='\t', row.names = FALSE)


