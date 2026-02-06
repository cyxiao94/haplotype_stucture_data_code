library(data.table)
library(poolSeq)

#setwd('~/my_project/SouthAFrica_hotfluc_hotcons/04slim/test6_fst/11Ne/')
#initial parameters
input_rep <- 15

#load data
args <- commandArgs(trailingOnly = TRUE)

popsize <- 1250 #default population size
if ("--popsize" %in% args) {
  popsize <- as.numeric(args[which(args == "--popsize") + 1])
}

input_file <- args[which(args == "--input_file") + 1]
#input_file <- "test_freq.txt"
treatment <- strsplit(input_file, "_freq_matrix.txt.gz")[[1]][1]
f_data <- fread(
  paste("./freq_matrix/", input_file, sep = ""),
  select = c(
    "pos",
    "p1_Rep1_F0", "p1_Rep1_F80",
    "p1_Rep2_F0", "p1_Rep2_F80",
    "p1_Rep3_F0", "p1_Rep3_F80",
    "p1_Rep4_F0", "p1_Rep4_F80",
    "p1_Rep5_F0", "p1_Rep5_F80",
    "p2_Rep1_F0", "p2_Rep1_F80",
    "p2_Rep2_F0", "p2_Rep2_F80",
    "p2_Rep3_F0", "p2_Rep3_F80",
    "p2_Rep4_F0", "p2_Rep4_F80",
    "p2_Rep5_F0", "p2_Rep5_F80",
    "p3_Rep1_F0", "p3_Rep1_F80",
    "p3_Rep2_F0", "p3_Rep2_F80",
    "p3_Rep3_F0", "p3_Rep3_F80",
    "p3_Rep4_F0", "p3_Rep4_F80",
    "p3_Rep5_F0", "p3_Rep5_F80"
  )
)
colnames(f_data) <- c("pos", paste("F", rep(c(0, 80), 15), ".R", rep(1:15, each =
                                                                       2), ".freq", sep = ""))


##1. remove 0 variance site and size contain NA--------------
#remove invariant site
freq=f_data[,2:ncol(f_data)]
var=apply(freq, 1, function(x) var(na.omit(x)))
ind=which(var==0)
if (length(ind) !=0)
  {f_data=f_data[-ind,]}



## estimate Ne---------
med_Ne=NULL
ne=NULL
for (i in 1:input_rep){
  af=as.data.frame(subset(f_data, select = paste("F",c(0,80),".R",i,".freq",sep="")))
  cov=data.frame("F0"=rep(100, nrow(af)),
                    "F80"=rep(100, nrow(af)))
  #get Ne
  colnames(af)<-c("F0","F80"); colnames(cov)<-c("F0","F80")
  NB_trials <- 100
  off <- 1
  nb_obs <- dim(af)[1]
  ne <- NULL
  for(j in 1:NB_trials){
    set.seed(off)
    ind <- sample(x = 1:nb_obs, size = 1000)
    ne <- rbind(ne, data.frame(trial = j,
                               ne = estimateNe(p0 = af[ind,"F0"], pt = af[ind,"F80"], cov0 = cov[ind,"F0"], covt = cov[ind,"F80"],
                                               t = 80, ploidy=2, truncAF=0.05, method=c("P.planI"), poolSize=rep(600, times=2), Ncensus=popsize)))
    off <- off+1
  }
  
  ind <- which(is.na(ne$ne) | ne$ne<0)
  if(length(ind)>0){ne <- ne[-ind, ]}
  Ne_rep <- round(median(ne$ne))
  med_Ne<-c(med_Ne,Ne_rep)
}

rm(af);rm(cov);rm(ne);rm(i);rm(ind);rm(j);rm(nb_obs);rm(NB_trials);rm(Ne_rep);rm(off)

#write file
ne_data<-data.frame(Ne=med_Ne,
                    group=rep(c("A","B","C"), each=5),
                    rep1=paste("r",rep(1:5,3),sep=""),
                    rep2=paste("r",1:15,sep=""),
                    treatment=treatment)
write.table(ne_data,file=paste("./04estimate_Ne/", treatment,".ne",sep=""), quote=F, sep='\t', row.names = FALSE)
