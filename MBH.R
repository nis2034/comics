# Process UK BioBank data
rm(list=ls())

library(tidyverse)
library(dplyr)
#library(RNOmni)

#######
#setwd("/home/thabib/Documents/UKB/")
setwd("/Users/tah4012/Documents/Research/UKB/")

#################################
# Read NMR data and metadata
#################################

nmr_dt <- read.table(file="NMR/nmr_results.tsv", header=T, sep="\t", check.names = F, row.names = 1)
nmr_metadata <- read.table(file="NMR/NMR_Pheno.txt", header=T, sep="\t", check.names = F)
col.order <- nmr_metadata$FieldID
nmr_dt <- subset(nmr_dt, select = nmr_metadata$FieldID)
nmr_dt <- nmr_dt[ , col.order]
nmr_metadata$Description <- gsub(" ", "_", nmr_metadata$Description)
colnames(nmr_dt) <- nmr_metadata$Description
nmr <- nmr_dt
rm(nmr_dt)

################################## 
# Read OLINK data
################################## 

olink <- read.table(file="OLINK/results_olink_instance_0.tsv", header=T, sep="\t", check.names = F, row.names = 1)
colnames(olink) <- toupper(gsub("olink_instance_0.","",colnames(olink)))

###############################
# Read CLINICAL data
###############################

clin <- read.table(file="CLINICAL/clinical_results.tsv", header=T, sep="\t", check.names = F, row.names = 1)
clin_info <- read.table(file="CLINICAL/clinical.txt", header=T, sep="\t", check.names = F)
col.order <- clin_info$FieldID
clin <- subset(clin, select = clin_info$FieldID)
clin <- clin[ , col.order]
clin_info$Description <- gsub(" ", "_", clin_info$Description)
colnames(clin) <- clin_info$Description

################################################
# Common participants between NMR and OLINK
################################################

common <- intersect(rownames(nmr), rownames(olink))
nmr_intersect <- nmr[common,] # give you common rows in data frame 1  
olink_intersect <- olink[common,] 

# common <- merge(nmr,olink, by = "row.names")
# rownames(common) <- common$Row.names
# common <- common[,c(-1)]

rm(common)
################################################
# Common participants between NMR and CLIN
################################################

common <- intersect(rownames(nmr), rownames(clin))
nmr_intersect <- nmr[common,] # give you common rows in data frame 1  
clin_intersect <- clin[common,] 

# common <- merge(nmr,clin, by = "row.names")
# rownames(common) <- common$Row.names
# common <- common[,c(-1)]
rm(common)
################################################
# Common participants between OLINK and CLIN
################################################

common <- intersect(rownames(olink), rownames(clin))
olink_intersect <- olink[common,] # give you common rows in data frame 1  
clin_intersect <- clin[common,] 

# common <- merge(olink,clin, by = "row.names")
# rownames(common) <- common$Row.names
# common <- common[,c(-1)]
rm(common)
###################################################################################
# Calculate correlation between different platforms and select Mutual Best Scores
###################################################################################

scor <- cor(nmr_intersect,olink_intersect, method = "spearman", use = "pairwise.complete.obs")
#scor <- cor(nmr_intersect,clin_intersect, method = "spearman", use = "pairwise.complete.obs")
#scor <- cor(olink_intersect,clin_intersect, method = "spearman", use = "pairwise.complete.obs")

scor[is.na(scor)] = 0
MBHfile = paste("mutual_best_hits_NMR_OLINK.tsv", sep="")
#MBHfile = paste("mutual_best_hits_NMR_CLIN.tsv", sep="")
#MBHfile = paste("mutual_best_hits_OLINK_CLIN.tsv", sep="")

cat(file=MBHfile, append = FALSE,
    "PLAT1","PLAT2","N","Spearman rho","P-value","BONF","TRAIT1","TRAIT2",sep="\t","\n")
n1 <- dim(scor)[1]
n2 <- dim(scor)[2]
# find the best correlation in PLAT2 starting from PLAT1
maxcorr1 = rep(NA, n1)
maxcorrval1 = rep(NA, n1)
for (j1 in c(1:n1)) {
  maxcorr1[j1] = which.max( scor[j1,] )
  maxcorrval1[j1] = scor[j1,maxcorr1[j1]]  
}

# find the best correlation in PLAT1 starting from PLAT2
maxcorr2 = rep(NA, n2)
maxcorrval2 = rep(NA, n2)
for (j2 in c(1:n2)) {
  maxcorr2[j2] = which.max( scor[,j2] )
  maxcorrval2[j2] = scor[maxcorr2[j2],j2]  
}

ix12 = which(maxcorr2[maxcorr1[]] == c(1:n1))
ix21 = which(maxcorr1[maxcorr2[]] == c(1:n2))

if(length(ix12) != length(ix21)) {stop("Incosistent results")}

bPLAT1 = rownames(scor)[ix12]
bPLAT2 = colnames(scor)[maxcorr1[ix12]]
bCorr1 = maxcorrval1[ix12]
bCorr2 = maxcorrval2[ix21]

for (k in c(1:length(ix12))) {
  # compute correlation again to obtain a p-value and the number of non-NA samples
  k1 = which(colnames(nmr_intersect) == bPLAT1[k])
  k2 = which(colnames(olink_intersect) == bPLAT2[k])
  var1 = nmr_intersect[,k1]
  var2 = olink_intersect[,k2]
  #res <- outer(var1,var2,Vectorize(function(a,b) cor.test(a,b, method = "spearman")$p.value))
  res = cor.test(var1,var2, method = "spearman", exact = FALSE)
  coragain = res$estimate
  if (abs(coragain-bCorr1[k])>1e-5) {stop("something went wrong")}
  pcorr = res$p.value
  ncorr = length(which(!is.na(var1+var2)))
  #cor(var1,var2, use = "pairwise.complete.obs", method = "spearman") 
  bonf = pcorr < (0.05/(dim(scor)[1]*dim(scor)[2]))
   cat(file = MBHfile, append = TRUE, 
       "NMR", "OLINK", ncorr, bCorr1[k], pcorr, bonf, bPLAT1[k], bPLAT2[k], sep = "\t", "\n")
       #"NMR", "CLIN", ncorr, bCorr1[k], pcorr, bonf, bPLAT1[k], bPLAT2[k], sep = "\t", "\n")
       #"OLINK", "CLIN", ncorr, bCorr1[k], pcorr, bonf, bPLAT1[k], bPLAT2[k], sep = "\t", "\n")
}
