# Process Covid data
rm(list=ls())

library(tidyverse)
library(dplyr)
library(magrittr)
library(readxl)


#################################
# Read Alamar data and metadata
#################################
alamar_dt <- readxl::read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "Alamar" )
alamar_dt <- alamar_dt[,-c(1,3)] %>% as.data.frame()
rownames(alamar_dt) = alamar_dt$Assay
alamar_dt = alamar_dt[,-1] %>% t() %>% as.data.frame()
colnames(alamar_dt) = paste(colnames(alamar_dt), "[ALAMAR]", sep = "")


################################## 
# Read OLINK data
################################## 

olink_dt <- readxl::read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "Olink" )
olink_dt <- olink_dt[,-c(1,3)] %>% as.data.frame()
rownames(olink_dt) = olink_dt$Assay
olink_dt = olink_dt[,-1] %>% t() %>% as.data.frame()
colnames(olink_dt) = paste(colnames(olink_dt), "[OLINK]", sep = "")

###############
# DIA
###############

dia_dt <- read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "DIA" )
dia_dt <- dia_dt[,-c(1,3)] %>% as.data.frame()
rownames(dia_dt) = dia_dt$Assay
dia_dt = dia_dt[,-1] %>% t() %>% as.data.frame()
colnames(dia_dt)
colnames(dia_dt) = paste(colnames(dia_dt), "[DIA]", sep = "")

MBHfile = paste("mutual_best_hits_covid.tsv", sep="")
cat(file=MBHfile, append = FALSE,
    "PLAT1","PLAT2","N","Spearman rho","P-value","BONF","TRAIT1","TRAIT2",sep="\t","\n")

MBH(alamar_dt,olink_dt, "ALAMAR","OLINK" )
MBH(olink_dt,dia_dt, "OLINK","DIA" )
MBH(alamar_dt,dia_dt, "ALAMAR","DIA" )

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,"MBH")
openxlsx::writeData(wb, "GGM", GGM_covid, rowNames = T, colNames=T)

openxlsx::saveWorkbook (wb, file="NETWORKS.xlsx", overwrite=TRUE)


MBH <- function(D1,D2, PLAT1, PLAT2){
  common <- intersect(rownames(D1), rownames(D2))
  D1_intersect <- D1[common,] # give you common rows in data frame 1  
  D2_intersect <- D2[common,] 
  
  rm(common)
  ###################################################################################
  # Calculate correlation between different platforms and select Mutual Best Scores
  ###################################################################################
  
  scor <- cor(D1_intersect,D2_intersect, method = "spearman", use = "pairwise.complete.obs")
  scor[is.na(scor)] = 0
  
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
    k1 = which(colnames(D1_intersect) == bPLAT1[k])
    k2 = which(colnames(D2_intersect) == bPLAT2[k])
    var1 = D1_intersect[,k1]
    var2 = D2_intersect[,k2]
    #res <- outer(var1,var2,Vectorize(function(a,b) cor.test(a,b, method = "spearman")$p.value))
    res = cor.test(var1,var2, method = "spearman", exact = FALSE)
    coragain = res$estimate
    if (abs(coragain-bCorr1[k])>1e-5) {stop("something went wrong")}
    pcorr = res$p.value
    ncorr = length(which(!is.na(var1+var2)))
    #cor(var1,var2, use = "pairwise.complete.obs", method = "spearman") 
    bonf = pcorr < (0.05/(dim(scor)[1]*dim(scor)[2]))
    cat(file = MBHfile, append = TRUE, 
        PLAT1, PLAT2, ncorr, bCorr1[k], pcorr, bonf, bPLAT1[k], bPLAT2[k], sep = "\t", "\n")
    
  }
}