rm(list=ls())

library(tidyverse)
library(dplyr)
library(RNOmni)
library(GeneNet)
library(magrittr)
library(readxl)

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

metadata <- readxl::read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "Metadata" )
###############
# Read Alamar data and metadata
###############

alamar_dt <- read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "Alamar" )
alamar_dt <- alamar_dt[,-c(1,3)] %>% as.data.frame()
rownames(alamar_dt) = alamar_dt$Assay
alamar_dt = alamar_dt[,-1] %>% t() %>% as.data.frame()
colnames(alamar_dt) = paste(colnames(alamar_dt), "[ALAMAR]", sep = "")

# NA removed
alamar_mt <- as.matrix(alamar_dt)
alamar_anno <- data.frame(TRAITID = colnames(alamar_dt), TRAIT = colnames(alamar_dt), 
                          PLAT = "ALAMAR", SHORTNAME = colnames(alamar_dt) )



###########################
# Graphical Gaussian Models (GeneNet)
# ggm.estimate.pcor
# gene names or metabolites in this case should be in columns and participants as row. 
# Typical gene expression matrix has gene names in rows and samples in columns. If the data is
# in this format then it needs to be transposed before applying GGM
###########################
# Compute Partial Correlations and Select Relevant Edges

pcor = ggm.estimate.pcor(alamar_mt, method="dynamic") 
alamar.edges = network.test.edges(pcor, direct = TRUE)

#use strongest 200 edges
#nmr.net = extract.network(nmr.edges, method.ggm = "number", cutoff.ggm = 200)
alamar.net = extract.network(alamar.edges)
node.labels = as.character(colnames(alamar_mt))

matching_indices_1 <- match(alamar.net$node1, row_number(node.labels))
matching_indices_2 <- match(alamar.net$node2, row_number(node.labels))

alamar.net$nodelabel1 <- node.labels[matching_indices_1]
alamar.net$nodelabel2 <- node.labels[matching_indices_2]

alamar.net$PLAT <- "ALAMAR"
colnames(alamar.net) <- c("COR", "TRAIT1","TRAIT2","PVALUE","QVAL","PROB","LOG.SPVAR",
                          "PVAL.DIR","QVAL.DIR","PROB.DIR","DIRECTION",
                          "TRAITID1","TRAITID2", "PLAT")
write.table(alamar.net, file="GGMs/Alamar_GGM_Output.txt", sep="\t", col.names = NA)
save(alamar.net, file="Alamar_GGM.rda")

###############
# Olink
###############

olink_dt <- read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "Olink" )
olink_dt <- olink_dt[,-c(1,3)] %>% as.data.frame()
rownames(olink_dt) = olink_dt$Assay
olink_dt = olink_dt[,-1] %>% t() %>% as.data.frame()
colnames(olink_dt) = paste(colnames(olink_dt), "[OLINK]", sep = "")

olink_mt <- as.matrix(olink_dt)
olink_anno <- data.frame(TRAITID = colnames(olink_dt), TRAIT = colnames(olink_dt), 
                          PLAT = "OLINK", SHORTNAME = colnames(olink_dt) )


pcor = ggm.estimate.pcor(olink_mt, method="dynamic") 
olink.edges = network.test.edges(pcor, direct = TRUE)


olink.net = extract.network(olink.edges)
node.labels = as.character(colnames(olink_mt))

matching_indices_1 <- match(olink.net$node1, row_number(node.labels))
matching_indices_2 <- match(olink.net$node2, row_number(node.labels))

olink.net$nodelabel1 <- node.labels[matching_indices_1]
olink.net$nodelabel2 <- node.labels[matching_indices_2]

olink.net$PLAT <- "OLINK"
colnames(olink.net) <- c("COR", "TRAIT1","TRAIT2","PVALUE","QVAL","PROB","LOG.SPVAR",
                         "PVAL.DIR","QVAL.DIR","PROB.DIR","DIRECTION",
                         "TRAITID1","TRAITID2", "PLAT")
write.table(olink.net, file="GGMs/Olink_GGM_Output.txt", sep="\t", col.names = NA)
save(olink.net, file="Olink_GGM.rda")

###############
# DIA
###############

dia_dt <- read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "DIA" )
dia_dt <- dia_dt[,-c(1,3)] %>% as.data.frame()
rownames(dia_dt) = dia_dt$Assay
dia_dt = dia_dt[,-1] %>% t() %>% as.data.frame()
colnames(dia_dt) = paste(colnames(dia_dt), "[DIA]", sep = "")

dia_mt <- as.matrix(dia_dt)

dia_anno <- data.frame(TRAITID = colnames(dia_dt), TRAIT = colnames(dia_dt), 
                         PLAT = "DIA", SHORTNAME = colnames(dia_dt) )

pcor = ggm.estimate.pcor(dia_mt, method="dynamic") 
dia.edges = network.test.edges(pcor, direct = TRUE)


dia.net = extract.network(dia.edges)
node.labels = as.character(colnames(dia_mt))

matching_indices_1 <- match(dia.net$node1, row_number(node.labels))
matching_indices_2 <- match(dia.net$node2, row_number(node.labels))

dia.net$nodelabel1 <- node.labels[matching_indices_1]
dia.net$nodelabel2 <- node.labels[matching_indices_2]

dia.net$PLAT <- "DIA"
colnames(dia.net) <- c("COR", "TRAIT1","TRAIT2","PVALUE","QVAL","PROB","LOG.SPVAR",
                       "PVAL.DIR","QVAL.DIR","PROB.DIR","DIRECTION",
                       "TRAITID1","TRAITID2", "PLAT")
write.table(dia.net, file="GGMs/DIA_GGM_Output.txt", sep="\t", col.names = NA)
save(dia.net, file="DIA_GGM.rda")

###### combine networks into one file #######
dia.net %>% str()
olink.net %>% str()
alamar.net %>% str()
dia.net$TRAIT1 = dia.net$TRAITID1
dia.net$TRAIT2 = dia.net$TRAITID2
dia.net.req = dia.net %>% select( "TRAITID1","TRAITID2" , "TRAIT1","TRAIT2","COR", "PVALUE", "PLAT")

olink.net$TRAIT1 = olink.net$TRAITID1
olink.net$TRAIT2 = olink.net$TRAITID2
olink.net.req = olink.net %>% select( "TRAITID1","TRAITID2" , "TRAIT1","TRAIT2","COR", "PVALUE", "PLAT")

alamar.net$TRAIT1 = alamar.net$TRAITID1
alamar.net$TRAIT2 = alamar.net$TRAITID2
alamar.net.req = alamar.net %>% select( "TRAITID1","TRAITID2" , "TRAIT1","TRAIT2","COR", "PVALUE", "PLAT")

GGM_covid = dia.net.req %>% rbind(olink.net.req) %>% rbind(alamar.net.req)

ANNO_covid = alamar_anno  %>% rbind(olink_anno ) %>% rbind(dia_anno ) %>% unique()


##############################################
# MBH 
##############################################
MBHfile = paste("mutual_best_hits_covid.tsv", sep="")
cat(file=MBHfile, append = FALSE,
    "PLAT1","PLAT2","N","Spearman rho","P-value","BONF","TRAIT1","TRAIT2",sep="\t","\n")

MBH(alamar_dt,olink_dt, "ALAMAR","OLINK" )
MBH(olink_dt,dia_dt, "OLINK","DIA" )
MBH(alamar_dt,dia_dt, "ALAMAR","DIA" )

MBH_covid <- read.table(file="mutual_best_hits_covid.tsv", header=T, check.names = F, sep="\t")
colnames(MBH_covid)
MBH_covid$TRAITID1 <- MBH_covid$TRAIT1
MBH_covid$TRAITID2 <- MBH_covid$TRAIT2
colnames(MBH_covid)[colnames(MBH_covid) == "P-value"] <- "PVALUE"
colnames(MBH_covid)[colnames(MBH_covid) == "Spearman rho"] <- "COR"

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,"GGM")
openxlsx::addWorksheet(wb,"ANNO")
openxlsx::addWorksheet(wb,"MBP")
openxlsx::writeData(wb, "GGM", GGM_covid, rowNames = T, colNames=T)
openxlsx::writeData(wb, "ANNO", ANNO_covid, rowNames = T, colNames=T)
openxlsx::writeData(wb, "MBP", MBH_covid, rowNames = T, colNames=T)
openxlsx::saveWorkbook (wb, file="NETWORKS.xlsx", overwrite=TRUE)


