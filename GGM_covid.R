# Process Covid study data
# clear the workspace
rm(list=ls())

library(tidyverse)
library(dplyr)
library(RNOmni)
library(GeneNet)
library(magrittr)
library(readxl)


metadata <- readxl::read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "Metadata" )
###############
# Read Alamar data and metadata
###############

alamar_dt <- read_excel("Data/covid_merge_7datasets_full.xlsx", sheet = "Alamar" )
alamar_dt <- alamar_dt[,-c(1,3)] %>% as.data.frame()
rownames(alamar_dt) = alamar_dt$Assay
alamar_dt = alamar_dt[,-1] %>% t() %>% as.data.frame()
# alamar_dt$master_id = rownames(alamar_dt)
# alamar_dt <- alamar_dt %>%
#   select(master_id, everything())
# rownames(alamar_dt) = NULL

# NA removed
alamar_mt <- as.matrix(alamar_dt)


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

olink_mt <- as.matrix(olink_dt)


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

dia_mt <- as.matrix(dia_dt)


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
dia.net$TRAIT1 =paste(dia.net$TRAITID1, "[",dia.net$PLAT,"]", sep = "")
dia.net$TRAIT2 =paste(dia.net$TRAITID2, "[",dia.net$PLAT,"]", sep = "")
dia.net.req = dia.net %>% select( "TRAITID1","TRAITID2" , "TRAIT1","TRAIT2","COR", "PVALUE", "PLAT")

olink.net$TRAIT1 =paste(olink.net$TRAITID1, "[",olink.net$PLAT,"]", sep = "")
olink.net$TRAIT2 =paste(olink.net$TRAITID2, "[",olink.net$PLAT,"]", sep = "")
olink.net.req = olink.net %>% select( "TRAITID1","TRAITID2" , "TRAIT1","TRAIT2","COR", "PVALUE", "PLAT")

alamar.net$TRAIT1 =paste(alamar.net$TRAITID1, "[",alamar.net$PLAT,"]", sep = "")
alamar.net$TRAIT2 =paste(alamar.net$TRAITID2, "[",alamar.net$PLAT,"]", sep = "")
alamar.net.req = alamar.net %>% select( "TRAITID1","TRAITID2" , "TRAIT1","TRAIT2","COR", "PVALUE", "PLAT")

GGM_covid = dia.net.req %>% rbind(olink.net.req) %>% rbind(alamar.net.req)
ANNO_covid = GGM_covid[,c("TRAITID1", "TRAIT1","PLAT", SHORTNAME = "TRAIT1")]
names(ANNO_covid) = c("TRAITID",	"TRAIT",	"PLAT",	"SHORTNAME")




wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb,"GGM")
openxlsx::addWorksheet(wb,"ANNO")
openxlsx::writeData(wb, "GGM", GGM_covid, rowNames = T, colNames=T)
openxlsx::writeData(wb, "ANNO", ANNO_covid, rowNames = T, colNames=T)
openxlsx::saveWorkbook (wb, file="NETWORKS.xlsx", overwrite=TRUE)
###############################################################################################################################################
# Below lines are not required, only for visuals
##################
# Plot GGM network
##################
library(graph)

gr = network.make.graph(alamar.net, node.labels, drop.singles = TRUE)
num.nodes(gr)

#' Number of directed ("forward") and undirected ("none") edges:
table(  edge.info(gr)$dir )
#' Nodes connected with many edges:
sort(node.degree(gr), decreasing=TRUE)[1:10]
library("Rgraphviz")
#' For a more beautiful plot of the network set node and edge parameters:

#' Set global node and edge attributes:
globalAttrs = list()
globalAttrs$edge = list(color = "black", lty = "solid", lwd = 1, arrowsize=1)
globalAttrs$node = list(fillcolor = gray(.95), shape = "ellipse", fixedsize = FALSE)

#' Set attributes of some particular nodes:
nodeAttrs = list()
nodeAttrs$fillcolor = c('570' = "red", "81" = "red") # highlight hub nodes

#' Set edge attributes:
edi = edge.info(gr) # edge directions and correlations
edgeAttrs = list()
edgeAttrs$dir =  edi$dir # set edge directions 
cutoff = quantile(abs(edi$weight), c(0.2, 0.8)) # thresholds for line width / coloring
edgeAttrs$lty = ifelse(edi$weight < 0, "dotted", "solid") # negative correlation
edgeAttrs$color = ifelse( abs(edi$weight <= cutoff[1]), "grey", "black") # lower 20% quantile
edgeAttrs$lwd = ifelse(abs(edi$weight >= cutoff[2]), 2, 1) # upper 20% quantile

#+ fig.width=8, fig.height=7
plot(gr, attrs = globalAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, "fdp")

#table(nmr_olink_cor$NMR == nmr_olink_cor_pvalue$NMR)




