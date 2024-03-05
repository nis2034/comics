# Process UK BioBank data
# clear the workspace
rm(list=ls())

library(tidyverse)
library(dplyr)
library(RNOmni)
library(GeneNet)

#######
#setwd("/home/thabib/Documents/UKB/")
setwd("/Users/tah4012/Documents/Research/UKB/")
###############
# Read NMR data and metadata
###############

nmr_dt <- read.table(file="Data/nmr_results.tsv", header=T, sep="\t", check.names = F, row.names = 1)
nmr_metadata <- read.table(file="Data/NMR_Pheno.txt", header=T, sep="\t", check.names = F)
col.order <- nmr_metadata$FieldID
nmr_dt <- subset(nmr_dt, select = nmr_metadata$FieldID)
nmr_dt <- nmr_dt[ , col.order]
nmr_metadata$Description <- gsub(" ", "_", nmr_metadata$Description)
colnames(nmr_dt) <- nmr_metadata$Description
nmr <- nmr_dt

head(nmr_dt[1:5,1:5])
# NA removed
nmr_df <- as.matrix(na.omit(nmr_dt))
###########################
# Graphical Gaussian Models (GeneNet)
# ggm.estimate.pcor
# gene names or metabolites in this case should be in columns and participants as row. 
# Typical gene expression matrix has gene names in rows and samples in columns. If the data is
# in this format then it needs to be transposed before applying GGM
###########################
# Compute Partial Correlations and Select Relevant Edges

pcor = ggm.estimate.pcor(nmr_df, method="dynamic") 
nmr.edges = network.test.edges(pcor, direct = TRUE)

#use strongest 200 edges
#nmr.net = extract.network(nmr.edges, method.ggm = "number", cutoff.ggm = 200)
nmr.net = extract.network(nmr.edges)
node.labels = as.character(colnames(nmr_df))

matching_indices_1 <- match(nmr.net$node1, row_number(node.labels))
matching_indices_2 <- match(nmr.net$node2, row_number(node.labels))

nmr.net$nodelabel1 <- node.labels[matching_indices_1]
nmr.net$nodelabel2 <- node.labels[matching_indices_2]

nmr.net$PLAT <- "NMR"
colnames(nmr.net) <- c("COR", "TRAIT1","TRAIT2","PVALUE","QVAL","PROB","LOG.SPVAR",
                       "PVAL.DIR","QVAL.DIR","PROB.DIR","DIRECTION",
                       "TRAITID1","TRAITID2", "PLAT")
write.table(nmr.net, file="NMR/NMR_GGM_Output.txt", sep="\t", col.names = NA)
save(nmr.net, file="NMR_GGM.rda")


################################## 
# OLINK
################################## 
rm(list=ls())
olink <- read.table(file="OLINK/results_olink_instance_0.tsv", header=T, sep="\t", check.names = F, row.names = 1)
colnames(olink) <- toupper(gsub("olink_instance_0.","",colnames(olink)))

# NA removed
olink_df <- as.matrix(na.omit(olink))
###########################
# Graphical Gaussian Models
# ggm.estimate.pcor
# gene names or metabolites in this case should be in columns and participants as row. Typical
# Typical gene expression matrix has gene names in rows and samples in columns. If the data is
# in this format then it needs to be transposed before applying GGM
###########################
# Compute Partial Correlations and Select Relevant Edges

pcor = ggm.estimate.pcor(olink_df, method="dynamic") 
olink.edges = network.test.edges(pcor, direct = TRUE)

#use strongest 200 edges
#olink.net = extract.network(nmr.edges, method.ggm = "number", cutoff.ggm = 200)
olink.net = extract.network(olink.edges)
node.labels = as.character(colnames(olink_df))

matching_indices_1 <- match(olink.net$node1, row_number(node.labels))
matching_indices_2 <- match(olink.net$node2, row_number(node.labels))

olink.net$nodelabel1 <- node.labels[matching_indices_1]
olink.net$nodelabel2 <- node.labels[matching_indices_2]

olink.net$PLAT <- "OLINK"
colnames(olink.net) <- c("COR", "TRAIT1","TRAIT2","PVALUE","QVAL","PROB","LOG.SPVAR",
                         "PVAL.DIR","QVAL.DIR","PROB.DIR","DIRECTION",
                         "TRAITID1","TRAITID2", "PLAT")
write.table(olink.net, file="OLINK_GGM_Output.txt", sep="\t", col.names = NA)
save(olink.net, file="OLINK_GGM.rda")

###############################
# CLINICAL
###############################
rm(list=ls())
clin <- read.table(file="CLINICAL/clinical_results.tsv", header=T, sep="\t", check.names = F, row.names = 1)
head(clin)
clin_info <- read.table(file="CLINICAL/clinical.txt", header=T, sep="\t", check.names = F)
col.order <- clin_info$FieldID
clin <- subset(clin, select = clin_info$FieldID)
clin <- clin[ , col.order]
clin_info$Description <- gsub(" ", "_", clin_info$Description)
colnames(clin) <- clin_info$Description
#Remove NA
clin_df <- as.matrix(na.omit(clin))

###########################
# Graphical Gaussian Models
# ggm.estimate.pcor
# gene names or metabolites in this case should be in columns and participants as row. Typical
# Typical gene expression matrix has gene names in rows and samples in columns. If the data is
# in this format then it needs to be transposed before applying GGM
###########################
# Compute Partial Correlations and Select Relevant Edges

pcor = ggm.estimate.pcor(clin_df, method="dynamic") 
clin.edges = network.test.edges(pcor, direct = TRUE)

#use strongest 200 edges
#olink.net = extract.network(nmr.edges, method.ggm = "number", cutoff.ggm = 200)
clin.net = extract.network(clin.edges)
node.labels = as.character(colnames(clin_df))

matching_indices_1 <- match(clin.net$node1, row_number(node.labels))
matching_indices_2 <- match(clin.net$node2, row_number(node.labels))

clin.net$nodelabel1 <- node.labels[matching_indices_1]
clin.net$nodelabel2 <- node.labels[matching_indices_2]

clin.net$PLAT <- "CLIN"
colnames(clin.net) <- c("COR", "TRAIT1","TRAIT2","PVALUE","QVAL","PROB","LOG.SPVAR",
                        "PVAL.DIR","QVAL.DIR","PROB.DIR","DIRECTION",
                        "TRAITID1","TRAITID2", "PLAT")
write.table(clin.net, file="CLIN_GGM_Output.txt", sep="\t", col.names = NA)

save(clin.net, file="CLIN_GGM.rda")

###############################################################################################################################################
# Below lines are not required, only for visuals
##################
# Plot GGM network
##################
library(graph)

gr = network.make.graph(nmr.net, node.labels, drop.singles = TRUE)
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


##################
# Plot GGM network
##################
library(graph)
gr = network.make.graph(olink.net, node.labels, drop.singles = TRUE)
num.nodes(gr)
edge.info(gr)$weight

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
nodeAttrs$fillcolor = c('IL19' = "red", "81" = "red") # highlight hub nodes

#' Set edge attributes:
edi = edge.info(gr) # edge directions and correlations
edgeAttrs = list()
edgeAttrs$dir =  edi$dir # set edge directions 
cutoff = quantile(abs(edi$weight), c(0.2, 0.8)) # thresholds for line width / coloring
edgeAttrs$lty = ifelse(edi$weight < 0, "dotted", "solid") # negative correlation
edgeAttrs$color = ifelse( abs(edi$weight <= cutoff[1]), "grey", "black") # lower 20% quantile
edgeAttrs$lwd = ifelse(abs(edi$weight >= cutoff[2]), 2, 1) # upper 20% quantile

#+ fig.width=8, fig.height=7
#plot(gr, attrs = globalAttrs, nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs, "fdp")
plot(gr, attrs = globalAttrs, nodeAttrs = makeNodeAttrs(gr, fontsize=30), 
     edgeAttrs = edgeAttrs, "fdp")

edges.results <- ggm.list.edges(pcor)
nrow(edges.results)

edges.results[50:70,]
idx1 <- match(edges.results$node1, row_number(node.labels))
idx2 <- match(edges.results$node2, row_number(node.labels))

edges.results$node1 <- node.labels[idx1]
edges.results$node2 <- node.labels[idx2]
network.make.dot(filename="olink_ggm.dot", edges.results, node.labels = node.labels, main = "OLINK graph")


##################
# Plot GGM network
##################
library(graph)
gr = network.make.graph(clin.net, node.labels, drop.singles = TRUE)
num.nodes(gr)
edge.info(gr)$weight

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
nodeAttrs$fillcolor = c('IL19' = "red", "81" = "red") # highlight hub nodes

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
plot(gr, attrs = globalAttrs, nodeAttrs = makeNodeAttrs(gr, fontsize=30), 
     edgeAttrs = edgeAttrs, "fdp")


