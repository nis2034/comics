rm(list=ls())
library("readxl")
library("tidyverse")
library("visNetwork")
library("shiny")
library("shinydashboard")
library("DT")
library("dplyr")

# Read Excel file and displays its sheets
## list all excel sheets 
sheets <- excel_sheets(path = "NETWORKS.xlsx")
#sheets <- excel_sheets(path = "NETWORKS_tanwir.xlsx")
sheets_ok <- c()
sheets_notok <- c()
sheets_flag <- c()
# loop through all the sheets for four columns
must_have_columns <- c("TRAITID1","TRAITID2","PVALUE","BETA", "COR")
annotation_columns <- c("TRAITID","SHORTNAME","PLAT")
# i=1
n=0
m=0
# or pre-allocate for slightly more efficiency
datalist = vector("list", length = length(sheets))

for (i in 1:length(sheets)) {
  #sheetContent = read_excel(path = "NETWORKS_tanwir.xlsx", sheet = sheets[i], .name_repair = ~make.unique(.x, sep = "_"))
  sheetContent = read_excel(path = "NETWORKS.xlsx", sheet = sheets[i], .name_repair = ~make.unique(.x, sep = "_"))
  musthave_cols <- unlist(lapply(must_have_columns, function(header) {
    grep(paste0("^",header,"$"), names(sheetContent), ignore.case = TRUE)
  }))
  annotation_cols <- unlist(lapply(annotation_columns, function(header) {
    grep(paste("^",header,"$",sep = ""), names(sheetContent), ignore.case = TRUE)
  }))
  if(length(musthave_cols) == 4){
    m=m+1
    idx <- grep("^COR$|^BETA$", names(sheetContent))
    dat <- data.frame(sheetContent %>% select(musthave_cols, idx[1]))
    
    
    dat$type <- rep(sheets[i])
    colnames(dat) <- c("to","from","pvalue","weight","type")
    dat$i <- i  # maybe you want to keep track of which iteration produced it?
    dat$id <- paste0(sheets[i],"_0",rownames(dat))
    datalist[[i]] <- dat # add it to your list
    sheets_ok[m] <- sheets[i]
  } else if(length(annotation_cols) == 3) {
    idx <- grep("^PLAT$", names(sheetContent))
    anno <- sheetContent %>% select(all_of(annotation_cols))
    
  } else {
    n=n+1
    columns_notfound <- toString(must_have_columns[!(must_have_columns %in% names(sheetContent))])
    sheets_flag[n] <- paste(sheets[i],": missing columns - ", columns_notfound)
    sheets_notok[n] <- sheets[i]
  }
}
datalist %>% head()
length(datalist)
class(datalist)

datalist = datalist[!(sapply(datalist, is.null))]

my_list <- list(1, NULL, 3, NULL, 5, NULL)

filtered_list <- my_list[!sapply(my_list, is.null)]

print(filtered_list)

all_edges = do.call(rbind, datalist)
all_edges$sign = sign(all_edges$weight)
all_edges$weight = abs(all_edges$weight)

all_nodes = data.frame(TRAITID = unique(sort(c(all_edges$from, all_edges$to))))

head(all_nodes)
tail(all_edges)
# get node annotations
all_nodes = left_join(all_nodes, anno)
#all_nodes$SHORTNAME <- paste0(all_nodes$PLAT,":",all_nodes$TRAITNAME)

# replace all_edges$to and $from ids with all_nodes$SHORTNAME
TRAITID2SHORTNAME = all_nodes$SHORTNAME
names(TRAITID2SHORTNAME) = all_nodes$TRAITID

SHORTNAME2TRAITID = all_nodes$TRAITID
names(SHORTNAME2TRAITID) = all_nodes$SHORTNAME

all_edges$from = TRAITID2SHORTNAME[all_edges$from]
all_edges$to = TRAITID2SHORTNAME[all_edges$to]
all_nodes$id = all_nodes$SHORTNAME
names(all_nodes)[which(names(all_nodes) == "PLAT")] = "plat"

# round the p-values and weights
all_edges$pvalue = signif(all_edges$pvalue, digits = 2)
all_edges$weight = signif(all_edges$weight, digits = 2)


######################################################################################################################
###############################################################################################################################################################################
# a data structure for a network
#########################################################

fullnet = list(
  edges = all_edges,
  nodes = all_nodes
)
head(all_nodes)
##################################################################
# function: neighbors - extract all neighbors of a given node list
# input: a node list
# output: a node list
##################################################################

##########################################################################
##########################################################################
# 
# # reduce the network to exclude STAT
# traitnetnodes = fullnet %>% network2nodes() 
# traitnetnodes = traitnetnodes[grep("STAT: ", traitnetnodes, invert = TRUE)]
# traitnet = nodes2network(traitnetnodes, fullnet) 
# 
# # get the STATnodes
# STATnetnodes = fullnet %>% network2nodes() 
# STATnetnodes = STATnetnodes[grep("STAT: ", STATnetnodes, invert = FALSE)]

##########################################################################
##########################################################################
# "cg19693031:chr1:144152909 TXNIP" is TXNIP
# 3485-28_2 is B2M
##########################################################################
if (FALSE) { # test code
  
  summary_network(fullnet)
  
  TXNIP = fullnet$nodes$id[grep("cg19693031", fullnet$nodes$id)]
  neighbors(TXNIP, fullnet)
  
  LEPR = fullnet$nodes$id[grep("Leptin receptor", fullnet$nodes$id)]
  neighbors(LEPR, fullnet) %>% maxneighbors_noSTAT(fullnet, limit = 20)
  neighbors(LEPR, fullnet) %>% maxneighbors_noSTAT(fullnet, limit = 200)
  
  neighbors("STAT: AGE", fullnet)
  neighbors("SOMA: C9 : Complement component C9", fullnet)
  neighbors("some error", fullnet)
  neighbors(c(), fullnet)
  
  maxneighbors(c("1","2"), fullnet)
  maxneighbors("cg19693031", fullnet, limit = 100)
  
  neighbors(c(TXNIP), fullnet) 
  neighbors(c(TXNIP), fullnet) %>% neighbors(fullnet) 
  neighbors(c(TXNIP), fullnet) %>% neighbors(fullnet) %>% neighbors(fullnet) 
  
  subnet = maxneighbors(c(TXNIP), fullnet, limit = 100) %>%  nodes2network(fullnet)
  summary_network(subnet)
  
  subnet = maxneighbors_noSTAT(c(TXNIP), fullnet, limit = 100) %>%  nodes2network(fullnet)
  summary_network(subnet)
  
  neighbors(TXNIP, subnet)
  neighbors(TXNIP, fullnet)
  
}
##########################################################################
##########################################################################

# add colors for platforms to network
# to be used as: PLATcols[fullnet$nodes$plat]
# 

PLATlist = unique(fullnet$nodes$plat)

PLATcols = rep("#999999",length(PLATlist))
names(PLATcols) = PLATlist


# num_plats <- length(PLATlist)
# colors <- rainbow(num_plats)
# PLATcols = colors
# names(PLATcols) = PLATlist



PLATshapes = rep("star",length(PLATlist))
# PLATshapes = rep("icon",length(PLATlist))
names(PLATshapes) = PLATlist

# # define PLAT colors manually
PLATcols["DNA"] = "#23bbee"
PLATcols["SOMA"] = "#a62281"
PLATcols["BRAIN"] = "#f2921f"
PLATcols["BM"] = "#ffc815"
PLATcols["LD"] = "#ffc815"
PLATcols["CPG"] = "#145da9"
PLATcols["CLIN"] = "#a0b6a8"
PLATcols["CM"] = "#57ba47"
PLATcols["IgA"] = "#e41d30"
PLATcols["RNA"] = "#5c2d83"
PLATcols["HD4"] = "#57ba47"
PLATcols["IgG"] = "#e41d30"
PLATcols["miRNA"] = "#5c2d83"
PLATcols["OLINK"] = "#a62281"
PLATcols["PGP"] = "#e41d30"
PLATcols["PM"] = "#57ba47"
PLATcols["SM"] = "#57ba47"
PLATcols["UM"] = "#57ba47"
PLATcols["STAT"] = "#EEEEEE"
PLATcols["GWAS"] = "#EEEEEE"

# define PLAT colors manually
PLATshapes["UM"] = "square"
PLATshapes["CM"] = "square"
PLATshapes["SM"] = "triangle"
PLATshapes["DNA"] = "diamond"
PLATshapes["RNA"] = "diamond"
PLATshapes["CPG"] = "diamond"

PLATshapes["STAT"] = "circle"
PLATshapes["GWAS"] = "dot"

fullnet$nodes$color = PLATcols[fullnet$nodes$plat]

fullnet$nodes$group = fullnet$nodes$plat

fullnet$nodes$shape = PLATshapes[fullnet$nodes$plat]

fullnet$edges$title = paste(fullnet$edges$type, ": ", fullnet$edges$id,
                            ", p=", fullnet$edges$pvalue,
                            ", beta=",
                            ifelse(fullnet$edges$sign >0, "", "-"), fullnet$edges$weight,
                            sep = "")
fullnet$edges$color = ifelse(fullnet$edges$sign > 0, "blue", "red")

#convert all pvalues == 0 to a very small number so -log10 won't throw an error
fullnet$edges$pvalue[fullnet$edges$pvalue == 0] <- 1E-100

edge_thickness = -log10(fullnet$edges$pvalue)
minVal = min(edge_thickness)-1
maxVal = max(edge_thickness)-1
numItv = 10
seq_range <- round(seq(minVal, maxVal, length.out=numItv))

fullnet$edges$width = cut(edge_thickness, seq_range, labels = F)

# list of ids to choose from
plat_list = c("ALL", fullnet$nodes$plat %>% unique %>% sort())

# limit the number of nodes here
max_nodes_list = c(1,20,40,60,80,100,150,200)