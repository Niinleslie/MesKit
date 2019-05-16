#' Prepare inputs for SCHISM
#' @description read tsv documents from PyClone and figure out the list of filtered clusters. According to those 
#' targeted clusters, filter targeted mutations and reorganize data as two outputs: clusterEstimates.tsv and 
#' mutation_to_cluster.tsv, which are the input of SCHISM.
#' 
#' @param dir.cluster.tsv specify the directory of cluster.tsv document
#' @param dir.loci.tsv specify the directory of loci.tsv document
#' @param dir.output specify the directory of two output for SCHISM
#' @return clusterEstimates.tsv, mutation_to_cluster.tsv
#' 
#' @examples
#' \dontrun{
#' prepareSchismInput(dir.cluster.tsv, dir.loci.tsv, dir.output)
#'}

## prepare Schism input
prepareSchismInput <- function(dir.cluster.tsv, dir.loci.tsv, dir.output){
  # read mutations in cluster.tsv from PyClone and get targeted clusters
  cluster.tsv = read.table(dir.cluster.tsv, sep="\t", stringsAsFactors=F, header=T)
  cluster_filter = cluster.tsv[!is.na(cluster.tsv$cluster_id) & cluster.tsv$size >= 5 & cluster.tsv$mean >= 0.1,]
  cluster_ls  = unique(cluster_filter$cluster_id)
  
  # read mutations within targeted clusters in loc.tsv file from PyClone 
  loci.tsv = read.table(dir.loci.tsv, sep="\t", stringsAsFactors=F, header=T)
  loci_filtered <- loci.tsv[which(loci.tsv$cluster_id %in% cluster_ls & loci.tsv$sample_id %in% sep3),]
  loci_filtered <- na.omit(loci_filtered)
  
  # generate two outputs for SCHISM: clusterEstimates.tsv, mutation_to_cluster.tsv
  clusterEstimates.tsv <- data.frame(loci_filtered$sample_id, loci_filtered$mutation_id, 
                                     loci_filtered$cellular_prevalence, loci_filtered$cellular_prevalence_std)
  names(clusterEstimates.tsv) = c("sampleID",	"mutationID",	"cellularity", "sd")
  
  mutation_to_cluster.tsv <- data.frame(loci_filtered$mutation_id, loci_filtered$cluster_id)
  names(mutation_to_cluster.tsv) = c("mutationID",	"clusterID")
  
  # set the directory of outputs and make two SCHISM-needed files 
  setwd(dir.output)
  write.table(clusterEstimates.tsv, file="W.clusterEstimates.tsv", sep = "\t", quote=F, row.names=F)
  write.table(mutation_to_cluster.tsv, file="W.mutation-to-cluster.tsv", sep = "\t", quote=F, row.names=F)
}

#' Get outputs from SCHISM and draw the fishplot
#' @description read cluster.cellularity and GA.consensusTree documents from SCHISM. Construct fraction table and
#' parents vector for drawing fishplot.
#' 
#' @param dir.cluster.cellularity specify the directory of cluster.cellularity document
#' @param dir.GA.consensusTree specify the directory of GA.consensusTree document
#' @return the fishplot showing the evolution relationship predicted by SCHISM
#'
#' @examples
#' \dontrun{
#' schism2Fishplot(dir.cluster.cellularity, dir.GA.consensusTree)
#'}

## dependencies of fishplot
library(png)
library(Hmisc)
library(plotrix)
library(fishplot)

## generate results for createFishPlotObjects
schism2Fishplot <- function(dir.cluster.cellularity, dir.GA.consensusTree){
  # get cellularity infomation
  cluster.cellularity = read.table(dir.cluster.cellularity, sep="\t", stringsAsFactors=F, header=T)
  
  # get the list of samples and clusters in this function
  ls.sample = unique(cluster.cellularity$sampleID)
  ls.cluster = unique(cluster.cellularity$clusterID)
  
  # build the cancer celluarity/fraction table for fishplot
  frac.c = c()
  for (sample_name in ls.sample) {
    sample <- cluster.cellularity[which(cluster.cellularity$sampleID == sample_name), ]$cellularity
    frac.c <- c(frac.c,sample)
  }
  frac.table = matrix(frac.c*100, ncol = length(ls.sample))
  rownames(frac.table) <- ls.cluster
  colnames(frac.table) <- ls.sample
  
  # read the evolution relationship of different subclones
  GA.consensusTree = read.table(dir.GA.consensusTree, sep="\t", stringsAsFactors=F, header=T)
  
  # figure out clusters which will be first nodes of the evolution tree
  ls.ends = unique(GA.consensusTree[which(!GA.consensusTree$parent %in% GA.consensusTree$child), 1])
  
  # make sure the first node of the tree would be first in parents
  if (!rownames(frac.table)[1] %in% ls.ends){
    frac.end = frac.table[which(rownames(frac.table) == ls.ends[1]),]
    frac.table <- frac.table[-which(rownames(frac.table) == ls.ends[1]),]
    frac.table <- rbind2(frac.end, frac.table)
    rownames(frac.table)[1] <- ls.ends[1]
  }
  
  # rearrange subclonal evolution relationship as a vector(parents)
  parents <- c(1:length(ls.cluster))
  for (cluster_name in ls.cluster){
    parent_cluster <- GA.consensusTree[which(GA.consensusTree$child == cluster_name), 1]
    if (cluster_name %in% ls.ends){
      parents[which(rownames(frac.table) == cluster_name)] <- 0
    } else{
      parents[which(rownames(frac.table) == cluster_name)] <- which(rownames(frac.table) == min(parent_cluster))
    }
  }
  
  # fishplot printing
  fish = createFishObject(frac.table, parents, timepoints = c(1:length(ls.sample)), fix.missing.clones=TRUE)
  if (length(ls.cluster) > 10) {
    fish = setCol(fish, as.character(ls.cluster))
  }
  fish = layoutClones(fish)
  fishPlot(fish, shape="spline", title.btm="PatientID", cex.title=0.5,
           vlines=seq(1, length(ls.sample)), vlab=ls.sample, pad.left=0.5)
}


# # directorys
# dir.cluster.tsv = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/hu.cluster.tsv"
# dir.loci.tsv = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/hu.loci.tsv"
# dir.output = "/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/Schism"
# 
# dir.GA.consensusTree = "/home/ninomoriaty/R_Project/EvolCancer/E2/sample-output/E2.GA.consensusTree"
# dir.cluster.cellularity = "/home/ninomoriaty/R_Project/EvolCancer/E2/sample-output/E2.cluster.cellularity"


# # test: separate the first 3 sample and the next 4 sample to genrate a more proper result
# sep3 <- c("WGC033386D", "WGC033387D", "WGC033388D")
# sep4 <- c("WGC033389D", "WGC033391D", "WGC033392D", "WGC033393D")
# & loci.tsv$sample_id %in% sep4
# #test: calculate common mutation
# mutation_ls <- loci_filtered[ ,1]
# mut.common <- c()
# for (mutation in mutation_ls){
#   if (length(unique(loci_filtered[which(loci_filtered$mutation_id == mutation),2])) == 4){
#     mut.common <- c(mut.common, mutation)
#   }
# }
# length(unique(mut.common))

# # test: what if the value is ok 
# frac.table[3,1] <- 30
# frac.table[2,4] <- 10
# frac.table[3,5] <- 49.60
# frac.table[2,5] <- 49.50
# frac.table[2,6] <- 49.70
# 
# # test: first 3 columns errors(couldn't make a combination of subclones)
# frac.table[3,1] <- 38.75 # minus 1
# 
# # test: last 4 columns errors
# frac.table[4,2] <- 49.72
# frac.table[3,4] <- 43.00
# 
# # tips: E2 test
# frac.table[which(frac.table != 0)] <- frac.table[which(frac.table != 0)] -1
# frac.table[1, ] <- 100 

### Method1: get the result from HT, but however, it's so strange that I couldn't find the relationship by cost.
## you need to know the method resolve this data, however, it is given by Schism..
# HT.cpov = read.table("/home/ninomoriaty/R_Project/EvolCancer/EvolCancer/Schism/0_W_results.HT.cpov", sep="\t", stringsAsFactors=F, header=T)
# parents_filter <- HT.cpov[which(HT.cpov$topologyCost == 0 & HT.cpov$parent.cluster != HT.cpov$child.cluster), c("parent.cluster", "child.cluster")]
# parents_freq <- as.data.frame(table(parents_filter$parent.cluster))
# order(as.data.frame(table(parents_filter$parent.cluster)), decreasing = T)
# parents_filter[order(, decreasing = T), ]

# 
# df.parents <- data.frame(parent.cluster = c(NA), child.cluster = c(NA))
# for (cluster_name in ls.cluster){
#   child_cluster <- parents_filter[which(parents_filter$parent.cluster == cluster_name), ]
#   if (length(child_cluster[,1]) == 0){
#     next()
#   } else {
#     df.parents <- df.parents[-which(df.parents$child.cluster %in% child_cluster$child.cluster),]
#   }
#   df.parents <- rbind(df.parents, child_cluster)
# }

## get parents
# ls.ends <- df.parents[-which(df.parents$parent.cluster %in% df.parents$child.cluster), ]$parent.cluster
# parents <- c(1:length(ls.cluster))
# 
# for (cluster_name in ls.cluster){
#   if (cluster_name %in% ls.ends){
#     parents[which(rownames(frac.table) == cluster_name)] <- 0
#   } else{
#     parents[which(rownames(frac.table) == cluster_name)] <- which(rownames(frac.table) == df.parents[which(df.parents$child.cluster == cluster_name), 1])
#   }
# }

