<<<<<<< HEAD
#' Class NJtree 
#' 
#' @param patientID patient name
#' @return NJtree object
#' 
#' @examples
#' maf <- read.Maf(patientID)
#' njtree <- read.NJtree(maf)
#' getNJtreenj(njtree)
#' getNJtreemut_sort(njtree)
#' getNJtreeSignature(njtree)


# import pkgs
library(plyr)
library(tidyr)
library(ape)
library(ggplot2)
library(dplyr)
# dependencies of deconstructSigs
library(reshape2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomeInfoDb)
library(grDevices)
library(graphics)
library(utils)
library(deconstructSigs)
# data frame needed
library(plyr)



# set NJtree object
read.NJtree <- function(maf, use.indel = FALSE, use.ccf = FALSE, mut.signature = TRUE, 
                        sig.min.mut.number = 50, ccf.mutation.id, ccf.mutation.sep){
  maf.dat <- maf@data
  ccf <- maf@ccf.loci
  patientID <- maf@patientID
  if(use.ccf){
    mut_sort <- mut_ccf_sort(maf.dat = maf.dat, ccf = ccf, use.indel, ccf.mutation.id, ccf.mutation.sep)
    if(is.null(ccf)){
      stop("Missing ccf file. Check whether maf@ccf.loci is NULL")
    }
  }
  else{
    mut_sort <- mut_binary_sort(maf.dat = maf.dat, use.indel = use.indel)
  }
  mat.nj = nj(dist.gene(t(mut_sort)))
  branch <- read.njtree(mat.nj)
  list.sig <- Mutational_sigs_tree(maf.dat = maf.dat, branch, patientID, mut.threshold = sig.min.mut.number)
  signature <- list.sig[[1]]
  mut.branches <- list.sig[[2]]
  njtree <- new('NJtree', nj = mat.nj, mut_sort = mut_sort, signature = signature, mut_branches = mut.branches)
  return(njtree)
}
#Prevent class 'phylo' from not existing
setClass('phylo')
#Class NJtree
setClass('NJtree', slots = c(nj = 'phylo', mut_sort = 'matrix', signature = 'data.frame', mut_branches = 'list'))
#extract mut_sort from NJtree object
setGeneric("getMutSort", function(x, signature){standardGeneric("getMutSort")})
setMethod("getMutSort", 'NJtree', function(x){x@mut_sort})
#extract nj from NJtree object
setGeneric("getPhyloTree", function(x, signature){standardGeneric("getPhyloTree")})
setMethod("getPhyloTree",'NJtree', function(x){x@nj})
#extract siganture from NJtree object
setGeneric("getNJtreeSignature", function(x, signature){standardGeneric("getNJtreeSignature")})
setMethod("getNJtreeSignature", 'NJtree', function(x){x@signature})


###### output test ######
patientID = "311252"
dat.dir = "./data/multi_lesion"
# NJtree object
maf <- read.Maf(patientID, dat.dir)
njtree <- read.NJtree(maf, use.indel = T, use.ccf = F)
getMutSort(njtree)
getPhyloTree(njtree)
getNJtreeSignature(njtree)

=======
#' Class NJtree 
#' 
#' @param patientID patient name
#' @return NJtree object
#' 
#' @examples
#' maf <- read.Maf(patientID)
#' njtree <- read.NJtree(maf)
#' getNJtreenj(njtree)
#' getNJtreemut_sort(njtree)
#' getNJtreeSignature(njtree)


# import pkgs
library(plyr)
library(tidyr)
library(ape)
library(ggplot2)
library(dplyr)
# dependencies of deconstructSigs
library(reshape2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomeInfoDb)
library(grDevices)
library(graphics)
library(utils)
library(deconstructSigs)
# data frame needed
library(plyr)



# set NJtree object
read.NJtree <- function(maf, use.indel = FALSE, use.ccf = FALSE, mut.signature = TRUE, 
                        sig.min.mut.number = 50, ccf.mutation.id, ccf.mutation.sep){
  maf.dat <- maf@data
  ccf <- maf@ccf.loci
  if(use.ccf){
    mut_sort <- mut_ccf_sort(maf.dat = maf.dat, ccf = ccf, use.indel, ccf.mutation.id, ccf.mutation.sep)
    if(is.null(ccf)){
      stop("Missing ccf file. Check whether maf@ccf.loci is NULL")
    }
  }
  else{
    mut_sort <- mut_binary_sort(maf.dat = maf.dat, use.indel = use.indel)
  }
  mat.nj = nj(dist.gene(t(mut_sort)))
  branch <- read.njtree(mat.nj)
  list.tree <- Mutational_sigs_tree(maf.dat = maf.dat, branch, mut.threshold = sig.min.mut.number)
  signature <- list.tree[[1]] 
  mut_branches <- list.tree[[2]]
  njtree <- new('NJtree', nj = mat.nj, mut_sort = mut_sort, signature = signature, mut_branches = mut_branches)
  return(njtree)
}
#Prevent class 'phylo' from not existing
setClass('phylo')
#Class NJtree
setClass('NJtree', slots = c(nj = 'phylo', mut_sort = 'matrix', signature = 'data.frame', mut_branches = 'data.frame'))
#extract mut_sort from NJtree object
setGeneric("getMutSort", function(x, signature){standardGeneric("getMutSort")})
setMethod("getMutSort", 'NJtree', function(x){x@mut_sort})
#extract nj from NJtree object
setGeneric("getPhyloTree", function(x, signature){standardGeneric("getPhyloTree")})
setMethod("getPhyloTree",'NJtree', function(x){x@nj})
#extract siganture from NJtree object
setGeneric("getNJtreeSignature", function(x, signature){standardGeneric("getNJtreeSignature")})
setMethod("getNJtreeSignature", 'NJtree', function(x){x@signature})


###### output test ######
patientID = "311252"
dat.dir = "./data/multi_lesion"
# NJtree object
maf <- read.Maf(patientID, dat.dir)
njtree <- read.NJtree(maf, use.indel = T, use.ccf = F)
getMutSort(njtree)
getPhyloTree(njtree)
getNJtreeSignature(njtree)

>>>>>>> master
