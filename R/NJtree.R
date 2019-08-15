#' Class NJtree 
#' 
#' @param maf return from readMaf()
#' @param use.indel Seclet SNP in Variant type
#' @param sig.min.mut.number minimum mutation number in each branch
#' @param ccf.mutation.id determine the id of mutation
#' @param ccf.mutation.sep sep of id
#' @return NJtree object
#' 
#' @examples
#' maf <- readMaf(patientID = "311252",mafDir = './inst/exdata/multi_lesion/maf/311252.maf', refBuild = "hg19")
#' njtree <- read.NJtree(maf)
#' getNJtreenj(njtree)
#' getNJtreemut_sort(njtree)
#' getNJtreeSignature(njtree)


# set NJtree object
NJtree <- function(maf, use.indel = FALSE, use.ccf = FALSE, mut.signature = TRUE, 
                   sig.min.mut.number = 50, ccf.mutation.id, ccf.mutation.sep){
  maf.dat <- maf@data
  ccf <- maf@ccf.loci
  patientID <- maf@patientID
  if(use.ccf){
    mut_sort <- mut_ccf_sort(maf.dat = maf.dat, ccf = ccf, use.indel, ccf.mutation.id, ccf.mutation.sep)
    if(is.null(ccf)){
      stop("Missing ccf file. Check whether maf@ccf.loci is NULL")
    }
  }else{
    mut_sort <- mut_binary_sort(maf.dat = maf.dat, use.indel = use.indel)
  }
  mat.nj = nj(dist.gene(t(mut_sort)))
  branch <- read.njtree(mat.nj)
  mut.branches <- treeMutationalBranches(maf, branch)
  njtree <- new('NJtree', nj = mat.nj, mut_sort = mut_sort, patientID = patientID, mut_branches = mut.branches)
  return(njtree)
}
#Prevent class 'phylo' from not existing
setClass('phylo')
#Class NJtree
setClass('NJtree', slots = c(nj = 'phylo', mut_sort = 'matrix', patientID = 'character', mut_branches = 'list'))
#extract mut_sort from NJtree object
setGeneric("getMutSort", function(x, signature){standardGeneric("getMutSort")})
setMethod("getMutSort", 'NJtree', function(x){x@mut_sort})
#extract nj from NJtree object
setGeneric("getPhyloTree", function(x, signature){standardGeneric("getPhyloTree")})
setMethod("getPhyloTree",'NJtree', function(x){x@nj})
#extract siganture from NJtree object
setGeneric("getNJtreeSignature", function(x, signature){standardGeneric("getNJtreeSignature")})
setMethod("getNJtreeSignature", 'NJtree', function(x){x@signature})

