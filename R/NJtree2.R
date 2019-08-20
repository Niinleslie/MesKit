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
#' maf.File <- system.file("extdata/multi_lesion/maf", "311252.maf", package = "Meskit")
#' sampleInfo.File <- system.file("extdata/multi_lesion", "sample_info.txt", package = "Meskit")
#' pyCloneCluster <- system.file("extdata/multi_lesion/ccf", "311252.cluster.tsv", package = "Meskit")
#' pyCloneLoci <- system.file("extdata/multi_lesion/ccf", "311252.loci.tsv", package = "Meskit")
#' maf <- readMaf(patientID = "311252", mafFile = maf.File, sampleInfo = sampleInfo.File, refBuild = "hg19")
#' njtree <- NJtree(maf)
#' getNJtreenj(njtree)
#' getNJtreemut_sort(njtree)
#' getNJtreeSignature(njtree)

# set NJtree object
NJtree <- function(maf, use.indel = FALSE, mut.signature = TRUE, 
                   sig.min.mut.number = 50, ccf.mutation.id, ccf.mutation.sep){
  maf.dat <- maf@data
  patientID <- maf@patientID
  refBuild <- paste("BSgenome.Hsapiens.UCSC.", maf@ref.build, sep = "")

  mut_sort <- mut_binary_sort(maf.dat = maf.dat, use.indel = use.indel)
  mat.nj = nj(dist.gene(t(mut_sort)))
  branch <- read.njtree(mat.nj)
  mut.branches <- treeMutationalBranches(maf, branch)
  ccf_sort <- matrix()
  if(is.null(maf@ccf.loci)){
    ccf_sort <- mut_ccf_sort(maf.dat = maf.dat, ccf = ccf, use.indel, ccf.mutation.id, ccf.mutation.sep)
  }

  njtree <- new('NJtree', patientID = patientID, nj = mat.nj, mut_sort = mut_sort, ccf_sort = ccf_sort, mut_branches = mut.branches, refBuild = maf@ref.build)
  return(njtree)
}
#Prevent class 'phylo' from not existing
setClass('phylo')
#Class NJtree
setClass('NJtree', slots = c(nj = 'phylo', patientID = 'character', mut_sort = 'matrix', ccf_sort = 'matrix', mut_branches = 'list', refBuild = 'character'))
#extract mut_sort from NJtree object
setGeneric("getMutSort", function(x, signature){standardGeneric("getMutSort")})
setMethod("getMutSort", 'NJtree', function(x){x@mut_sort})
#extract nj from NJtree object
setGeneric("getPhyloTree", function(x, signature){standardGeneric("getPhyloTree")})
setMethod("getPhyloTree",'NJtree', function(x){x@nj})
