#' Class NJtree
#' 
#' 
#' @param maf return from readMaf()
#' @param use.indel Seclet SNP in Variant type
#' @param ccf.mutation.id manually specify which columns could be joint by ccf.mutation.sep to get the same format of mutation id in ccfy
#' @param ccf.mutation.sep manually specify the separator character.Values on each line of the ccf.mutation.id are separated by this character. I
#' @return NJtree object
#' 
#' @exportClass NJtree
#' @export getNJtree
#' 
#' @examples
#' maf.File <- system.file("extdata/multi_lesion/maf", "311252.maf", package = "Meskit")
#' sampleInfo.File <- system.file("extdata/multi_lesion", "sample_info.txt", package = "Meskit")
#' pyCloneCluster <- system.file("extdata/multi_lesion/ccf", "311252.cluster.tsv", package = "Meskit")
#' pyCloneLoci <- system.file("extdata/multi_lesion/ccf", "311252.loci.tsv", package = "Meskit")
#' maf <- readMaf(patientID="311252", mafFile=maf.File, sampleInfoFile=sampleInfoFile, refBuild="hg19")
#' njtree <- getNJtree(maf)


# set NJtree object
getNJtree <- function(maf, use.indel = FALSE, 
                   ccf.mutation.id = c("Hugo_Symbol","Chromosome","Start_Position"),
                   ccf.mutation.sep = ":"){
  maf.dat <- maf@data
  patientID <- maf@patientID
  refBuild <- paste("BSgenome.Hsapiens.UCSC.", maf@ref.build, sep = "")
  mut_sort.id <- mut_binary_sort(maf.dat = maf.dat, use.indel = use.indel)
  mut_sort <- as.matrix(mut_sort.id[, -1])
  mat.nj = nj(dist.gene(t(mut_sort)))
  branchAlias <- read.njtree(mat.nj)
  mut_branches <- .treeMutationalBranches(maf, branchAlias, mut_sort.id)
  ccf_sort <- matrix()
  if(!is.null(maf@ccf.loci)){
    ccf_sort <- mut_ccf_sort(maf.dat = maf.dat, ccf = maf@ccf.loci, use.indel, ccf.mutation.id, ccf.mutation.sep)
  }

  njtree <- new('NJtree', patientID = patientID, nj = mat.nj, mut_sort = mut_sort, mut_branches = mut_branches, ccf_sort = ccf_sort, refBuild = maf@ref.build)
  return(njtree)
}
#Prevent class 'phylo' from not existing
setClass('phylo')
#Class NJtree
setClass('NJtree', slots = c(nj = 'phylo', patientID = 'character', mut_sort = 'matrix', mut_branches = 'list', ccf_sort = 'matrix', refBuild = 'character'))
#extract mut_sort from NJtree object
setGeneric("getMutSort", function(x, signature){standardGeneric("getMutSort")})
setMethod("getMutSort", 'NJtree', function(x){x@mut_sort})
#extract nj from NJtree object
setGeneric("getPhyloTree", function(x, signature){standardGeneric("getPhyloTree")})
setMethod("getPhyloTree",'NJtree', function(x){x@nj})
