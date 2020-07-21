#' Maf class
#' @description Maf class.
#' @slot data data.table of MAF file containing somatic mutations.
#' @slot sample.info data.frame of sample information per patient.
#' @slot nonSyn.vc list of variant classifications which are considered as non-silent.
#' Default NULL, use Variant Classifications with "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"
#' @slot ref.build human reference genome version. Default: 'hg19'. Optional: 'hg18' or 'hg38'.
#' @rdname Maf-class
#' @exportClass Maf
Maf <- setClass(
  Class = "Maf",
  slots = c(
    data = 'data.table',
    sample.info = 'data.frame',
    nonSyn.vc = 'character',
    ref.build = 'character'
    
  )
)

#' MafList class
#' @description S4 class for storing a list of Maf objects.
#' @slot .Data a list of \code{\link{Maf}} objects.
#' @section Constructor:\describe{
#'                  \item{\code{MafList (...)}}{combine multiple Maf 
#'                  objects supplied in ... into a MafList object.}
#'                }
#' 
#' @rdname MafList-class
#' @exportClass MafList
MafList <- setClass(
  Class = "MafList",
  contains = "list"
)

## Prevent class 'phylo' from not existing
setClass('phylo')


#' phyloTree class
#' @name phyloTree-class
#' @aliases phyloTree
#' @description S4 class for storing informations about phylogenetic tree.
#' @slot patientID patient ID.
#' @slot tree a object of class "phylo".
#' @slot bootstrap.value a numeric vector of bootstrap values.
#' @slot method approach to construct a phylogenetic tree.
#' @slot binary.matrix a presense/absent binary matrix of mutations.
#' @slot ccf.matrix a ccf matrix of mutations.
#' @slot mut.branches a data.frame of mutations per trunk/branch.
#' @slot branch.type a data.frame of trunk/branch types based on shared pattern.
#' @slot ref.build human reference genome version. Default: 'hg19'. Optional: 'hg18' or 'hg38'.
#' @rdname phyloTree-class
#' @export
setClass(Class = 'phyloTree', slots = c(
  patientID = 'character', 
  tree = 'phylo',
  bootstrap.value = 'numeric',    
  method = 'character', 
  binary.matrix = 'matrix', 
  ccf.matrix = 'matrix', 
  mut.branches = 'data.frame', 
  branch.type = 'data.frame',
  ref.build = 'character'
  
))


#' phyloTreeList class
#' @description S4 class for storing a list of phyloTree objects.
#' @slot .Data a list of \code{\link{phyloTree}} objects.
#' @section Constructor:\describe{
#'                  \item{\code{phyloTreeList (...)}}{combine multiple phyloTree 
#'                  objects supplied in ... into a phyloTreeList object.}
#'                }
#' 
#' @rdname phyloTreeList-class
#' @exportClass phyloTreeList
setClass(Class = 'phyloTreeList', contains = "list")