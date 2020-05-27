#' Set Maf class
#' @description Maf class.
#' @name Maf
#' @aliases Maf
#' @docType class
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
#' @description S4 class for storing multi-patients Maf object.
#' @name MafList
#' @aliases MafList
#' @docType class
#' @rdname MafList-class
#' @exportClass MafList
MafList <- setClass(
  Class = "MafList",
  contains = "list"
)

## Prevent class 'phylo' from not existing
setClass('phylo')


#' phyloTree class
#' @description S4 class for storing informations about phylogenetic tree .
#' @name phyloTree
#' @aliases phyloTree
#' @docType phyloTree
#' @rdname phyloTree-class
#' @exportClass phyloTree
setClass('phyloTree', slots = c(
  patientID = 'character', 
  tree = 'phylo',
  bootstrap.value = 'numeric',    
  method = 'character', 
  binary.matrix = 'matrix', 
  ccf.matrix = 'matrix', 
  mut.branches = 'data.frame', 
  branch.type = 'data.frame',
  refBuild = 'character'
  
))


#' phyloTreeList class
#' @description S4 class for storing multi-patients phyloTree.
#' @name phyloTreeList
#' @aliases phyloTreeList
#' @docType phyloTreeList
#' @rdname phyloTreeList-class
#' @exportClass phyloTreeList
setClass('phyloTreeList', contains = "list")