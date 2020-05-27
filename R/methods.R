#' @name showMafData
#' @rdname showMafData
#' @param object An object of Maf
#' @return Maf data
#' @exportMethod showMafData
#' @examples 
#' showMafData(maf)
setGeneric(name = "showMafData", function(object) standardGeneric("showMafData"))

#' @rdname showMafData
#' @aliases showMafData
setMethod(f = "showMafData",signature = "Maf", function(object)object@data)


#' @name showSampleInfo
#' @rdname showSampleInfo
#' @param object An object of Maf
#' @return sample information
#' @exportMethod showSampleInfo
#' @examples 
#' showSampleInfo(maf)
setGeneric(name = "showSampleInfo", function(object) standardGeneric("showSampleInfo"))

#' @rdname showSampleInfo
#' @aliases showSampleInfo
setMethod(f = "showSampleInfo",signature = "Maf", function(object)object@sample.info)

#' @name showNonSyn_vc
#' @rdname showNonSyn_vc
#' @param object An object of Maf
#' @return a list of Variant classifications which are considered as non-silent.
#' @exportMethod showNonSyn_vc
#' @examples 
#' showNonSyn_vc(maf)
setGeneric(name = "showNonSyn_vc", function(object) standardGeneric("showNonSyn_vc"))

#' @rdname showNonSyn_vc
#' @aliases showNonSyn_vc
setMethod(f = "showNonSyn_vc",signature = "Maf", function(object)object@nonSyn.vc)

#' @name showRefBuild
#' @rdname showRefBuild
#' @param object An object of Maf
#' @return Human reference genome versions of Maf
#' @exportMethod showRefBuild
#' @examples 
#' showRefBuild(maf)
setGeneric(name = "showRefBuild", function(object) standardGeneric("showRefBuild"))

#' @rdname showRefBuild
#' @aliases showRefBuild
setMethod(f = "showRefBuild",signature = "Maf", function(object)object@ref.build)

#' @name subsetMafList
#' @rdname subsetMafList
#' @param object An object of Maf
#' @return Human reference genome versions of Maf
#' @exportMethod subsetMafList
#' @examples 
#' subsetMafList(maf)
setGeneric(name = "subsetMafList", function(object,patient.id = NULL) standardGeneric("subsetMafList"))


#' @rdname subsetMafList
#' @aliases subsetMafList
setMethod(f = "subsetMafList",signature = "MafList", 
          function(object,patient.id = NULL){
  if(!is.null(patient.id)){
    patient.setdiff <- setdiff(patient.id, names(object))
    if(length(patient.setdiff) > 0){
      stop(paste0("Patient ", patient.setdiff, " can not be found in your data"))
    }
    object <- object[names(object)  %in% patient.id]
  }
  return(object)
}) 


#' @name showBinaryMatrix
#' @rdname showBinaryMatrix
#' @param object An object of phyloTree
#' @return binary matrix of phyloTree
#' @exportMethod showBinaryMatrix
#' @examples 
#' showBinaryMatrix(phyloTree)
setGeneric(name = "showBinaryMatrix", function(object) standardGeneric("showBinaryMatrix"))

#' @rdname showBinaryMatrix
#' @aliases showBinaryMatrix
setMethod(f = "showBinaryMatrix",signature = "phyloTree", function(object)object@binary.matrix)

#' @name showPatientID
#' @rdname showPatientID
#' @param object An object of phyloTree
#' @return patientID of phyloTree
#' @exportMethod showPatientID
#' @examples 
#' showPatientID(phyloTree)
setGeneric(name = "showPatientID", function(object) standardGeneric("showPatientID"))

#' @rdname showPatientID
#' @aliases showPatientID
setMethod(f = "showPatientID",signature = "phyloTree", function(object)object@patientID)

#' @name showBootstrapValue
#' @rdname showBootstrapValue
#' @param object An object of phyloTree
#' @return bootstrap value of phyloTree
#' @exportMethod showBootstrapValue
#' @examples 
#' showBootstrapValue(phyloTree)
setGeneric(name = "showBootstrapValue", function(object) standardGeneric("showBootstrapValue"))

#' @rdname showBootstrapValue
#' @aliases showBootstrapValue
setMethod(f = "showBootstrapValue",signature = "phyloTree", function(object)object@bootstrap.value)

#' @name showTreeMethod
#' @rdname showTreeMethod
#' @param object An object of phyloTree
#' @return tree construction method of phyloTree
#' @exportMethod showTreeMethod
#' @examples 
#' showTreeMethod(phyloTree)
setGeneric(name = "showTreeMethod", function(object) standardGeneric("showTreeMethod"))

#' @rdname showTreeMethod
#' @aliases showTreeMethod
setMethod(f = "showTreeMethod",signature = "phyloTree", function(object)object@method)

#' @name showCCFMatrix
#' @rdname showCCFMatrix
#' @param object An object of phyloTree
#' @return ccf matrix of phyloTree
#' @exportMethod showCCFMatrix
#' @examples 
#' showCCFMatrix(phyloTree)
setGeneric(name = "showCCFMatrix", function(object) standardGeneric("showCCFMatrix"))

#' @rdname showCCFMatrix
#' @aliases showCCFMatrix
setMethod(f = "showCCFMatrix",signature = "phyloTree", function(object)object@ccf.matrix)

#' @name showMutBranches
#' @rdname showMutBranches
#' @param object An object of phyloTree
#' @return branches mutation of phyloTree
#' @exportMethod showMutBranches
#' @examples 
#' showMutBranches(phyloTree)
setGeneric(name = "showMutBranches", function(object) standardGeneric("showMutBranches"))

#' @rdname showMutBranches
#' @aliases showMutBranches
setMethod(f = "showMutBranches",signature = "phyloTree", function(object)object@mut.branches)

#' @name showBranchType
#' @rdname showBranchType
#' @param object An object of phyloTree
#' @return branch type of phyloTree
#' @exportMethod showBranchType
#' @examples 
#' showBranchType(phyloTree)
setGeneric(name = "showBranchType", function(object) standardGeneric("showBranchType"))

#' @rdname showBranchType
#' @aliases showBranchType
setMethod(f = "showBranchType",signature = "phyloTree", function(object)object@branch.type)




#' @name subsetPhyloTreeList
#' @rdname subsetPhyloTreeList
#' @param object An object of phyloTreeList
#' @return select the specific patients in phyloTreeList
#' @exportMethod subsetPhyloTreeList
#' @examples 
#' subsetPhyloTreeList(phyloTree)
setGeneric(name = "subsetPhyloTreeList", function(object,patient.id = NULL) standardGeneric("subsetPhyloTreeList"))

#' @rdname subsetPhyloTreeList
#' @aliases subsetPhyloTreeList
setMethod(f = "subsetPhyloTreeList",signature = "phyloTreeList", function(object,patient.id = NULL){
    if(!is.null(patient.id)){
        patient.setdiff <- setdiff(patient.id, names(object))
        if(length(patient.setdiff) > 0){
            stop(paste0(patient.setdiff, " can not be found in your data"))
        }
        object<- object[names(object)  %in% patient.id] 
    }
    return(object)
}) 
