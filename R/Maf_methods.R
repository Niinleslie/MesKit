#' @name showMafData
#' @rdname showMafData
#' @param object An object of classMaf
#' @return classMaf data
#' @exportMethod showMafData
#' @examples 
#' showMafData(maf)
setGeneric(name = "showMafData", function(object) standardGeneric("showMafData"))

#' @rdname showMafData
#' @aliases showMafData
setMethod(f = "showMafData",signature = "classMaf", function(object)object@data)


#' @name showSampleInfo
#' @rdname showSampleInfo
#' @param object An object of classMaf
#' @return sample information
#' @exportMethod showSampleInfo
#' @examples 
#' showSampleInfo(maf)
setGeneric(name = "showSampleInfo", function(object) standardGeneric("showSampleInfo"))

#' @rdname showSampleInfo
#' @aliases showSampleInfo
setMethod(f = "showSampleInfo",signature = "classMaf", function(object)object@sample.info)

#' @name showNonSyn_vc
#' @rdname showNonSyn_vc
#' @param object An object of classMaf
#' @return a list of Variant classifications which are considered as non-silent.
#' @exportMethod showNonSyn_vc
#' @examples 
#' showNonSyn_vc(maf)
setGeneric(name = "showNonSyn_vc", function(object) standardGeneric("showNonSyn_vc"))

#' @rdname showNonSyn_vc
#' @aliases showNonSyn_vc
setMethod(f = "showNonSyn_vc",signature = "classMaf", function(object)object@nonSyn.vc)

#' @name showRefBuild
#' @rdname showRefBuild
#' @param object An object of classMaf
#' @return Human reference genome versions of classMaf
#' @exportMethod showRefBuild
#' @examples 
#' showRefBuild(maf)
setGeneric(name = "showRefBuild", function(object) standardGeneric("showRefBuild"))

#' @rdname showRefBuild
#' @aliases showRefBuild
setMethod(f = "showRefBuild",signature = "classMaf", function(object)object@ref.build)

#' @name subsetMaf_list
#' @rdname subsetMaf_list
#' @param object An object of classMaf
#' @return Human reference genome versions of classMaf
#' @exportMethod subsetMaf_list
#' @examples 
#' subsetMaf_list(maf)
setGeneric(name = "subsetMaf_list", function(object,patient.id = NULL) standardGeneric("subsetMaf_list"))


#' @rdname subsetMaf_list
#' @aliases subsetMaf_list
setMethod(f = "subsetMaf_list",signature = "classMaf_list", 
          function(object,patient.id = NULL){
  if(!is.null(patient.id)){
    patient.setdiff <- setdiff(patient.id, names(object@patient.list))
    if(length(patient.setdiff) > 0){
      stop(paste0("Patient ", patient.setdiff, " can not be found in your data"))
    }
    object@patient.list <- object@patient.list[names(object@patient.list)  %in% patient.id]
  }
  return(object@patient.list)
}) 

