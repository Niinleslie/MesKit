#--- Check phyloTree and filter patient 
subPhyloTree <- function(object, patient.id = NULL){
  if(is(object, "phyloTree")){
    phylotree_list <- list(object)
    names(phylotree_list) <- unique(getPhyloTreePatient(object))
    return(phylotree_list) 
  }else if(is(object, "phyloTreeList")){
    if(!is.null(patient.id)){
      patient.setdiff <- setdiff(patient.id, names(object))
      if(length(patient.setdiff) > 0){
        stop(paste0(patient.setdiff, " can not be found in phyloTreeList"))
      }
      object<- object[names(object)  %in% patient.id] 
    }
    return(object)
  }else{
    stop("Input should be phyloTree or phyloTreeList object generated by function getPhyloTree!")
  }
}