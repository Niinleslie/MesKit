#' @title  getMafData
#' @rdname getMafData
#' @param object An object of Maf
#' @return Maf data
#' @exportMethod getMafData
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' getMafData(maf$HCC5647)
setGeneric(name = "getMafData", function(object) standardGeneric("getMafData"))

#' @rdname getMafData
#' @aliases getMafData
setMethod(f = "getMafData", signature = "Maf", function(object)object@data)


#' @title getSampleInfo
#' @rdname getSampleInfo
#' @param object An object of Maf
#' @return Sample information
#' @exportMethod getSampleInfo
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' getSampleInfo(maf$HCC5647)
setGeneric(name = "getSampleInfo", function(object) standardGeneric("getSampleInfo"))

#' @rdname getSampleInfo
#' @aliases getSampleInfo
setMethod(f = "getSampleInfo", signature = "Maf", function(object)object@sample.info)

#' @title getNonSyn_vc
#' @rdname getNonSyn_vc
#' @param object An object of Maf
#' @return A list of Variant classifications which are considered as non-silent.
#' @exportMethod getNonSyn_vc
#' @examples
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' getNonSyn_vc(maf$HCC5647)
setGeneric(name = "getNonSyn_vc", function(object) standardGeneric("getNonSyn_vc"))

#' @rdname getNonSyn_vc
#' @aliases getNonSyn_vc
setMethod(f = "getNonSyn_vc", signature = "Maf", function(object)object@nonSyn.vc)

#' @title getMafRef
#' @rdname getMafRef
#' @param object An object of Maf
#' @return Human reference genome versions of Maf
#' @exportMethod getMafRef
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' getMafRef(maf$HCC5647)
setGeneric(name = "getMafRef", function(object) standardGeneric("getMafRef"))

#' @rdname getMafRef
#' @aliases getMafRef
setMethod(f = "getMafRef", signature = "Maf", function(object)object@ref.build)


#' @title getMafPatient
#' @rdname getMafPatient
#' @param object An object of Maf
#' @return Human reference genome versions of Maf
#' @exportMethod getMafPatient
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' getMafPatient(maf$HCC5647)
setGeneric(name = "getMafPatient", function(object) standardGeneric("getMafPatient"))

#' @rdname getMafPatient
#' @aliases getMafPatient
setMethod(f = "getMafPatient", signature = "Maf", function(object)unique(getMafData(object)$Patient_ID))


#' @title getTree
#' @rdname getTree
#' @param object An object of phyloTree
#' @return Tree object of phyloTree
#' @exportMethod getTree
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getTree(phyloTree$HCC5647)
setGeneric(name = "getTree", function(object) standardGeneric("getTree"))

#' @rdname getTree
#' @aliases getTree
setMethod(f = "getTree", signature = "phyloTree", function(object)object@tree)


#' @title getBinaryMatrix
#' @rdname getBinaryMatrix
#' @param object An object of phyloTree
#' @return Binary matrix of phyloTree
#' @exportMethod getBinaryMatrix
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getBinaryMatrix(phyloTree$HCC5647)
setGeneric(name = "getBinaryMatrix", function(object) standardGeneric("getBinaryMatrix"))

#' @rdname getBinaryMatrix
#' @aliases getBinaryMatrix
setMethod(f = "getBinaryMatrix", signature = "phyloTree", function(object)object@binary.matrix)

#' @title getPhyloTreePatient
#' @rdname getPhyloTreePatient
#' @param object An object of phyloTree
#' @return patientID of phyloTree
#' @exportMethod getPhyloTreePatient
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getPhyloTreePatient(phyloTree$HCC5647)
setGeneric(name = "getPhyloTreePatient", function(object) standardGeneric("getPhyloTreePatient"))

#' @rdname getPhyloTreePatient
#' @aliases getPhyloTreePatient
setMethod(f = "getPhyloTreePatient", signature = "phyloTree", function(object)object@patientID)

#' @title getBootstrapValue
#' @rdname getBootstrapValue
#' @param object An object of phyloTree
#' @return Bootstrap value of phyloTree
#' @exportMethod getBootstrapValue
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getBootstrapValue(phyloTree$HCC5647)
setGeneric(name = "getBootstrapValue", function(object) standardGeneric("getBootstrapValue"))

#' @rdname getBootstrapValue
#' @aliases getBootstrapValue
setMethod(f = "getBootstrapValue", signature = "phyloTree", function(object)object@bootstrap.value)

#' @title getTreeMethod
#' @rdname getTreeMethod
#' @param object An object of phyloTree
#' @return Tree construction method of phyloTree
#' @exportMethod getTreeMethod
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getTreeMethod(phyloTree$HCC5647)
setGeneric(name = "getTreeMethod", function(object) standardGeneric("getTreeMethod"))

#' @rdname getTreeMethod
#' @aliases getTreeMethod
setMethod(f = "getTreeMethod", signature = "phyloTree", function(object)object@method)

#' @title getCCFMatrix
#' @rdname getCCFMatrix
#' @param object An object of phyloTree
#' @return CCF matrix of phyloTree
#' @exportMethod getCCFMatrix
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getCCFMatrix(phyloTree$HCC5647)
setGeneric(name = "getCCFMatrix", function(object) standardGeneric("getCCFMatrix"))

#' @rdname getCCFMatrix
#' @aliases getCCFMatrix
setMethod(f = "getCCFMatrix", signature = "phyloTree", function(object)object@ccf.matrix)

#' @title getMutBranches
#' @rdname getMutBranches
#' @param object An object of phyloTree
#' @return Branches mutation of phyloTree
#' @exportMethod getMutBranches
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getMutBranches(phyloTree$HCC5647)
setGeneric(name = "getMutBranches", function(object) standardGeneric("getMutBranches"))

#' @rdname getMutBranches
#' @aliases getMutBranches
setMethod(f = "getMutBranches", signature = "phyloTree", function(object)object@mut.branches)

#' @title getBranchType
#' @rdname getBranchType
#' @param object An object of phyloTree
#' @return Branch type of phyloTree
#' @exportMethod getBranchType
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getBranchType(phyloTree$HCC5647)
setGeneric(name = "getBranchType", function(object) standardGeneric("getBranchType"))

#' @rdname getBranchType
#' @aliases getBranchType
setMethod(f = "getBranchType", signature = "phyloTree", function(object)object@branch.type)

#' @title getPhyloTreeRef
#' @rdname getPhyloTreeRef
#' @param object An object of phyloTree
#' @return Reference genome versions of phyloTree
#' @exportMethod getPhyloTreeRef
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getPhyloTreeRef(phyloTree$HCC5647)
setGeneric(name = "getPhyloTreeRef", function(object) standardGeneric("getPhyloTreeRef"))

#' @rdname getPhyloTreeRef
#' @aliases getPhyloTreeRef
setMethod(f = "getPhyloTreeRef", signature = "phyloTree", function(object)object@ref.build)



#' @title getPhyloTreeRef
#' @rdname getPhyloTreeTsbLabel
#' @param object An object of phyloTree
#' @return relationship between Tumor_Sample_Barcode and Tumor_Label
#' @exportMethod getPhyloTreeTsbLabel
#' @examples 
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' clin.File <- system.file("extdata/", "HCC_LDC.clin.txt", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, clinicalFile = clin.File, ccfFile=ccf.File, refBuild="hg19")
#' phyloTree <- getPhyloTree(maf)
#' getPhyloTreeTsbLabel(phyloTree$HCC5647)
setGeneric(name = "getPhyloTreeTsbLabel", function(object) standardGeneric("getPhyloTreeTsbLabel"))

#' @rdname getPhyloTreeRef
#' @aliases getPhyloTreeRef
setMethod(f = "getPhyloTreeTsbLabel", signature = "phyloTree", function(object)object@tsb.label)
