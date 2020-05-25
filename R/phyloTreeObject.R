byMP <- function(mut_dat){
   matTree <- nj(dist.gene(mut_dat))
   tree_dat <- phangorn::as.phyDat(mut_dat, type="USER", levels = c(0, 1))
   tree_pars <- suppressMessages(phangorn::optim.parsimony(matTree, tree_dat,method = "sankoff",trace = F)) 
   matTree <- phangorn::acctran(tree_pars, tree_dat)
   return(matTree)
}

byML <- function(mut_dat){
   matTree <- nj(dist.gene(mut_dat))
   tree_dat <- phangorn::as.phyDat(mut_dat, type="USER", levels = c(0, 1))
   fitJC <- phangorn::pml(matTree, tree_dat)
   fitJC <- try(phangorn::optim.pml(fitJC,control = phangorn::pml.control(trace = F))) 
   matTree <- fitJC$tree
   return(matTree)
}

treeMutationalBranches <- function(maf_data, branch.id, binary.matrix){
   binary.matrix <- as.data.frame(binary.matrix)
   binary.matrix$mut_id <- rownames(binary.matrix)
   ## get mutationalSigs-related  infomation
   mut_id <- dplyr::select(tidyr::unite(maf_data, "mut_id", 
                                       Hugo_Symbol, Chromosome, 
                                       Start_Position, 
                                       Reference_Allele, Tumor_Seq_Allele2, 
                                       sep=":"), mut_id)
   mutSigRef <- data.frame(Branch_ID = as.character(maf_data$Tumor_Sample_Barcode), 
                           Branch_Tumor_ID = as.character(maf_data$Tumor_ID),
                           chr = as.character(maf_data$Chromosome),
                           pos = maf_data$Start_Position, 
                           pos_end = maf_data$End_Position,
                           ref = maf_data$Reference_Allele, 
                           alt = maf_data$Tumor_Seq_Allele2,
                           Hugo_Symbol = maf_data$Hugo_Symbol, 
                           mut_id = mut_id,
                           stringsAsFactors=FALSE)
   
   ## get branch infomation
   branchChar <- as.character(branch.id$Branch_ID)
   ls.branch <- branchChar[order(nchar(branchChar), branchChar)]
   branches <- strsplit(ls.branch, split='∩')
   
   ## generate mutational intersections for each branch
   mutBranchesOutput <- list()
   for (branch in branches){
      ## generate intersection's mut_id and get the mutation information in mutSigRef
      branch <- unlist(branch)
      ## generate the branch name
      branchName <- paste(branch, collapse="∩")
      # branch.id <- append(branch, "mut_id")
      unbranch <- names(binary.matrix)[!names(binary.matrix) %in% branch]
      unbranch <- unbranch[unbranch!="mut_id"]
      ## get branch tumor type
      types <- unique(maf_data[Tumor_Sample_Barcode %in% branch, ]$Tumor_ID)
      tsbs <- unique(maf_data[Tumor_Sample_Barcode %in% branch, ]$Tumor_Sample_Barcode)
      tsbs.all <- unique(maf_data$Tumor_Sample_Barcode)
      
      if(length(tsbs.all) == length(tsbs)){
         Branch_Tumor_ID <- "Public"
      }else if(length(types) > 1){
         if(length(unique(maf_data$Tumor_ID)) == 1){
             Branch_Tumor_ID <- "Shared"
         }else{
             Branch_Tumor_ID <- paste0("Shared_",paste(types,collapse = "_"))
         }
      }else if(length(types) == 1 & length(branch) > 1){
          if(length(unique(maf_data$Tumor_ID)) == 1){
              Branch_Tumor_ID <- "Shared"
          }else{
              Branch_Tumor_ID <- paste0("Shared_",paste(types,collapse = "_"))
          }
      }else if(length(types) == 1 & length(branch) == 1){
          if(length(unique(maf_data$Tumor_ID)) == 1){
              Branch_Tumor_ID <- "Private"
          }else{
              Branch_Tumor_ID <- paste0("Private_",paste(types,collapse = "_"))
          }
      }
      
      ## generate mutation intersection for specific branch
      branch.intersection <- dplyr::intersect(
         binary.matrix %>% dplyr::filter_at(branch, dplyr::all_vars(. == 1)), 
         binary.matrix %>% dplyr::filter_at(unbranch, dplyr::all_vars(. == 0)))
      # print(branch)
      # print(branch.intersection)
      ## special situation: branch.intersection NULL
      if (nrow(branch.intersection) == 0){
         # message(paste(branchName, ": There are no private mutations for branch ", sep=""))
         branch.mut <- data.frame(Branch_ID=branchName,
                                  Branch_Tumor_ID = Branch_Tumor_ID,
                                  chr=NA,
                                  pos=NA,
                                  pos_end=NA,
                                  ref=NA,
                                  alt=NA,
                                  Hugo_Symbol=NA,
                                  mut_id=NA
                                  # Alias=NA
                                  )
         mutBranchesOutput[[branchName]] <- branch.mut
         next()
      }
      
      branch.mut.id <- branch.intersection$mut_id
      ## data duplication
      branch.mut <- mutSigRef[which(mutSigRef$mut_id %in% branch.mut.id), ]
      branch.mut$Branch_ID <- branchName
      
      ## Tumor_ID 
      branch.mut$Branch_Tumor_ID <- Branch_Tumor_ID
      
      # branch.mut$Branch_Tumor_ID <- paste(types,collapse = "_")
      branch.mut <- branch.mut[!duplicated(branch.mut),]
      # branch.mut$Alias <- as.character(branch.id[which(branch.id$Branch == branchName), ]$Alias)
      
      ## generate branch mutation list
      mutBranchesOutput[[branchName]] <- branch.mut
   }
   
   # print(binary.matrix[!binary.matrix$mut_id %in% mutBranchesOutput$mut_id,])
   
   mutBranchesOutput <- plyr::rbind.fill(mutBranchesOutput) %>% 
       dplyr::select(-mut_id)
   
   branch.type <- mutBranchesOutput %>% 
       dplyr::select(Branch_ID, Branch_Tumor_ID)
   
   mut.branches <- mutBranchesOutput %>% 
       dplyr::filter(!is.na(chr))
   
   return(list(
       mut.branches = mut.branches,
       branch.type = branch.type
   ))
}



#Prevent class 'phylo' from not existing
setClass('phylo')
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




#' @name subset_phyloTree_list
#' @rdname subset_phyloTree_list
#' @param object An object of phyloTree_list
#' @return select the specific patients in phyloTree_list
#' @exportMethod subset_phyloTree_list
#' @examples 
#' subset_phyloTree_list(phyloTree)
setGeneric(name = "subset_phyloTree_list", function(object,patient.id = NULL) standardGeneric("subset_phyloTree_list"))

#' @rdname subset_phyloTree_list
#' @aliases subset_phyloTree_list
setMethod(f = "subset_phyloTree_list",signature = "phyloTree_list", function(object,patient.id = NULL){
    if(!is.null(patient.id)){
        patient.setdiff <- setdiff(patient.id, names(object@patient.list))
        if(length(patient.setdiff) > 0){
            stop(paste0(patient.setdiff, " can not be found in your data"))
        }
        object@patient.list <- object@patient.list[names(object@patient.list)  %in% patient.id] 
    }
    return(object@patient.list)
}) 
