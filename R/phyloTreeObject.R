doGetPhyloTree <- function(patient.dat = NULL,
                           refBuild = NULL,
                           method = "NJ",
                           min.vaf = 0.02, 
                           min.CCF = NULL,
                           bootstrap.rep.num = 100){
   patientID <- unique(patient.dat$Patient_ID)
   if(!is.null(min.CCF)){
      
      if("CCF" %in% colnames(patient.dat)){
         patient.dat <- patient.dat[CCF > min.CCF, ]
      }
      else{
         warnings("min.CCF argument only works when CCF data is provided!")
         
      }
   }
   patient.dat <- patient.dat[which(patient.dat$VAF > min.vaf), ]
   ## information input
   binary.matrix <- getMutMatrix(patient.dat, use.ccf = FALSE)
   if("CCF" %in% colnames(patient.dat)){
      ccf.matrix <- getMutMatrix(patient.dat, use.ccf = TRUE)
   }else{
      ccf.matrix <- matrix() 
   }
   mut_dat <- t(binary.matrix)
   if(method == "NJ"){
      matTree <- nj(dist.gene(mut_dat))
      bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e)nj(dist.gene(e)),B = bootstrap.rep.num,quiet = T)/(bootstrap.rep.num)*100
   }else if(method == "MP"){
      matTree <- byMP(mut_dat)
      bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e){byMP(e)},B = bootstrap.rep.num,quiet = T)/(bootstrap.rep.num)*100 
   }else if(method == "ML"){
      matTree <- byML(mut_dat)
      bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e)byML(e),B = bootstrap.rep.num,quiet = T)/(bootstrap.rep.num)*100
   }else if(method == "FASTME.bal"){
      matTree <- ape::fastme.bal(dist.gene(mut_dat))
      bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e) ape::fastme.bal(dist.gene(e)),B = bootstrap.rep.num,quiet = T)/(bootstrap.rep.num)*100
   }else if(method == "FASTME.ols"){
      matTree <- ape::fastme.ols(dist.gene(mut_dat))
      bootstrap.value <- ape::boot.phylo(matTree, mut_dat, function(e) ape::fastme.ols(dist.gene(e)),B = bootstrap.rep.num,quiet = T)/(bootstrap.rep.num)*100
   }
   branch.id <- readPhyloTree(matTree)
   
   mut.branches_types <- treeMutationalBranches(patient.dat, branch.id, binary.matrix)
   mut.branches <- mut.branches_types$mut.branches
   branch.type <- mut.branches_types$branch.type
   
   phylo.tree <- new('phyloTree',
                     patientID = patientID, tree = matTree, 
                     binary.matrix = binary.matrix, ccf.matrix = ccf.matrix, 
                     mut.branches = mut.branches, branch.type = branch.type,
                     refBuild = refBuild,
                     bootstrap.value = bootstrap.value, method = method)
   return(phylo.tree)
}

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

treeMutationalBranches <- function(patient.dat, branch.id, binary.matrix){
   binary.matrix <- as.data.frame(binary.matrix)
   binary.matrix$mut.id <- rownames(binary.matrix)
   ## get mutationalSigs-related  infomation
   mafData <- patient.dat
   branchChar <- as.character(branch.id$Branch_ID)
   mutId <- dplyr::select(tidyr::unite(mafData, "mut.id", 
                                       Hugo_Symbol, Chromosome, 
                                       Start_Position, 
                                       Reference_Allele, Tumor_Seq_Allele2, 
                                       sep=":"), mut.id)
   mutSigRef <- data.frame(as.character(mafData$Tumor_Sample_Barcode), 
                           as.character(mafData$Tumor_ID),
                           paste("chr", as.character(mafData$Chromosome), sep=""),
                           mafData$Start_Position, 
                           mafData$End_Position,
                           mafData$Reference_Allele, 
                           mafData$Tumor_Seq_Allele2,
                           mafData$Hugo_Symbol, 
                           mutId,
                           stringsAsFactors=FALSE)
   colnames(mutSigRef) <- c("Branch_ID", 
                            "Branch_Tumor_ID",
                            "chr", "pos", 
                            "pos_end", "ref", 
                            "alt", "Hugo_Symbol", 
                            "mut_id")
   
   ## get branch infomation
   ls.branch <- branchChar[order(nchar(branchChar), branchChar)]
   branches <- strsplit(ls.branch, split='∩')
   
   ## generate mutational intersections for each branch
   mutBranchesOutput <- list()
   for (branch in branches){
      ## generate intersection's mut.id and get the mutation information in mutSigRef
      branch <- unlist(branch)
      ## generate the branch name
      branchName <- paste(branch, collapse="∩")
      branch.id <- append(branch, "mut.id")
      unbranch <- names(binary.matrix)[which(!(names(binary.matrix) %in% branch.id))]
      
      ## get branch tumor type
      types <- unique(patient.dat[Tumor_Sample_Barcode %in% branch, ]$Tumor_ID)
      tsbs <- unique(patient.dat[Tumor_Sample_Barcode %in% branch, ]$Tumor_Sample_Barcode)
      tsbs.all <- unique(patient.dat$Tumor_Sample_Barcode)
      
      if(length(tsbs.all) == length(tsbs)){
         Branch_Tumor_ID <- "Public"
      } 
      else if(length(types) > 1){
         Branch_Tumor_ID <- paste0("Shared_",paste(types,collapse = "_"))
      }
      else if(length(types) == 1 & length(branch) > 1){
         Branch_Tumor_ID <- paste0("Shared_",paste(types,collapse = "_"))
      }else if(length(types) == 1 & length(branch) == 1){
         Branch_Tumor_ID <- paste0("Private_",paste(types,collapse = "_"))
      }
      
      ## generate mutation intersection for specific branch
      branch.intersection <- dplyr::intersect(
         binary.matrix %>% dplyr::filter_at(branch, dplyr::all_vars(. == 1)), 
         binary.matrix %>% dplyr::filter_at(unbranch, dplyr::all_vars(. == 0))) 
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
      
      branch.mut.id <- branch.intersection$mut.id
      
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

