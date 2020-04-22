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
   branchAlias <- readPhyloTree(matTree)
   mut.branches <- treeMutationalBranches(patient.dat, branchAlias, binary.matrix)
   phylo.tree <- new('phyloTree', patientID = patientID, tree = matTree, 
                     binary.matrix = binary.matrix, ccf.matrix = ccf.matrix, 
                     mut.branches = mut.branches, refBuild = refBuild,
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

treeMutationalBranches <- function(maf.dat, branchAlias, binary.matrix){
   binary.matrix <- as.data.frame(binary.matrix)
   binary.matrix$mut.id <- rownames(binary.matrix)
   ## get mutationalSigs-related  infomation
   maf_input <- maf.dat
   branchChar <- as.character(branchAlias$Branch)
   datChr <- data.frame(chr=as.character(maf_input$Chromosome), stringsAsFactors=FALSE)
   datChr$chr <- paste("chr", datChr$chr, sep="")
   datMutgene <-  maf_input$Hugo_Symbol
   mutId <- dplyr::select(tidyr::unite(maf_input, "mut.id", 
                                       Hugo_Symbol, Chromosome, 
                                       Start_Position, 
                                       Reference_Allele, Tumor_Seq_Allele2, 
                                       sep=":"), mut.id)
   mutSigRef <- data.frame(as.character(maf_input$Tumor_Sample_Barcode), 
                           as.character(maf_input$Tumor_Type),
                           datChr,
                           maf_input$Start_Position, 
                           maf_input$End_Position,
                           maf_input$Reference_Allele, 
                           maf_input$Tumor_Seq_Allele2,
                           datMutgene, 
                           mutId,
                           stringsAsFactors=FALSE)
   colnames(mutSigRef) <- c("Branch_ID", 
                            "Branch_Tumor_Type",
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
      types <- unique(maf.dat[Tumor_Sample_Barcode %in% branch, ]$Tumor_Type)
      tsbs <- unique(maf.dat[Tumor_Sample_Barcode %in% branch, ]$Tumor_Sample_Barcode)
      tsbs.all <- unique(maf.dat$Tumor_Sample_Barcode)
      
      if(length(tsbs.all) == length(tsbs)){
         Branch_Tumor_Type <- "Public"
      } 
      else if(length(types) > 1){
         Branch_Tumor_Type <- paste0("Shared_",paste(types,collapse = "_"))
      }
      else if(length(types) == 1 & length(branch) > 1){
         Branch_Tumor_Type <- paste0("Shared_",paste(types,collapse = "_"))
      }else if(length(types) == 1 & length(branch) == 1){
         Branch_Tumor_Type <- paste0("Private_",paste(types,collapse = "_"))
      }
      
      ## generate mutation intersection for specific branch
      branch.intersection <- dplyr::intersect(
         binary.matrix %>% dplyr::filter_at(branch, dplyr::all_vars(. == 1)), 
         binary.matrix %>% dplyr::filter_at(unbranch, dplyr::all_vars(. == 0))) 
      ## special situation: branch.intersection NULL
      if (nrow(branch.intersection) == 0){
         # message(paste(branchName, ": There are no private mutations for branch ", sep=""))
         branch.mut <- data.frame(Branch_ID=branchName, 
                                  Branch_Tumor_Type = Branch_Tumor_Type,
                                  chr=NA,
                                  pos=NA,
                                  pos_end=NA,
                                  ref=NA,
                                  alt=NA,
                                  Hugo_Symbol=NA,
                                  mut_id="NoSigTag",
                                  Alias=as.character(branchAlias[which(branchAlias$Branch == branchName), ]$Alias))
         mutBranchesOutput[[branchName]] <- branch.mut
         next()
      }
      
      branch.mut.id <- branch.intersection$mut.id
      
      ## data duplication
      branch.mut <- mutSigRef[which(mutSigRef$mut_id %in% branch.mut.id), ]
      branch.mut$Branch_ID <- branchName
      
      ## Tumor_Type 
      branch.mut$Branch_Tumor_Type <- Branch_Tumor_Type
      
      # branch.mut$Branch_Tumor_Type <- paste(types,collapse = "_")
      branch.mut <- branch.mut[!duplicated(branch.mut),]
      branch.mut$Alias <- as.character(branchAlias[which(branchAlias$Branch == branchName), ]$Alias)
      
      ## generate branch mutation list
      mutBranchesOutput[[branchName]] <- branch.mut
   }
   return(mutBranchesOutput)
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
    mut.branches = 'list', 
    refBuild = 'character'
))

