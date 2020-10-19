byMP <- function(mut_dat){
   matTree <- nj(dist.gene(mut_dat))
   tree_dat <- phangorn::as.phyDat(mut_dat, type="USER", levels = c(0, 1))
   tree_pars <- suppressMessages(phangorn::optim.parsimony(matTree, tree_dat, method = "sankoff", trace = FALSE)) 
   matTree <- phangorn::acctran(tree_pars, tree_dat)
   return(matTree)
}

byML <- function(mut_dat){
   matTree <- nj(dist.gene(mut_dat))
   tree_dat <- phangorn::as.phyDat(mut_dat, type="USER", levels = c(0, 1))
   fitJC <- phangorn::pml(matTree, tree_dat)
   fitJC <- try(phangorn::optim.pml(fitJC,control = phangorn::pml.control(trace = FALSE))) 
   matTree <- fitJC$tree
   return(matTree)
}

treeMutationalBranches <- function(maf_data, branch.id, binary.matrix){
   binary.matrix <- as.data.frame(binary.matrix)
   binary.matrix$mut_id <- rownames(binary.matrix)
   ## get mutationalSigs-related  infomation
   mut_id <- dplyr::select(tidyr::unite(maf_data, "mut_id", 
                                       "Hugo_Symbol", "Chromosome", 
                                       "Start_Position", 
                                       "Reference_Allele", "Tumor_Seq_Allele2", 
                                       sep=":"), "mut_id")
   mutSigRef <- data.frame(Branch_ID = as.character(maf_data$Tumor_Sample_Barcode), 
                           Mutation_Type = as.character(maf_data$Tumor_ID),
                           Tumor_ID = as.character(maf_data$Tumor_ID),
                           Hugo_Symbol = maf_data$Hugo_Symbol, 
                           Chromosome = as.character(maf_data$Chromosome),
                           Start_Position = maf_data$Start_Position, 
                           End_Position = maf_data$End_Position,
                           Reference_Allele = maf_data$Reference_Allele, 
                           Tumor_Allele  = maf_data$Tumor_Seq_Allele2,
                           mut_id = mut_id,
                           stringsAsFactors=FALSE)
   
   # if("Tumor_Sample_Label" %in% colnames(maf_data)){
   #    mutSigRef$Branch_Label <- maf_data$Tumor_Sample_Label
   # }
   
   ## get branch infomation
   branchChar <- as.character(branch.id$Branch_ID)
   ls.branch <- branchChar[order(nchar(branchChar), branchChar)]
   branches <- strsplit(ls.branch, split='&')
   
   ## generate mutational intersections for each branch
   mutBranchesOutput <- list()
   for (branch in branches){
      ## generate intersection's mut_id and get the mutation information in mutSigRef
      branch <- unlist(branch)
      ## generate the branch name
      branchName <- paste(branch, collapse="&")
      # branch.id <- append(branch, "mut_id")
      unbranch <- names(binary.matrix)[!names(binary.matrix) %in% branch]
      unbranch <- unbranch[unbranch!="mut_id"]
      ## get branch tumor type
      ids_df <- maf_data %>% 
         dplyr::filter(.data$Tumor_Sample_Barcode %in% branch) %>% 
         dplyr::distinct(.data$Tumor_Sample_Barcode, .keep_all = TRUE)
      tumor_ids <- ids_df$Tumor_ID
      mutation_type_ids <- unique(maf_data[maf_data$Tumor_Sample_Barcode %in% branch, ]$Tumor_ID)
      tsbs <- unique(maf_data[maf_data$Tumor_Sample_Barcode %in% branch, ]$Tumor_Sample_Barcode)
      tsbs.all <- unique(maf_data$Tumor_Sample_Barcode)
      
      if(length(tsbs.all) == length(tsbs)){
         Mutation_Type <- "Public"
         # tumor_ids <- paste(ids, collapse = "&")
      }else if(length(mutation_type_ids) > 1){
         Mutation_Type <- paste0("Shared_", paste(mutation_type_ids,collapse = "_"))
         # tumor_ids <- paste(ids, collapse = "&")
      }else if(length(mutation_type_ids) == 1 & length(branch) > 1){
          if(length(unique(maf_data$Tumor_ID)) == 1){
              Mutation_Type <- "Shared"
          }else{
              Mutation_Type <- paste0("Shared_",paste(mutation_type_ids,collapse = "_"))
          }
         # tumor_ids <- ids
      }else if(length(mutation_type_ids) == 1 & length(branch) == 1){
          if(length(unique(maf_data$Tumor_ID)) == 1){
              Mutation_Type <- "Private"
          }else{
              Mutation_Type <- paste0("Private_",paste(mutation_type_ids,collapse = "_"))
          }
         # tumor_ids <- ids
      }
      tumor_ids <- paste(tumor_ids, collapse = "&")
      ## initialize 
      . = NULL
      ## generate mutation intersection for specific branch
      branch.intersection <- dplyr::intersect(
         binary.matrix %>% dplyr::filter_at(branch, dplyr::all_vars(. == 1)), 
         binary.matrix %>% dplyr::filter_at(unbranch, dplyr::all_vars(. == 0)))
      # print(branch)
      # print(branch.intersection)
      ## special situation: branch.intersection NULL
      # if (nrow(branch.intersection) == 0){
      #    left_binary_matrix <- binary.matrix %>% dplyr::filter_at(branch, dplyr::all_vars(. == 1))
      #    left_mut_ids <- left_binary_matrix$mut_id
      #    branch.mut <- mutSigRef[which(mutSigRef$mut_id %in% left_mut_ids), ]
      #    branch.mut$Branch_ID <- paste0(branchName, "(NA)") 
      #    
      #    ## Tumor_ID 
      #    branch.mut$Mutation_Type <- paste0(Mutation_Type, "(NA)") 
      #    # message(paste(branchName, ": There are no private mutations for branch ", sep=""))
      #    mutBranchesOutput[[branchName]] <- branch.mut
      #    next()
      # }
      
      if (nrow(branch.intersection) == 0){
         # message(paste(branchName, ": There are no private mutations for branch ", sep=""))
         branch.mut <- data.frame(Branch_ID=branchName,
                                  Mutation_Type = Mutation_Type,
                                  Tumor_ID = tumor_ids,
                                  Chromosome=NA,
                                  Start_Position=NA,
                                  End_Position=NA,
                                  Reference_Allele=NA,
                                  Tumor_Allele=NA,
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
      branch.mut$Mutation_Type <- Mutation_Type
      branch.mut$Tumor_ID <- tumor_ids
      
      # branch.mut$Mutation_Type <- paste(ids,collapse = "_")
      branch.mut <- branch.mut[!duplicated(branch.mut),]
      # branch.mut$Alias <- as.character(branch.id[which(branch.id$Branch == branchName), ]$Alias)
      
      ## generate branch mutation list
      mutBranchesOutput[[branchName]] <- branch.mut
   }
   
   # print(binary.matrix[!binary.matrix$mut_id %in% mutBranchesOutput$mut_id,])
   
   mut.branches <- dplyr::bind_rows(mutBranchesOutput)%>% 
      dplyr::filter(!is.na(.data$Chromosome))
   
   mut_left <- mutSigRef %>%
      dplyr::filter(!.data$mut_id %in% mut.branches$mut_id) %>% 
      dplyr::mutate(Branch_ID = paste0(.data$Branch_ID,"(NA)"),
                    Mutation_Type = paste0(.data$Mutation_Type,"(NA)")) %>% 
      dplyr::select(-"mut_id")
   
   
   
   mut.branches <- dplyr::select(mut.branches, -"mut_id") %>% 
      dplyr::bind_rows(mut_left)
   
   branch.type <- mut.branches %>% 
       dplyr::select("Branch_ID", "Mutation_Type") %>% 
       dplyr::distinct(.data$Branch_ID, .keep_all = TRUE) %>% 
      dplyr::filter(!grepl("\\(NA\\)", .data$Branch_ID))
   
   
   return(list(
       mut.branches = mut.branches,
       branch.type = branch.type
   ))
}
