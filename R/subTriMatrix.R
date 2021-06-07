subTriMatrix <- function(phyloTree_list, CT = FALSE, level = 2){
  
  level <- as.character(level)
  level <- match.arg(level,
                     choices = c(1, 2, 3, 4, 5, 6),
                     several.ok = FALSE)
  
  bases <- c("A","C","G","T")
  if(CT){
    types <- c("C>A","C>G","C>T at CpG","C>T other","T>A","T>C","T>G")
  }else{
    types <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  }
  
  ## 96 trinucleotide 
  
  all_tri <- lapply(
    types,
    function(type){
      lapply(
        bases,
        function(base.up){
          lapply(
            bases,
            function(base.down){
              if(type == "C>T at CpG"){
                if(base.down != "G"){
                  return(NULL)
                }
              }
              tri <- paste(base.up,"[",type,"]",base.down,sep = "")
              return(tri)
            }
          )
        }
      )
    }
  ) %>% unlist()
  
  ref64_seq <- lapply(
    bases,
    function(base.mid){
      base.mid1 <- base.mid
      r1 <- lapply(
        bases,
        function(base.up){
          lapply(
            bases,
            function(base.down){
              tri <- paste(base.up, base.mid, base.down, sep = "")
              return(tri)
            }
          )
        }
      )
      return(r1)
    }
  ) %>% unlist()
  
  ref64_type <- lapply(
    bases,
    function(base.mid){
      base.mid1 <- base.mid
      r1 <- lapply(
        bases,
        function(base.up){
          lapply(
            bases,
            function(base.down){
              if(base.mid == "G"){
                base.mid1 <- "C"
              }
              else if(base.mid == "A"){
                base.mid1 <- "T"
              }
              n <- paste(base.up,base.mid1,base.down,sep = "")
              return(n)
            }
          )
        }
      )
      return(r1)
    }
  ) %>% unlist()
  
  names(ref64_type) <- ref64_seq
  ref64 <- ref64_type
  
  processTriMat <- function(phyloTree){
    patient <- getPhyloTreePatient(phyloTree)
    ## check reference
    refBuild <- getPhyloTreeRef(phyloTree)
    ref.options = c('hg18', 'hg19', 'hg38')
    if(!refBuild %in% ref.options){
      stop("'refBuild' can only be one of 'hg18', 'hg19' or 'hg38'")
    }else {
      refBuild <- paste("BSgenome.Hsapiens.UCSC.", refBuild, sep = "")
    }
    mut.branches <- phyloTree@mut.branches
    if(nrow(mut.branches) == 0){
      stop("There are not enough mutations in ",patient)
    }
    
    if(grepl("chr",mut.branches$Chromosome[1])){
      origin_context <- Biostrings::getSeq(get(refBuild),
                                           S4Vectors::Rle(mut.branches$Chromosome),
                                           mut.branches$Start_Position-1,
                                           mut.branches$Start_Position+1)
    }else{
      origin_context <- Biostrings::getSeq(get(refBuild),
                                           S4Vectors::Rle(paste("chr",mut.branches$Chromosome,sep = "")),
                                           mut.branches$Start_Position-1,
                                           mut.branches$Start_Position+1)
    }
    
    origin_context <- as.character(origin_context)
    
    muts <- paste(mut.branches$Reference_Allele, mut.branches$Tumor_Allele, sep = ">")
    
    mut_types <- gsub('G>T', 'C>A', muts)
    mut_types <- gsub('G>C', 'C>G', mut_types)
    mut_types <- gsub('G>A', 'C>T', mut_types)
    mut_types <- gsub('A>T', 'T>A', mut_types)
    mut_types <- gsub('A>G', 'T>C', mut_types)
    mut_types <- gsub('A>C', 'T>G', mut_types)
    
    idx <- which(muts != mut_types)
    c <- origin_context[idx]
    c <- IRanges::reverse(chartr("ATGC", "TACG", c))
    origin_context[idx] <- c
    
    mut.branches$mut_type <- mut_types
    mut.branches$origin_context <- origin_context
    
    if(!CT){
      mut.branches <- mut.branches %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(context = paste0(strsplit(.data$origin_context,"")[[1]][1],
                                       "[", .data$mut_type, "]",
                                       strsplit(.data$origin_context,"")[[1]][3])) %>% 
        as.data.frame()
    }else{
      CpG = c("ACG", "CCG", "TCG", "GCG")
      mut.branches <- mut.branches %>% 
        dplyr::rowwise() %>% 
        dplyr::mutate(context = dplyr::case_when(
          .data$mut_type == "C>T" & .data$origin_context %in% CpG ~  
            paste0(strsplit(.data$origin_context,"")[[1]][1],"[C>T at CpG]",strsplit(.data$origin_context,"")[[1]][3]),
          .data$mut_type == "C>T" & !.data$origin_context %in% CpG ~  
            paste0(strsplit(.data$origin_context,"")[[1]][1],"[C>T other]",strsplit(.data$origin_context,"")[[1]][3]),
          .data$mut_type != "C>T" ~ 
            paste0(strsplit(.data$origin_context,"")[[1]][1],"[", .data$mut_type, "]",strsplit(.data$origin_context,"")[[1]][3])
        )) %>% 
        as.data.frame()
    }
    
    ## initialize
    . <- NULL
    
    if(level == 1){
      mut.branches <- tidyr::separate_rows(mut.branches, .data$Branch_ID, sep = "&")
      
      branch_data_list <- list(mut.branches)
      names(branch_data_list) <- patient
    }else if(level == 2){
        mut.branches <- mut.branches %>% 
          tidyr::separate_rows(.data$Tumor_ID, sep = "&")
      branch_data_list <- split(mut.branches, mut.branches$Tumor_ID)
      
    }else if(level == 3){
      mut.branches <- tidyr::separate_rows(mut.branches, .data$Branch_ID, sep = "&") %>% 
        dplyr::mutate(Tumor_Sample_Barcode = .data$Branch_ID)
      
      branch_data_list <- split(mut.branches, mut.branches$Tumor_Sample_Barcode)
    }else if(level == 4){
      mut.branches <- mut.branches[mut.branches$Tree_Mut == TRUE,]
      branch_data_list <- split(mut.branches, mut.branches$Branch_ID)
      
    }else if(level == 5){
      ## sort mutation type by Public Shared Private
      mut.branches <- tidyr::separate_rows(mut.branches, .data$Branch_ID, sep = "&")
      
      mutation_type <- unique(mut.branches$Mutation_Type)
      public <- unique(mutation_type)[grep("Public", unique(mutation_type))] 
      shared <- sort(unique(mutation_type)[grep("Shared", unique(mutation_type))]) 
      private <- sort(unique(mutation_type)[grep("Private", unique(mutation_type))])
      mutation_type_level <- c(public, shared, private)
      mut.branches$Mutation_Type <- factor(mut.branches$Mutation_Type, levels = mutation_type_level)
      branch_data_list <- split(mut.branches, mut.branches$Mutation_Type)
    }else if(level == 6){
      mut.branches <- mut.branches[mut.branches$Tree_Mut == TRUE,]
      ## define Trunk
      branch_names <- unique(mut.branches$Branch_ID)
      branch_sample_num <- lapply(branch_names,function(x){
        s <- strsplit(x,"&")[[1]]
        num <- length(s)
      })
      
      trunk_name <- branch_names[which.max(branch_sample_num)]
      ## label the Trunk
      if (length(trunk_name) == 0){
        warning(paste0("Patient ", patient,": no trunk mutations were detected!"))
        return(NA)
      }
      
      mut.branches <- dplyr::mutate(
        mut.branches,
        "Branch_ID" = dplyr::if_else(
          .data$Branch_ID == trunk_name,
          "Trunk",
          "Branches"
        )
      )
      mut.branches$Branch_ID <- factor(mut.branches$Branch_ID, levels = c("Trunk","Branches"))
      branch_data_list <- split(mut.branches, mut.branches$Branch_ID)
    }

    tri_matrix_list <- lapply(
      branch_data_list,
      function(branch_data){
        branch_count <- table(branch_data$context)
        branch_count <- branch_count[names(branch_count) %in% all_tri] 
        
        m <- branch_count[all_tri]
        m[is.na(m)] <- 0
        names(m) <- all_tri
        
        branch_matrix <- matrix(m, ncol = length(all_tri), nrow = 1)
        colnames(branch_matrix) <- as.character(all_tri)
        return(as.data.frame(branch_matrix))
      }
    )
    
    
    tri_matrix <- dplyr::bind_rows(tri_matrix_list)
    row.names(tri_matrix) <- names(branch_data_list)
    
    tsb.label <- getPhyloTreeTsbLabel(phyloTree)
    
    if(nrow(tsb.label) == 0){
      return(as.matrix(tri_matrix))  
    }else{
      return(list(
        tri_matrix = as.matrix(tri_matrix),
        tsb.label = tsb.label
      ))  
    }
  }
  
  result <- lapply(phyloTree_list, processTriMat)
  
  return(result)
}