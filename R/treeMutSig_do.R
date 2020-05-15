doTreeMutSig <- function(phyloTree,
                         geneList=NULL,
                         min.mut.count=15,
                         signaturesRef="cosmic_v2",
                         associated = c(), 
                         signature.cutoff = 0.1,
                         withinTumor = FALSE){
   ## get branches information from phyloTree object
   mutSigRef <- phyloTree@mut.branches
   
   patientID <- phyloTree@patientID
   refBuild <- phyloTree@refBuild
   ref.options = c('hg18', 'hg19', 'hg38')
   if(!refBuild %in% ref.options){
      stop("refBuild can only be either 'hg18', 'hg19' or 'hg38'")
   }else {
      refBuild <- paste("BSgenome.Hsapiens.UCSC.", refBuild, sep = "")
   }
   
   
   ## Select mutations in selected genes
   if(!is.null(geneList)){
      mutSigRef <- mutSigRef %>% 
         dplyr::filter(Hugo_Symbol %in% geneList)
   }
   
   ## get the mutational signature of the branch 
   mutSigsOutput <- data.frame()
   spectrum_df <- data.frame()
   
   ## count 96 substitution typs in each branches
   sigsInput <- countTriplet(mutSigRef = mutSigRef,
                             withinTumor = withinTumor,
                             refBuild = refBuild,
                             patientID = patientID,
                             CT = FALSE)
   df.aetiology <- NULL
   if(class(signaturesRef) == 'character'){
       signaturesRef <- match.arg(signaturesRef,
                                  choices = c("cosmic_v2","nature2013","genome_cosmic_v3","exome_cosmic_v3"),
                                  several.ok = FALSE)
       signatures.aetiology <- readRDS(file = system.file("extdata", "signatures.aetiology.rds", package = "MesKit")) 
       if (signaturesRef == "cosmic_v2"){
           sigsRef <- readRDS(file = system.file("extdata", "signatures.cosmic.rds", package = "MesKit"))
           df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v2$aetiology,
                                      sig = rownames(signatures.aetiology$cosmic_v2))
       }else if(signaturesRef == "nature2013"){
           sigsRef <- readRDS(file = system.file("extdata", "signatures.nature2013.rds", package = "MesKit"))
           df.aetiology <- data.frame(aeti = signatures.aetiology$nature2013$aetiology,
                                      sig = rownames(signatures.aetiology$nature2013))
       }else if(signaturesRef == "genome_cosmic_v3"){
           sigsRef <- readRDS(file = system.file("extdata", "signatures.genome.cosmic.v3.may2019.rds", package = "MesKit"))
           df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v3$aetiology,
                                      sig = rownames(signatures.aetiology$cosmic_v3))
       }else if(signaturesRef == "exome_cosmic_v3"){
           sigsRef <- readRDS(file = system.file("extdata", "signatures.exome.cosmic.v3.may2019.rds", package = "MesKit"))
           df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v3$aetiology,
                                      sig = rownames(signatures.aetiology$cosmic_v3))
       }
   }else if(class(signaturesRef)!= 'data.frame'){
       stop('Input signature reference needs to be a data frame')
   }else{
       sigsRef <- signaturesRef 
   }
   
   ## subset of the signatures reference
   if(!is.null(associated)){
       signature.setdiff <- setdiff(associated, rownames(sigsRef))
       if(length(signature.setdiff) > 0){
           stop(paste0(signature.setdiff, " can not be found in signature reference"))
       }
       sigsRef <- sigsRef[rownames(sigsRef) %in% associated, ]
   }
   
   ## cal cos sim signature and branch 
   cos_sim.mat <- calSim(sigsInput = sigsInput[which(rowSums(sigsInput)!=0),], sigsRef = sigsRef)
   
   if(withinTumor){
      branchesName <- unique(mutSigRef$Branch_Tumor_ID) 
   }else{
       branchesName <- unique(mutSigRef$Branch_ID) 
   }
   
   for (branchCounter in length(branchesName):1){
      ## generate a single branch
      branchName <- branchesName[branchCounter]
      # print(branchName)
      if(withinTumor){
         branch.mut.num <- length(which(mutSigRef$Branch_Tumor_ID == branchName))
      }else{
         branch.mut.num <-  length(which(mutSigRef$Branch_ID == branchName))
      }
      ## get the mutational signature of the branch
      if (branch.mut.num < min.mut.count){
         sig_names <- "noMapSig"
         sig_weights <- 0
         message(paste0("Warnings: ",
                        "mutation number of Branch " , branchName,
                        " is less than the min.mut.count argument! ",
                        "This branch will be skipped")
      )
         # if (any(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$mut_id == "NoSigTag")) {
         #    mut.count <- 0
         # }
        if(withinTumor){
           mut.count <- length(mutSigRef[which(mutSigRef$Branch_Tumor_ID == branchName), 1])
        }
        else{
           mut.count <- length(mutSigRef[which(mutSigRef$Branch_ID == branchName), 1])
        }
      }else{
          
          ## get mutation signature
         fit <- fitSignature(sigsInput = sigsInput[branchName,], sigsRef = sigsRef)
         fit$Reconstructed <- fit$Reconstructed/sum(fit$Reconstructed)
         fit$weight <- fit$weight/sum(fit$weight)

         idx <- which(fit$weight[1,] > signature.cutoff)
         sig_names <- colnames(fit$weight)[idx]
         sig_names <- gsub('[.]', ' ', sig_names)
         sig_weights <- fit$weight[,idx]
         sig_weights <- as.numeric(sig_weights)
         
         if(withinTumor){
            mut.count <- length(mutSigRef[which(mutSigRef$Branch_Tumor_ID == branchName), 1])
         }else{
            mut.count <- length(mutSigRef[which(mutSigRef$Branch_ID == branchName), 1])
         }
         
         ## get spectrum of reconstructed
         reconstructed_spectrum <- as.data.frame(fit$Reconstructed)
         reconstructed_spectrum$Branch <- as.character(row.names(fit$Reconstructed)) 
         reconstructed_spectrum <- tidyr::pivot_longer(reconstructed_spectrum ,
                                                   -Branch,
                                                   names_to = "group",
                                                   values_to = "sig_weights") %>% 
                              as.data.frame()
         reconstructed_spectrum$spectrum.type <- "Reconstructed"
         
         ## get Original fractions of 96 types
         original_spectrum <- as.data.frame(sigsInput[branchName,]/sum(sigsInput[branchName,]))
         original_spectrum$Branch <- as.character(row.names(fit$Reconstructed)) 
         original_spectrum <- tidyr::pivot_longer(original_spectrum ,
                                           -Branch,
                                           names_to = "group",
                                           values_to = "sig_weights") %>% 
                       as.data.frame()
         original_spectrum$spectrum.type <- "Original"
         
         reconstructed_spectrum <- rbind(reconstructed_spectrum,original_spectrum)
         
         # if(!withinTumor){
         #     reconstructed_spectrum$Alias <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Alias))
         # }
         # 
         
         ## sort signature by weights
         sig_order <- order(sig_weights,decreasing = T)
         sig_weights <- sig_weights[sig_order]
         sig_names <- sig_names[sig_order]
         ##  title
         t <- ''
         for(i in  1:length(sig_weights)){
             if(i == 1){
                 t <- paste(sig_names[i], ": ", round(sig_weights[i],3), sep = "")
                 next
             }
            t <- paste(t, " & ", sig_names[i], ": ", round(sig_weights[i],3), sep = "")
         }
         reconstructed_spectrum$Signature <- t
         spectrum_df <- rbind(spectrum_df, reconstructed_spectrum)
      }
      
      mutSigsBranch <- data.frame(
         sig = sig_names, 
         weight = format(round(sig_weights, digits = 4)))
      
      if(withinTumor){
         mutSigsBranch$mut.count <- mut.count
         mutSigsBranch$Branch_Tumor_ID <- as.character(branchName)
         mutSigsBranch <- dplyr::select(mutSigsBranch,
                                        Branch_Tumor_ID,
                                        sig, 
                                        mut.count,
                                        weight)
         
      }else{
         mutSigsBranch$branch <- c(branchName)
         # mutSigsBranch$Alias <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Alias))
         mutSigsBranch$mut.count <- mut.count
         mutSigsBranch$Branch_Tumor_ID <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Branch_Tumor_ID))
         
         mutSigsBranch <- dplyr::select(mutSigsBranch,
                                        branch, sig, mut.count,
                                        weight, Branch_Tumor_ID)
         
      }
      
      
      ## collect branches' mutataional signature information
      mutSigsOutput <- rbind(mutSigsOutput, mutSigsBranch)
   }
   mutSigsOutput$weight <- as.numeric(as.vector(mutSigsOutput$weight))
   
   ## calculation process(maybe could be replaced by lapply)
   if(nrow(spectrum_df) > 0){
       # if(withinTumor){
       #     spectrum_df$Alias <- spectrum_df$Branch
       # }
       if(withinTumor){
           branch_tumor_ids <- unique(spectrum_df$Branch) 
           public <- branch_tumor_ids[grep("Public", branch_tumor_ids)] 
           shared <- branch_tumor_ids[grep("Shared", branch_tumor_ids)] 
           private <- branch_tumor_ids[grep("Private", branch_tumor_ids)]
           type.level <- c(public, shared, private)
           spectrum_df$Patient_ID <- patientID
           colnames(spectrum_df) <- c("Branch_Tumor_ID", "Group", "Mutation_Probability","spectrum.type","Signature","Patient_ID")
           rownames(spectrum_df) <- 1:nrow(spectrum_df)
       }else{
           spectrum_df$Patient_ID <- patientID
           colnames(spectrum_df) <- c("Branch", "Group", "Mutation_Probability","spectrum.type", "Signature","Patient_ID")
           rownames(spectrum_df) <- 1:nrow(spectrum_df)
       }
   }
   
   cos_sim.mat <- as.matrix(cos_sim.mat)
   # if(withinTumor){
   #     Alias <- unique(mutSigsOutput$Branch_Tumor_ID)
   #     names(Alias) <- unique(mutSigsOutput$Branch_Tumor_ID)
   # }else{
   #     Alias <- unique(mutSigsOutput$Alias)
   #     names(Alias) <- unique(mutSigsOutput$branch)
   # }
   # rownames(cos_sim.mat) <- Alias[rownames(cos_sim.mat)] 
   message(paste0("Patient ", phyloTree@patientID, ": mutational signatures identified successfully!"))
   treeMSOutput <- list(
      sigsInput=sigsInput, 
      spectrum_df=spectrum_df,
      cos_sim.mat = cos_sim.mat,
      mutSigsOutput=mutSigsOutput, 
      df.aetiology=df.aetiology, 
      patientID = patientID)
   return(treeMSOutput)
}

countTriplet <- function(mutSigRef, withinTumor, refBuild, patientID, CT){
    
    if(nrow(mutSigRef) == 0){
        stop("Error: There are not enough mutations in ",patientID)
    }
    
    bases <- c("A","C","G","T")
    if(CT){
        types <- c("C>A","C>G","C>T at CpG","C>T other","T>A","T>C","T>G")
    }else{
        types <- c("C>A","C>G","C>T","T>A","T>C","T>G")
    }
    ## 96 substitution types
    triplet96 <- c()
    seq96 <- c()
    for(type in types){
        for(base.up in bases){
            for(base.down in bases){
                tri <- paste(base.up,"[",type,"]",base.down,sep = "")
                triplet96 <- append(triplet96,tri)
                if(type %in% types[1:(round(length(types)/2))]){
                    seqname <- paste(base.up,"C",base.down,sep = "")
                    seq96 <- append(seq96,seqname)
                }
                else{
                    seqname <- paste(base.up,"T",base.down,sep = "")
                    seq96 <- append(seq96,seqname)
                }
            }
        }
    }
    names(triplet96) <- seq96
    
    ref64.type <- c()
    ref64.seq <- c()
    for(base.mid in bases){
        base.mid1 <- base.mid
        for(base.up in bases){
            for(base.down in bases){
                tri <- paste(base.up,base.mid,base.down,sep = "")
                ref64.seq <- append(ref64.seq,tri)
                if(base.mid == "G"){
                    base.mid1 <- "C"
                }
                else if(base.mid == "A"){
                    base.mid1 <- "T"
                }
                n <- paste(base.up,base.mid1,base.down,sep = "")
                ref64.type <- append(ref64.type,n)
            }
        }
    }
    names(ref64.type) <- ref64.seq
    ref64 <- ref64.type
    
    if(withinTumor){
        ref.list <- split(mutSigRef, mutSigRef$Branch_Tumor_ID)
    }else{
        ref.list <- split(mutSigRef, mutSigRef$Branch_ID)
    }
    
    sigsInput <- lapply(ref.list, countTripletBranch,
                        triplet96 = triplet96,
                        ref64 = ref64,
                        refBuild = refBuild,
                        CT = CT)
    sigsInput <- plyr::rbind.fill(sigsInput)
    rownames(sigsInput) <- names(ref.list)
    
    return(sigsInput)
}

countTripletBranch <- function(ref, triplet96, ref64, refBuild, CT){
    
    CpG = c("ACG", "CCG", "TCG", "GCG")

    context <- Biostrings::getSeq(get(refBuild),Rle(ref$chr),ref$pos-1, ref$pos+1) 
    context <- as.character(context)
    context1 <- context
    context <- ref64[context]
    
    mut.types <- paste(ref$ref, ref$alt, sep = ">")
    mut.types = gsub('G>T', 'C>A', mut.types)
    mut.types = gsub('G>C', 'C>G', mut.types)
    mut.types = gsub('G>A', 'C>T', mut.types)
    mut.types = gsub('A>T', 'T>A', mut.types)
    mut.types = gsub('A>G', 'T>C', mut.types)
    mut.types = gsub('A>C', 'T>G', mut.types)
    
    count.mat <- matrix(0, ncol = length(triplet96), nrow = 1)
    colnames(count.mat) <- as.character(triplet96)
    
    for(i in 1:length(context)){
        
        tris <- triplet96[which(names(triplet96) == context[i])]
        type <- mut.types[i]
        
        if(CT & type == "C>T"){
            if(context1[i] %in% CpG){
                type <- "C>T at CpG"
            }else{
                type <- "C>T other"
            }
        }
        
        pos <- which(grepl(type,tris))
        tri <- tris[pos]
        
        count.mat[,tri] <- count.mat[,tri] + 1
    }
    
    count.df <- as.data.frame(count.mat)
}

calSim <- function(sigsInput, sigsRef){
    
    m1 <- t(sigsInput)
    m2 <- t(sigsRef)
    n1 <- ncol(m1)
    n2 <- ncol(m2)
    cos_sim.mat <- matrix(nrow = n1,ncol = n2)
    rownames(cos_sim.mat) <- colnames(m1)
    colnames(cos_sim.mat) <- colnames(m2)
    
    for(i in 1:n1){
        x <- m1[,i]
        x <- x/sum(x)
        for( j in 1:n2){
            y <- m2[,j]
            s <- as.numeric(x %*% y / (sqrt(x %*% x) * sqrt(y %*% y)))
            cos_sim.mat[i,j] <- s
        }
    }
    
    return(as.data.frame(cos_sim.mat) )
    
}


fitSignature <- function(sigsInput, sigsRef){
    
    type.n <- ncol(sigsInput)
    sample.n <- nrow(sigsInput)
    sig.n <- nrow(sigsRef)
    
    ref <- t(as.matrix(sigsRef))
    
    weight.mat <- matrix(1, ncol = sample.n, nrow = sig.n)
    Reconstructed.mat <- matrix(1, ncol = sample.n, nrow = type.n)
    
    ## solve nonnegative least-squares constraints.
    
    for(i in 1:sample.n){
        type.count <- as.numeric(sigsInput[i,]) 
        lsq <- pracma::lsqnonneg(ref, type.count)
        weight.mat[,i] <- lsq$x
        Reconstructed.mat[,i] <- ref %*% as.matrix(lsq$x)
    }
    
    colnames(weight.mat) <- rownames(sigsInput)
    rownames(weight.mat) <- rownames(sigsRef)
    
    colnames(Reconstructed.mat) <- rownames(sigsInput)
    rownames(Reconstructed.mat) <- colnames(sigsInput)
    
    return(list(weight = t(weight.mat), Reconstructed = t(Reconstructed.mat)))
    
}

## Provide a summary data frame for signatures in different branches.
doMutSigSummary <- function(treeMSOutput, withinTumor){
   mutSigsOutput <- treeMSOutput$mutSigsOutput
   df.aetiology <- treeMSOutput$df.aetiology
   
   ## list of branch names
   ls.branchesName <- mutSigsOutput$branch
   ls.aeti <- c()
   
   if(!is.null(df.aetiology)){
       for(i in 1:nrow(mutSigsOutput)){
           signature <- as.character(mutSigsOutput$sig[i]) 
           aetiology <- as.character(df.aetiology[which(df.aetiology$sig == signature), ]$aeti)
           ls.aeti <- c(ls.aeti, aetiology)
       }
       mutSigsOutput$aeti <- ls.aeti  
   }
   
   if(withinTumor){
      ## order branch tumor type 
      branch_tumor_ids <- unique(mutSigsOutput$Branch_Tumor_ID) 
      public <- branch_tumor_ids[grep("Public", branch_tumor_ids)] 
      shared <- branch_tumor_ids[grep("Shared", branch_tumor_ids)] 
      private <- branch_tumor_ids[grep("Private", branch_tumor_ids)]
      type.level <- c(public, shared, private)
      mutSigsOutput$Branch_Tumor_ID <- factor(mutSigsOutput$Branch_Tumor_ID, levels = type.level)
      
      if(!is.null(df.aetiology)){
          colnames(mutSigsOutput) <- c("Branch_Tumor_ID", "Signature","Mutation_number", "Signature_weight", "Aetiology")
      }else{
          colnames(mutSigsOutput) <- c("Branch_Tumor_ID", "Signature","Mutation_number", "Signature_weight")
      }

      mutSigsOutput <- mutSigsOutput %>% 
         dplyr::arrange(Branch_Tumor_ID,Signature_weight)
   }
   else{
      ## rename column
       if(!is.null(df.aetiology)){
           colnames(mutSigsOutput) <- c("Branch","Signature","Mutation_number","Signature_weight","Branch_Tumor_ID", "Aetiology")
       }else{
           colnames(mutSigsOutput) <- c("Branch","Signature","Mutation_number","Signature_weight","Branch_Tumor_ID")
       }
      
      mutSigsOutput <- mutSigsOutput %>% 
         dplyr::select(-Branch_Tumor_ID) %>% 
         dplyr::arrange(plyr::desc(Branch), plyr::desc(Signature_weight))
   }
   mutSigsOutput <- mutSigsOutput %>% 
      dplyr::filter(Signature != "noMapSig")
   
   if(nrow(mutSigsOutput) == 0){
       mutSigsOutput <- NA
       message(paste0("Warning: There is no signature on ",treeMSOutput$patientID) )
   }
   
   return(mutSigsOutput)
}
