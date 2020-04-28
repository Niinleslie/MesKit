doTreeMutSig <- function(phyloTree,
                         geneList=NULL,
                         min.mut.count=15,
                         signaturesRef="cosmic_v2",
                         tri.counts.method = "default",
                         withinType = FALSE){
   
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
   
   
   
   # if(use.shiny){
   #     incProgress(amount=1)
   #     setProgress(message = paste('Generating ', "mutation signatures - ", patientID, sep=""))
   # }
   
   ## Select mutations in selected genes
   if(!is.null(geneList)){
      selectedGenes <- geneList
      mutSigRef <- mutSigRef  %>%
         dplyr::rowwise() %>%
         dplyr::filter(any(strsplit(Hugo_Symbol, ",|;")[[1]] %in% selectedGenes)) %>%
         as.data.frame()
   }
   
   ## get the mutational signature of the branch 
   mutSigsOutput <- data.frame()
   sig.spectrum <- data.frame()
   
   ## count 96 substitution typs in each branches
   sigsInput <- countTriplet(mutSigRef = mutSigRef,
                             withinType = withinType,
                             refBuild = refBuild,
                             patientID = patientID,
                             CT = FALSE)
   

   signatures.aetiology <- readRDS(file = system.file("extdata", "signatures.aetiology.rds", package = "MesKit")) 
   if (signaturesRef == "cosmic_v2") {
       sigsRef <- readRDS(file = system.file("extdata", "signatures.cosmic.rds", package = "MesKit"))
       df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v2$aetiology,
                                  sig = rownames(signatures.aetiology$cosmic_v2))
   } else if (signaturesRef == "nature2013") {
      sigsRef <- sigsRef <- readRDS(file = system.file("extdata", "signatures.nature2013.rds", package = "MesKit"))
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
   
   ## cal cos sim signature and branch 
   cos_sim.mat <- calSim(sigsInput = sigsInput[which(rowSums(sigsInput)!=0),], sigsRef = sigsRef)
   
   if(withinType){
      branchesName <- unique(mutSigRef$Branch_Tumor_Type) 
   }else{
       branchesName <- unique(mutSigRef$Branch_ID) 
   }
   
   signatures.aetiology <- readRDS(file = system.file("extdata", "signatures.aetiology.rds", package = "MesKit")) 
   for (branchCounter in length(branchesName):1){
      ## generate a single branch
      branchName <- branchesName[branchCounter]
      
      if(withinType){
         branch.mut.num <- length(mutSigRef[which(mutSigRef$Branch_Tumor_Type == branchName), 1])
      }else{
         branch.mut.num <-  length(mutSigRef[which(mutSigRef$Branch_ID == branchName), 1])
      }
      ## get the mutational signature of the branch
      if (branch.mut.num < min.mut.count){
         sigs.branch.name <- "noMapSig"
         sigs.prob <- 0
         message(paste0("Warnings: ",
                        "mutation number of Branch " , branchName,
                        " is less than the min.mut.count argument! ",
                        "This branch will be skipped")
         )
         # if (any(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$mut_id == "NoSigTag")) {
         #    mut.count <- 0
         # }
        if(withinType){
           mut.count <- length(mutSigRef[which(mutSigRef$Branch_Tumor_Type == branchName), 1])
        }
        else{
           mut.count <- length(mutSigRef[which(mutSigRef$Branch_ID == branchName), 1])
        }
      }else{
          ## get mutation signature
         fit <- fitSignature(sigsInput = sigsInput[branchName,], sigsRef = sigsRef)
         fit$Reconstructed <- fit$Reconstructed/sum(fit$Reconstructed)
         fit$contribution <- fit$contribution/sum(fit$contribution)
         sigsWhich <- fit$contribution
         
         # sigs.branch <- sigsWhich[["weights"]][which(sigsWhich[["weights"]] != 0)]
         sigs.branch <- sigsWhich[,which(sigsWhich > 0.1)]
         sigs.branch.name <- colnames(sigsWhich)[which(sigsWhich > 0.1)]
         sigs.branch.name <- gsub('[.]', ' ', sigs.branch.name)
         sigs.prob <- as.numeric(sigs.branch)
         
         if(withinType){
            mut.count <- length(mutSigRef[which(mutSigRef$Branch_Tumor_Type == branchName), 1])
         }
         else{
            mut.count <- length(mutSigRef[which(mutSigRef$Branch_ID == branchName), 1])
         }
         
         ## get spectrum of sigs output 
         sig.rec <- as.data.frame(fit$Reconstructed)
         sig.rec$Branch <- as.character(row.names(fit$Reconstructed)) 
         sig.branch.spectrum <- tidyr::pivot_longer(sig.rec ,
                                                   -Branch,
                                                   names_to = "group",
                                                   values_to = "sigs.prob") %>% 
                              as.data.frame()
         sig.branch.spectrum$value.type <- "Reconstructed"
         
         ## get Original fractions of 96 types
         sig.Original <- as.data.frame(sigsInput[branchName,]/sum(sigsInput[branchName,]))
         sig.Original$Branch <- as.character(row.names(fit$Reconstructed)) 
         sig.Original <- tidyr::pivot_longer(sig.Original ,
                                           -Branch,
                                           names_to = "group",
                                           values_to = "sigs.prob") %>% 
                       as.data.frame()
         sig.Original$value.type <- "Original"
         
         sig.branch.spectrum <- rbind(sig.branch.spectrum,sig.Original)
         if(!withinType){
             sig.branch.spectrum$Alias <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Alias))
         }
         
         ## sort signature by proportion
         sig.order <- order(sigs.branch,decreasing = T)
         sigs.branch <- sigs.branch[sig.order]
         sigs.branch.name <- sigs.branch.name[sig.order]
         ##  title
         t <- ''
         for(i in  1:length(sigs.branch)){
             if(i == 1){
                 t <- paste(sigs.branch.name[i], ": ", round(sigs.branch[i],3), sep = "")
                 next
             }
            t <- paste(t, " & ", sigs.branch.name[i], ": ", round(sigs.branch[i],3), sep = "")
         }
         sig.branch.spectrum$Signature <- t
         sig.spectrum <- rbind(sig.spectrum, sig.branch.spectrum)
      }
      
      mutSigsBranch <- data.frame(
         sig = sigs.branch.name, 
         sig.prob=format(round(sigs.prob, digits = 4)))
      
      if(withinType){
         mutSigsBranch$mut.count <- mut.count
         mutSigsBranch$Branch_Tumor_Type <- as.character(branchName)
         mutSigsBranch <- dplyr::select(mutSigsBranch,
                                        Branch_Tumor_Type,
                                        sig, 
                                        mut.count,
                                        sig.prob)
         
      }else{
         mutSigsBranch$branch <- c(branchName)
         mutSigsBranch$Alias <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Alias))
         mutSigsBranch$mut.count <- mut.count
         mutSigsBranch$Branch_Tumor_Type <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Branch_Tumor_Type))
         
         mutSigsBranch <- dplyr::select(mutSigsBranch,
                                        branch, Alias, sig, mut.count,
                                        sig.prob, Branch_Tumor_Type)
         
      }
      
      
      ## collect branches' mutataional signature information
      mutSigsOutput <- rbind(mutSigsOutput, mutSigsBranch)
   }
   mutSigsOutput$sig.prob <- as.numeric(as.vector(mutSigsOutput$sig.prob))
   
   ## calculation process(maybe could be replaced by lapply)
   if(nrow(sig.spectrum) > 0){
       if(withinType){
           sig.spectrum$Alias <- sig.spectrum$Branch
       }
       if(withinType){
           all.types <- unique(sig.spectrum$Alias) 
           public <- all.types[grep("Public", all.types)] 
           shared <- all.types[grep("Shared", all.types)] 
           private <- all.types[grep("Private", all.types)]
           type.level <- c(public, shared, private)
           sig.spectrum$Alias <- factor(sig.spectrum$Alias, levels = type.level)
           sig.spectrum$Patient_ID <- patientID
           colnames(sig.spectrum) <- c("Branch_Tumor_Type", "Group", "Mutation_Probability","value.type","Signature","Alias","Patient_ID")
           rownames(sig.spectrum) <- 1:nrow(sig.spectrum)
       }else{
           sig.spectrum$Patient_ID <- patientID
           colnames(sig.spectrum) <- c("Branch", "Group", "Mutation_Probability","value.type","Alias", "Signature","Patient_ID")
       }
   }
   
   cos_sim.mat <- as.matrix(cos_sim.mat)
   if(withinType){
       Alias <- unique(mutSigsOutput$Branch_Tumor_Type)
       names(Alias) <- unique(mutSigsOutput$Branch_Tumor_Type)
   }else{
       Alias <- unique(mutSigsOutput$Alias)
       names(Alias) <- unique(mutSigsOutput$branch)
   }
   rownames(cos_sim.mat) <- Alias[rownames(cos_sim.mat)] 
   
   message(paste0("Patient ", phyloTree@patientID, ": mutational signatures identified successfully!"))
   
   treeMSOutput <- list(
      sigsInput=sigsInput, 
      sig.spectrum=sig.spectrum,
      cos_sim.mat = cos_sim.mat,
      mutSigsOutput=mutSigsOutput, 
      df.aetiology=df.aetiology, 
      patientID = patientID)
   return(treeMSOutput)
}

countTriplet <- function(mutSigRef, withinType, refBuild, patientID, CT){
    
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
    
    if(withinType){
        ref.list <- split(mutSigRef, mutSigRef$Branch_Tumor_Type)
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
    
    contribution.mat <- matrix(1, ncol = sample.n, nrow = sig.n)
    Reconstructed.mat <- matrix(1, ncol = sample.n, nrow = type.n)
    
    ## solve nonnegative least-squares constraints.
    
    for(i in 1:sample.n){
        type.count <- as.numeric(sigsInput[i,]) 
        lsq <- pracma::lsqnonneg(ref, type.count)
        contribution.mat[,i] <- lsq$x
        Reconstructed.mat[,i] <- ref %*% as.matrix(lsq$x)
    }
    
    colnames(contribution.mat) <- rownames(sigsInput)
    rownames(contribution.mat) <- rownames(sigsRef)
    
    colnames(Reconstructed.mat) <- rownames(sigsInput)
    rownames(Reconstructed.mat) <- colnames(sigsInput)
    
    return(list(contribution = t(contribution.mat), Reconstructed = t(Reconstructed.mat)))
    
}

## Provide a summary data frame for signatures in different branches.
doMutSigSummary <- function(treeMSOutput, withinType){
   mutSigsOutput <- treeMSOutput$mutSigsOutput
   df.aetiology <- treeMSOutput$df.aetiology
   
   ## list of branch names
   ls.branchesName <- mutSigsOutput$branch
   ls.aeti <- c()
   
   for(i in 1:nrow(mutSigsOutput)){
      signature <- as.character(mutSigsOutput$sig[i]) 
      aetiology <- as.character(df.aetiology[which(df.aetiology$sig == signature), ]$aeti)
      ls.aeti <- c(ls.aeti, aetiology)
   }
   
   mutSigsOutput$aeti <- ls.aeti
   if(withinType){
      
      ## order branch tumor type 
      all.types <- unique(mutSigsOutput$Branch_Tumor_Type) 
      public <- all.types[grep("Public", all.types)] 
      shared <- all.types[grep("Shared", all.types)] 
      private <- all.types[grep("Private", all.types)]
      type.level <- c(public, shared, private)
      mutSigsOutput$Branch_Tumor_Type <- factor(mutSigsOutput$Branch_Tumor_Type, levels = type.level)
      
      colnames(mutSigsOutput) <- c("Branch_Tumor_Type", "Signature","Mutation_number", "Signature_weight", "Aetiology")
      
      mutSigsOutput <- mutSigsOutput %>% 
         dplyr::arrange(Branch_Tumor_Type,Signature_weight)
   }
   else{
      
      ## rename column
      colnames(mutSigsOutput) <- c("Branch", "Alias",  "Signature","Mutation_number","Signature_weight","Branch_Tumor_Type", "Aetiology")
      
      mutSigsOutput <- mutSigsOutput %>% 
         dplyr::select(-Branch_Tumor_Type) %>% 
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
