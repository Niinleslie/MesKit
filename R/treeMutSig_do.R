doTreeMutSig <- function(phyloTree,
                         geneList=NULL,
                         min.mut.count=15,
                         signaturesRef="cosmic_v2",
                         tri.counts.method = "default",
                         withinType = FALSE,
                         MTB = FALSE){
   
   refBuild <- phyloTree@refBuild
   ref.options = c('hg18', 'hg19', 'hg38')
   if(!refBuild %in% ref.options){
      stop("refBuild can only be either 'hg18', 'hg19' or 'hg38'")
   }else {
      refBuild <- paste("BSgenome.Hsapiens.UCSC.", refBuild, sep = "")
   }
   
   ## get branches information from phyloTree object
   mutBranches <- phyloTree@mut.branches
   
   patientID <- phyloTree@patientID
   branchesName <- names(mutBranches)
   branchesNameList <- strsplit(branchesName, split='∩')
   
   ## regain the data frame of all branches with Branch_ID
   mutSigRef <- data.frame()
   for (branchName in branchesName){
      mutBranch <- mutBranches[[branchName]]
      mutSigRef <- rbind(mutBranch, mutSigRef)
   }
   colnames(mutSigRef) <- c("Branch_ID","Branch_Tumor_Type","chr", "pos", "pos_end", 
                            "ref", "alt", "Hugo_Symbol", "mut_id", 
                            "alias")
   
   ## Select mutations in selected genes
   if(!is.null(geneList)){
      selectedGenes <- geneList
      mutSigRef <- mutSigRef  %>%
         dplyr::rowwise() %>%
         dplyr::filter(any(strsplit(Hugo_Symbol, ",|;")[[1]] %in% selectedGenes)) %>%
         as.data.frame()
      branchesNameList <- strsplit(unique(mutSigRef$Branch_ID), split='∩')
      branchesName <- unique(mutSigRef$Branch_ID)
      
   }
   
   ## get the mutational signature of the branch 
   mutSigsOutput <- data.frame()
   sig.product <- c()
   lsPicName <- c()
   sig.product <- data.frame()
   
   ## count 96 substitution typs in each branches
   sigsInput <- countTriplet(mutSigRef = mutSigRef,
                             withinType = withinType,
                             refBuild = refBuild)
   
   branchesName <- rownames(sigsInput)
   
   if(MTB){
      if(withinType){
         trunkName <- unique(mutSigRef[which(mutSigRef$alias == "Public"), ]$Branch_ID) 
      }
      else{
         trunkName <- unique(mutSigRef[which(mutSigRef$alias == "T"), ]$Branch_ID) 
      }
      treeMSOutput <- list(
         sigsInput = sigsInput, 
         trunkName = trunkName,
         patientID = patientID)
      return(treeMSOutput)
   }

   signatures.aetiology <- readRDS(file = system.file("extdata", "signatures.aetiology.rds", package = "MesKit")) 
   if (signaturesRef == "cosmic_v2") {
       sigsRef <- deconstructSigs::signatures.cosmic
       df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v2$aetiology,
                                  sig = rownames(signatures.aetiology$cosmic_v2))
   } else if (signaturesRef == "nature2013") {
      sigsRef <- deconstructSigs::signatures.nature2013
      df.aetiology <- data.frame(aeti = signatures.aetiology$nature2013$aetiology,
                                 sig = rownames(signatures.aetiology$nature2013))
   }else if(signaturesRef == "genome_cosmic_v3"){
       sigsRef <- deconstructSigs::signatures.genome.cosmic.v3.may2019
       df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v3$aetiology,
                                  sig = rownames(signatures.aetiology$cosmic_v3))
   }else if(signaturesRef == "exome_cosmic_v3"){
       sigsRef <- deconstructSigs::signatures.exome.cosmic.v3.may2019
       df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v3$aetiology,
                                  sig = rownames(signatures.aetiology$cosmic_v3))
   }
   
   ## cal cos sim signature and branch 
   cos_sim.mat <- calSim(sigsInput = sigsInput, sigsRef = sigsRef)
   
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
         if (any(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$mut_id == "NoSigTag")) {
            mut.count <- 0
         }
         else{
            if(withinType){
               mut.count <- length(mutSigRef[which(mutSigRef$Branch_Tumor_Type == branchName), 1])
            }
            else{
               mut.count <- length(mutSigRef[which(mutSigRef$Branch_ID == branchName), 1])
            }
         }
      }else{
          ## get mutation signature
         fit <- fitSignature(sigsInput = sigsInput[branchName,], sigsRef = sigsRef)
         fit$reconstruct <- fit$reconstruct/sum(fit$reconstruct)
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
         
         ## get product of sigs output 
         sig.rec <- as.data.frame(fit$reconstruct)
         sig.rec$Branch <- as.character(row.names(fit$reconstruct)) 
         sig.branch.product <- tidyr::pivot_longer(sig.rec ,
                                                   -Branch,
                                                   names_to = "group",
                                                   values_to = "sigs.prob") %>% 
            as.data.frame()
         
         if(!withinType){
             sig.branch.product$alias <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$alias))
         }
         
         ##  title
         t <- paste(sigs.branch.name[1],": ", round(sigs.branch[1],3), sep = "")
         if(length(sigs.branch) > 1){
            for(i in  2:length(sigs.branch)){
               t <- paste(t, " & ", sigs.branch.name[i], ": ", round(sigs.branch[i],3), sep = "")
            }
         }
         sig.branch.product$Signature <- t
         sig.product <- rbind(sig.product, sig.branch.product)
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
         mutSigsBranch$alias <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$alias))
         mutSigsBranch$mut.count <- mut.count
         mutSigsBranch$Branch_Tumor_Type <- as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Branch_Tumor_Type))
         
         mutSigsBranch <- dplyr::select(mutSigsBranch,
                                        branch, alias, sig, mut.count,
                                        sig.prob, Branch_Tumor_Type)
         
      }
      
      
      ## collect branches' mutataional signature information
      mutSigsOutput <- rbind(mutSigsOutput, mutSigsBranch)
   }
   
   if(withinType){
       sig.product$alias <- sig.product$Branch
   }

   mutSigsOutput$sig.prob <- as.numeric(as.vector(mutSigsOutput$sig.prob))

   ## calculation process(maybe could be replaced by lapply)
   if(nrow(sig.product) > 0){
       if(withinType){
           all.types <- unique(sig.product$alias) 
           public <- all.types[grep("Public", all.types)] 
           shared <- all.types[grep("Shared", all.types)] 
           private <- all.types[grep("Private", all.types)]
           type.level <- c(public, shared, private)
           sig.product$alias <- factor(sig.product$alias, levels = type.level)
           sig.product$Patient_ID <- patientID
           colnames(sig.product) <- c("Branch_Tumor_Type", "Group", "Mutation_Probability","Signature","Alias","Patient_ID")
       }else{
           sig.product$Patient_ID <- patientID
           colnames(sig.product) <- c("Branch", "Group", "Mutation_Probability", "Alias", "Signature","Patient_ID")
       }
   }
   
   rownames(sig.product) <- 1:nrow(sig.product)
   
   cos_sim.mat <- as.matrix(cos_sim.mat)
   if(withinType){
       alias <- unique(mutSigsOutput$Branch_Tumor_Type)
       names(alias) <- unique(mutSigsOutput$Branch_Tumor_Type)
   }else{
       alias <- unique(mutSigsOutput$alias)
       names(alias) <- unique(mutSigsOutput$branch)
   }
   rownames(cos_sim.mat) <- alias[rownames(cos_sim.mat)] 
   
   message(paste0("Patient ", phyloTree@patientID, ": mutational signatures identified successfully!"))
   
   treeMSOutput <- list(
      sigsInput=sigsInput, 
      sig.product=sig.product,
      cos_sim.mat = cos_sim.mat,
      mutSigsOutput=mutSigsOutput, 
      df.aetiology=df.aetiology, 
      patientID = patientID)
   return(treeMSOutput)
}

countTriplet <- function(mutSigRef, withinType, refBuild){
    
    mutSigRef <- mutSigRef[mutSigRef$mut_id!="NoSigTag",]
    
    bases <- c("A","C","G","T")
    types <- c("C>A","C>G","C>T","T>A","T>C","T>G")
    
    ## 96 substitution types
    triplet96 <- c()
    seq96 <- c()
    for(type in types){
        for(base.up in bases){
            for(base.down in bases){
                tri <- paste(base.up,"[",type,"]",base.down,sep = "")
                triplet96 <- append(triplet96,tri)
                if(type %in% types[1:3]){
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
                        refBuild = refBuild)
    sigsInput <- plyr::rbind.fill(sigsInput)
    rownames(sigsInput) <- names(ref.list)
    
    return(sigsInput)
}

countTripletBranch <- function(ref, triplet96, ref64, refBuild){
    
    types <- c("C>A","C>G","C>T","T>A","T>C","T>G")
    
    
    context <- Biostrings::getSeq(get(refBuild),Rle(ref$chr),ref$pos-1, ref$pos+1) 
    context <- as.character(context)
    
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
    reconstruct.mat <- matrix(1, ncol = sample.n, nrow = type.n)
    
    ## solve nonnegative least-squares constraints.
    
    for(i in 1:sample.n){
        type.count <- as.numeric(sigsInput[i,]) 
        lsq <- pracma::lsqnonneg(ref, type.count)
        contribution.mat[,i] <- lsq$x
        reconstruct.mat[,i] <- ref %*% as.matrix(lsq$x)
    }
    
    colnames(contribution.mat) <- rownames(sigsInput)
    rownames(contribution.mat) <- rownames(sigsRef)
    
    colnames(reconstruct.mat) <- rownames(sigsInput)
    rownames(reconstruct.mat) <- colnames(sigsInput)
    
    return(list(contribution = t(contribution.mat), reconstruct = t(reconstruct.mat)))
    
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


# ## Visualize the mutational signatures of trunk/branches in a phylogenetic tree.
# doPlotMutSig <- function(tree.mutSig) {
#     
#     sig.product <- tree.mutSig$sig.product
#     
#     if(nrow(sig.product) == 0){
#         return(NA)
#     }
#     
#     sig.product$Signature <- as.character(sig.product$Signature)
#     ls.mutationType <-as.character(sig.product$Group)
#     ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
#     
#     ## generate Mutation Type for every column
#     sig.product$Group <- as.character(sig.product$Group)
#     for (mutationGroup in ls.mutationGroup) {
#         sig.product$Group[which(grepl(mutationGroup, sig.product$Group))] <- mutationGroup
#     }
#     
#     ## specific the label order of x axis
#     orderlist <- c(ls.mutationType)
#     sig.product$Type <- factor(ls.mutationType, levels = unique(orderlist) )
#     sig.product$Group <- factor(sig.product$Group, levels = ls.mutationGroup)
#     
#     if(nrow(sig.product) == 0){
#         warning(paste0("Patient ", tree.mutSig$patientID, ": There is no enough eligible mutations can be used."))
#         return(NA)
#     }
#     if("Branch" %in% colnames(sig.product)){
#         sig.product$Alias <- factor(sig.product$Alias, levels = unique(sig.product$Alias))
#     }
#     
#     # sig.product <- dplyr::distinct(sig.product, Branch, .keep_all = TRUE)
#     
#     group.colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF",
#                       "#3C5488FF", "#F39B7FFF", "#8491B4FF")
#     
#     # set the limit of y axis
#     yMax <- max(max(sig.product$Mutation_Probability), round(max(sig.product$Mutation_Probability) + 0.1, 1))
# 
#     pic <- ggplot(sig.product, 
#                   aes(x=Type, y=Mutation_Probability, group=Group, fill=Group)
#     ) + 
#         geom_bar(stat="identity") +
#         theme(panel.grid=element_blank(), 
#               panel.border=element_blank(), 
#               panel.background = element_blank(), 
#               legend.position='none', 
#               # axis.text.x=element_text(size=3, angle = 45, hjust = 1, vjust = 1), 
#               plot.title = element_text(size = 13, face = "bold", hjust = 0.5,vjust = 0),
#               axis.text.x= element_blank(), 
#               axis.ticks.x=element_blank(),
#               axis.line.y = element_blank(),
#               axis.ticks.length.y = unit(0.2, "cm"),
#               axis.text.y=element_text(size=6, color = "black"),
#               strip.background = element_blank(),
#               strip.text = element_blank(),
#         ) +
#         annotate(geom = "segment", x=-1, xend = -1, y = 0, 
#                  yend = yMax, size = 0.6)+
#         
#         ## side bar
#         geom_rect(aes(xmin = 96.5, xmax = 99, ymin = 0, ymax = yMax), fill = "#C0C0C0",alpha = 0.06) +
#         geom_text(data = dplyr::distinct(sig.product,Alias,.keep_all = TRUE),
#                   aes(x=98, y= yMax/2, label = Alias), 
#                   angle = 270, color = "black", size = 3, fontface = "plain") + 
#         
#         
#         ## background colors
#         geom_rect(aes(xmin=0, xmax=16.5, ymin=0, ymax=yMax),
#                   fill="#fce7e4", alpha=0.15) + 
#         geom_rect(aes(xmin=16.5, xmax=32.5, ymin=0, ymax=yMax),
#                   fill="#ecf8fa", alpha=0.25) + 
#         geom_rect(aes(xmin=32.5, xmax=48.5, ymin=0, ymax= yMax),
#                   fill="#dbfff9", alpha=0.05) + 
#         geom_rect(aes(xmin=48.5, xmax=64.5, ymin=0, ymax= yMax),
#                   fill="#e4e8f3", alpha=0.08) + 
#         geom_rect(aes(xmin=64.5, xmax=80.5, ymin=0, ymax= yMax),
#                   fill="#fdefeb", alpha=0.15) + 
#         geom_rect(aes(xmin=80.5, xmax=96.5, ymin=0, ymax= yMax),
#                   fill="#e5e8ef", alpha=0.1) + 
#         ## barplot
#         geom_bar(stat="identity") + 
#         ## combine different results
#         facet_grid(Alias ~ .) + 
#         ## color setting
#         scale_fill_manual(values=group.colors) +
#         
#         ## axis setting
#         xlab("Mutational type") + 
#         ylab("Mutation probability") + 
#         # ggtitle(paste0("Mutational signatures of ",tree.mutSig$patientID,"'s phylogenetic tree ") ) +
#         # scale_x_discrete(breaks = c(10,27,43,59,75,91),
#         #                  labels = c("C>A","C>G","C>T","T>A","T>C","T>G"))+
#         scale_y_continuous(limits=c(0, yMax+ 0.03), 
#                            breaks=seq(0, yMax, 0.1)) +
#         
#         ## signature notes and text parts
#         geom_text(data = dplyr::distinct(sig.product,Alias,.keep_all = TRUE) , 
#                   aes(x=0, y=yMax, label= Signature), 
#                   hjust = 0, vjust = 1.5, colour="#2B2B2B", size=3.5) + 
#         
#         ## x axis bar
#         geom_rect(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]),
#                   aes(xmin=0, xmax=16.5, ymin=yMax, ymax=yMax + 0.03),
#                   fill=group.colors[1], alpha=1) + 
#         geom_rect(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]),
#                   aes(xmin=16.5, xmax=32.5, ymin=yMax, ymax=yMax + 0.03),
#                   fill=group.colors[2], alpha=0.25) + 
#         geom_rect(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]),
#                   aes(xmin=32.5, xmax=48.5, ymin=yMax, ymax=yMax + 0.03),
#                   fill=group.colors[3], alpha=0.05) + 
#         geom_rect(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]),
#                   aes(xmin=48.5, xmax=64.5, ymin=yMax, ymax=yMax + 0.03),
#                   fill=group.colors[4], alpha=0.08) + 
#         geom_rect(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]),
#                   aes(xmin=64.5, xmax=80.5, ymin=yMax, ymax=yMax + 0.03),
#                   fill=group.colors[5], alpha=0.15) + 
#         geom_rect(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]),
#                   aes(xmin=80.5, xmax=96.5, ymin=yMax, ymax=yMax + 0.03),
#                   fill=group.colors[6], alpha=0.1) +
#         
#         # Mutational Type Labels
#         geom_text(data = subset(sig.product, Alias == levels(sig.product$Alias)[1] & Group == "C>A")[1,],
#                   aes(x = 10, y= yMax + 0.02, label = "C>A"),size = 3,hjust = 1) +
#         geom_text(data = subset(sig.product, Alias == levels(sig.product$Alias)[1] & Group == "C>G")[1,],
#                   aes(x = 26, y= yMax + 0.02, label = "C>G"),size = 3,hjust = 1) +
#         geom_text(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]& Group == "C>T")[1,],
#                   aes(x = 42, y= yMax + 0.02, label = "C>T"),size = 3,hjust = 1) +
#         geom_text(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]& Group == "T>A")[1,],
#                   aes(x = 58, y= yMax + 0.02, label = "T>A"),size = 3,hjust = 1) +
#         geom_text(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]& Group == "T>C")[1,],
#                   aes(x = 74, y= yMax + 0.02, label = "T>C"),size = 3,hjust = 1) +
#         geom_text(data = subset(sig.product, Alias == levels(sig.product$Alias)[1]& Group == "T>G")[1,],
#                   aes(x = 90, y= yMax + 0.02, label = "T>G"),size = 3,hjust = 1)
#     #pic <- pic + annotate(geom = "segment", x=10, xend = 10, y = 0, 
#     #yend = yMax, size = 2)
#     #yend = yMax, size = 2)
#     
#     message("Mutational signature plot generation done!")
#     return(pic)
# }