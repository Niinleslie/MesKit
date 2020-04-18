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
   ## deconstructSigs
   if(withinType){
      sigsInput <- suppressWarnings(
         deconstructSigs::mut.to.sigs.input(mut.ref=mutSigRef, 
                                            sample.id="Branch_Tumor_Type", 
                                            chr="chr", 
                                            pos="pos", 
                                            ref="ref", 
                                            alt="alt",
                                            bsg=get(refBuild)))
      branchesName <- unique(mutSigRef$Branch_Tumor_Type)
      
   }else{
      sigsInput <- suppressWarnings(
         deconstructSigs::mut.to.sigs.input(mut.ref=mutSigRef, 
                                            sample.id="Branch_ID", 
                                            chr="chr", 
                                            pos="pos", 
                                            ref="ref", 
                                            alt="alt",
                                            bsg=get(refBuild)))
   }
   
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
   
   #signatures.aetiology <- readRDS(file = system.file("extdata", "signatures.aetiology.rds", package = "MesKit")) 
   signatures.aetiology <- readRDS(file = system.file("extdata", "signatures.aetiology.rds", package = "MesKit")) 
   for (branchCounter in length(branchesName):1){
      ## generate a single branch
      # branch <- Filter(Negate(is.na), 
      #                  branchesNameList[[branchCounter]])
      branchName <- branchesName[branchCounter]
      
      if(withinType){
         branch.mut.num <- length(mutSigRef[which(mutSigRef$Branch_Tumor_Type == branchName), 1])
      }
      else{
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
         if (signaturesRef == "cosmic_v2") {
            sigsWhich <- deconstructSigs::whichSignatures(tumor.ref=sigsInput, 
                                                          signatures.ref=deconstructSigs::signatures.cosmic, 
                                                          sample.id=branchName,
                                                          contexts.needed=TRUE,
                                                          tri.counts.method = tri.counts.method)
         } else if (signaturesRef == "nature2013") {
            sigsWhich <- deconstructSigs::whichSignatures(tumor.ref=sigsInput, 
                                                          signatures.ref=deconstructSigs::signatures.nature2013, 
                                                          sample.id=branchName,
                                                          contexts.needed=TRUE,
                                                          tri.counts.method = tri.counts.method)
         }
         else if(signaturesRef == "genome_cosmic_v3"){
            sigsWhich <- deconstructSigs::whichSignatures(tumor.ref=sigsInput, 
                                                          signatures.ref= deconstructSigs::signatures.genome.cosmic.v3.may2019, 
                                                          sample.id=branchName,
                                                          contexts.needed=TRUE,
                                                          tri.counts.method = tri.counts.method)
         }
         else if(signaturesRef == "exome_cosmic_v3"){
            sigsWhich <- deconstructSigs::whichSignatures(tumor.ref=sigsInput, 
                                                          signatures.ref= deconstructSigs::signatures.exome.cosmic.v3.may2019, 
                                                          sample.id=branchName,
                                                          contexts.needed=TRUE,
                                                          tri.counts.method = tri.counts.method)
         }
         
         sigs.branch <- sigsWhich[["weights"]][which(sigsWhich[["weights"]] != 0)]
         sigs.branch.name <- colnames(sigs.branch)
         sigs.branch.name <- gsub('[.]', ' ', sigs.branch.name)
         sigs.prob <- as.numeric(sigs.branch)
         
         if(withinType){
            mut.count <- length(mutSigRef[which(mutSigRef$Branch_Tumor_Type == branchName), 1])
         }
         else{
            mut.count <- length(mutSigRef[which(mutSigRef$Branch_ID == branchName), 1])
         }
         
         ## get product of sigs output 
         # print(sigsWhich$product)
         sigsWhich$product <- as.data.frame(sigsWhich$product)
         sigsWhich$product$Branch <- as.character(row.names(sigsWhich$product)) 
         sig.branch.product <- tidyr::pivot_longer(sigsWhich$product,
                                                   -Branch,
                                                   names_to = "group",
                                                   values_to = "sigs.prob") %>% 
            as.data.frame()
         # print(sigs.branch.product)
         # sig.branch.product <- reshape2::melt(sigsWhich$product)
         # print(sigsWhich$product)
         # print(sig.branch.product)
         
         if(withinType){
            sig.branch.product$alias <- branchName
         }
         else{
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
         ## get mutational signature with max weight
         # sigsMax <- sigsWhich[["weights"]][which.max(sigsWhich[["weights"]])]
         # sigsMaxName <- colnames(sigsMax)
         # sigsMaxName <- gsub('[.]', ' ', sigsMaxName)
         # sigsMaxProb <- sigsMax[,1]
         
      }
      
      # mutSigsBranch <- data.frame(
      #     branch=c(branchName), 
      #     alias=as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$alias)), 
      #     sig=sigsMaxName, 
      #     mut.count=mut.count, 
      #     sig.prob=format(round(sigsMaxProb, digits = 3), nsmall = 3),
      #     Branch_Tumor_Type = as.character(unique(mutSigRef[which(mutSigRef$Branch_ID == branchName), ]$Branch_Tumor_Type)))
      # 
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
   
   mutSigsOutput$sig.prob <- as.numeric(as.vector(mutSigsOutput$sig.prob))
   ## calculation process(maybe could be replaced by lapply)
   if (signaturesRef =="cosmic_v2") {
      ## Aetiology from https://cancer.sanger.ac.uk/cosmic/signatures_v2 emm actually the additional feature may matter
      df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v2$aetiology,
                                 sig = rownames(signatures.aetiology$cosmic_v2))
   } else if (signaturesRef =="nature2013") {
      ## Aetiology from https://www.nature.com/articles/nature12477#s1
      df.aetiology <- data.frame(aeti = signatures.aetiology$nature2013$aetiology,
                                 sig = rownames(signatures.aetiology$nature2013))
      
   }
   else if(signaturesRef == "genome_cosmic_v3"){
      df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v3$aetiology,
                                 sig = rownames(signatures.aetiology$cosmic_v3))
   }
   else if(signaturesRef == "exome_cosmic_v3"){
      df.aetiology <- data.frame(aeti = signatures.aetiology$cosmic_v3$aetiology,
                                 sig = rownames(signatures.aetiology$cosmic_v3))
   }
   
   if(nrow(sig.product) > 0){
       if(withinType){
           all.types <- unique(sig.product$alias) 
           public <- all.types[grep("Public", all.types)] 
           shared <- all.types[grep("Shared", all.types)] 
           private <- all.types[grep("Private", all.types)]
           type.level <- c(public, shared, private)
           sig.product$alias <- factor(sig.product$alias, levels = type.level)
           colnames(sig.product) <- c("Branch_Tumor_Type", "Group", "Mutation_Probability", "Alias", "Signature")
       }else{
           colnames(sig.product) <- c("Branch", "Group", "Mutation_Probability", "Alias", "Signature")
       }
   }
   
   message(paste0("Sample ", phyloTree@patientID, ": mutational signatures identified successfully!"))
   treeMSOutput <- list(
      sigsInput=sigsInput, 
      sig.product=sig.product,
      mutSigsOutput=mutSigsOutput, 
      df.aetiology=df.aetiology, 
      patientID = patientID)
   return(treeMSOutput)
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
   # for (branch in ls.branchesName) {
   #     signature <- as.character(mutSigsOutput[which(mutSigsOutput$branch == branch), ]$sig)
   #     aetiology <- as.character(df.aetiology[which(df.aetiology$sig == signature), ]$aeti)
   #     ls.aeti <- c(ls.aeti, aetiology)
   # }
   
   ## rearrange the order of columns
   # mutSigsOutput <- data.frame(branch=mutSigsOutput$branch,
   #                             alias=mutSigsOutput$alias, 
   #                             mut.count=mutSigsOutput$mut.count, 
   #                             sig=mutSigsOutput$sig, 
   #                             sig.prob=mutSigsOutput$sig.prob, 
   #                             Branch_Tumor_Type = mutSigsOutput$Branch_Tumor_Type,
   #                             aeti=ls.aeti)
   
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


## Visualize the mutational signatures of trunk/branches in a phylogenetic tree.
doPlotMutSig <- function(tree.mutSig, withinType) {
   sigsInput <- tree.mutSig$sigsInput
   mutSigsOutput <- tree.mutSig$mutSigsOutput
   df.aetiology <- tree.mutSig$df.aetiology
   sig.product <- tree.mutSig$sig.product
   
   if(nrow(sig.product) == 0){
       return(NA)
   }
   
   sig.product$Signature <- as.character(sig.product$Signature)
   ls.mutationType <-as.character(sig.product$Group)
   ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
   
   ## generate Mutation Type for every column
   sig.product$Group <- as.character(sig.product$Group)
   for (mutationGroup in ls.mutationGroup) {
      sig.product$Group[which(grepl(mutationGroup, sig.product$Group))] <- mutationGroup
   }
   
   ## specific the label order of x axis
   orderlist <- c(ls.mutationType)
   sig.product$Type <- factor(ls.mutationType, levels = unique(orderlist) )
   
   if(nrow(sig.product) == 0){
      warning(paste0("Patient ", tree.mutSig$patientID, ": There is no enough eligible mutations can be used."))
      return(NA)
   }
   # sig.product <- dplyr::distinct(sig.product, Branch, .keep_all = TRUE)
   
   CA <- grid::textGrob(expression(bold("C > A")),
                        gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
   CG <- grid::textGrob(expression(bold("C > G")),
                        gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
   CT <- grid::textGrob(expression(bold("C > T")),
                        gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
   TA <- grid::textGrob(expression(bold("T > A")),
                        gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
   TC <- grid::textGrob(expression(bold("T > C")),
                        gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
   TG <- grid::textGrob(expression(bold("T > G")),
                        gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
   
   group.colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF",
                     "#3C5488FF", "#F39B7FFF", "#8491B4FF")
   
   # set the limit of y axis
   yMax <- max(max(sig.product$Mutation_Probability), round(max(sig.product$Mutation_Probability) + 0.1, 1))
   
   pic <- ggplot(sig.product, 
                 aes(x=Type, y=Mutation_Probability, group=Group, fill=Group)
   ) + 
      geom_bar(stat="identity") +
      theme(panel.grid=element_blank(), 
            panel.border=element_blank(), 
            panel.background = element_blank(), 
            legend.position='none', 
            # axis.text.x=element_text(size=3, angle = 45, hjust = 1, vjust = 1), 
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5,vjust = 0),
            axis.text.x= element_blank(), 
            axis.ticks.x=element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.length.y = unit(0.2, "cm"),
            axis.text.y=element_text(size=6, color = "black"),
            strip.background = element_blank(),
            strip.text = element_blank(),
      ) +
      annotate(geom = "segment", x=-1, xend = -1, y = 0, 
               yend = yMax, size = 0.6)+
      
      ## side bar
      geom_rect(aes(xmin = 96.7, xmax = 99, ymin = 0, ymax = yMax), fill = "#C0C0C0",alpha = 0.06) +
      geom_text(data = dplyr::distinct(sig.product,Alias,.keep_all = TRUE),
                aes(x=98, y= yMax/2, label = Alias), 
                angle = 270, color = "black", size = 3, fontface = "plain") + 
      
      
      ## background colors
      geom_rect(aes(xmin=0, xmax=16.5, ymin=0, ymax=yMax),
                fill="#fce7e4", alpha=0.15) + 
      geom_rect(aes(xmin=16.5, xmax=32.5, ymin=0, ymax=yMax),
                fill="#ecf8fa", alpha=0.25) + 
      geom_rect(aes(xmin=32.5, xmax=48.5, ymin=0, ymax= yMax),
                fill="#dbfff9", alpha=0.05) + 
      geom_rect(aes(xmin=48.5, xmax=64.5, ymin=0, ymax= yMax),
                fill="#e4e8f3", alpha=0.08) + 
      geom_rect(aes(xmin=64.5, xmax=80.5, ymin=0, ymax= yMax),
                fill="#fdefeb", alpha=0.15) + 
      geom_rect(aes(xmin=80.5, xmax=96.5, ymin=0, ymax= yMax),
                fill="#e5e8ef", alpha=0.1) + 
      ## barplot
      geom_bar(stat="identity") + 
      ## combine different results
      facet_grid(Alias ~ .) + 
      ## color setting
      scale_fill_manual(values=group.colors) +
      
      ## axis setting
      xlab("Mutational type") + 
      ylab("Mutation probability") + 
      ggtitle(paste0("Mutational signatures of ",tree.mutSig$patientID,"'s phylogenetic tree ") )+
      # scale_x_discrete(breaks = c(10,27,43,59,75,91),
      #                  labels = c("C>A","C>G","C>T","T>A","T>C","T>G"))+
      scale_y_continuous(limits=c(-0.03, yMax), 
                         breaks=seq(0, yMax, 0.1)) +
      
      ## signature notes and text parts
      geom_text(data = dplyr::distinct(sig.product,Alias,.keep_all = TRUE) , 
                aes(x=0, y=yMax, label= Signature), 
                hjust = 0, vjust = 1.5, colour="#2B2B2B", size=3.5) + 
      
      # Mutational Type Labels
      annotation_custom(grob = CA,  xmin = 10, xmax = 10, ymin = -0.065, ymax = -0) +
      annotation_custom(grob = CG,  xmin = 27, xmax = 27, ymin = -0.065, ymax = -0) +
      annotation_custom(grob = CT,  xmin = 43, xmax = 43, ymin = -0.065, ymax = -0) +
      annotation_custom(grob = TA,  xmin = 59, xmax = 59, ymin = -0.065, ymax = -0) +
      annotation_custom(grob = TC,  xmin = 75, xmax = 75, ymin = -0.065, ymax = -0) +
      annotation_custom(grob = TG,  xmin = 91, xmax = 91, ymin = -0.065, ymax = -0) +
      ## x axis bar
      geom_rect(aes(xmin=0, xmax=16.405, ymin=-0.01, ymax=-0.005),
                fill=group.colors[1], alpha=1) + 
      geom_rect(aes(xmin=16.595, xmax=32.405, ymin=-0.01, ymax=-0.005),
                fill=group.colors[2], alpha=0.25) + 
      geom_rect(aes(xmin=32.595, xmax=48.405, ymin=-0.01, ymax=-0.005),
                fill=group.colors[3], alpha=0.05) + 
      geom_rect(aes(xmin=48.595, xmax=64.405, ymin=-0.01, ymax=-0.005),
                fill=group.colors[4], alpha=0.08) + 
      geom_rect(aes(xmin=64.595, xmax=80.405, ymin=-0.01, ymax=-0.005),
                fill=group.colors[5], alpha=0.15) + 
      geom_rect(aes(xmin=80.595, xmax=96.5, ymin=-0.01, ymax=-0.005),
                fill=group.colors[6], alpha=0.1)
   message("Mutational signature plot generation done!")
   return(pic)
}

# ## Visualize the mutational signatures of trunk/branches in a phylogenetic tree.
# doPlotMutSig <- function(tree.mutSig) {
#     sigsInput <- tree.mutSig$sigsInput
#     mutSigsOutput <- tree.mutSig$mutSigsOutput
#     df.aetiology <- tree.mutSig$df.aetiology
#     sig.product <- tree.mutSig$sig.product
#     
#     ## calculate the Mutation Probability
#     sigsInputSum <- as.data.frame(apply(sigsInput, 1, function(x) sum(x)))
#     sigsInputTrans <- as.data.frame(t(sigsInput))
#     ls.branchesName <- as.character(rownames(sigsInput))
#     ls.mutationType <-as.character(rownames(sigsInputTrans))
#     ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
#     sigsInputBranches <- data.frame()
#     
#     ## calculation process(maybe could be replaced by lapply)
#     for (branch in ls.branchesName) {
#         sigsInputTrans[branch] <- sigsInputTrans[branch]/sigsInputSum[branch, ]
#         signature <- as.character(mutSigsOutput[which(mutSigsOutput$branch == branch), ]$sig)
#         alias <- as.character(mutSigsOutput[which(mutSigsOutput$branch == branch), ]$alias)
#         sigsWeight <- mutSigsOutput[which(mutSigsOutput$branch == branch), ]$sig.prob
#         aetiology <- as.character(df.aetiology[which(df.aetiology$sig == signature), ]$aeti)
#         sigsInputBranch <- data.frame(sigsInputTrans[branch], 
#                                       rep(branch, length(sigsInputTrans[branch])), 
#                                       rep(alias, length(sigsInputTrans[branch])), 
#                                       rep(signature, length(sigsInputTrans[branch])), 
#                                       rep(sigsWeight, length(sigsInputTrans[branch])), 
#                                       rep(aetiology , length(sigsInputTrans[branch])), 
#                                       stringsAsFactors = FALSE)
#         colnames(sigsInputBranch) <- c("Mutation_Probability", "Branch", "Alias", "Signature", "SigsWeight", "Aetiology")
#         sigsInputBranches <- rbind(sigsInputBranches, sigsInputBranch)
#     }
#     
#     df.sigsInputTrans <- data.frame(Mutational_Type=ls.mutationType, 
#                                     Group=ls.mutationType, 
#                                     sigsInputBranches, 
#                                     stringsAsFactors = FALSE)
#     
#     ## generate Mutation Type for every column
#     for (mutationGroup in ls.mutationGroup) {
#         df.sigsInputTrans$Group[which(grepl(mutationGroup, df.sigsInputTrans$Group))] <- mutationGroup
#     }
#     
#     ## specific the label order of x axis
#     orderlist <- c(ls.mutationType)
#     df.sigsInputTrans <- transform(df.sigsInputTrans, Mutational_Type = factor(Mutational_Type, levels = orderlist))
#     df.sigsInputTrans <- df.sigsInputTrans[which(df.sigsInputTrans$Signature != "noMapSig"), ]
#     if(nrow(df.sigsInputTrans) == 0){
#         warning(paste0("Patient ", tree.mutSig$patientID, ": There is no enough eligible mutations can be used."))
#         return(NA)
#     }
#     df.sigsInputText <- dplyr::distinct(df.sigsInputTrans, Branch, .keep_all = TRUE)
#     
#     CA <- grid::textGrob(expression(bold("C > A")),
#                          gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
#     CG <- grid::textGrob(expression(bold("C > G")),
#                          gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
#     CT <- grid::textGrob(expression(bold("C > T")),
#                          gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
#     TA <- grid::textGrob(expression(bold("T > A")),
#                          gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
#     TC <- grid::textGrob(expression(bold("T > C")),
#                          gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
#     TG <- grid::textGrob(expression(bold("T > G")),
#                          gp=grid::gpar(fontsize=7, fontface="bold"), vjust=0,hjust=1)
#     
#     group.colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF",
#                       "#3C5488FF", "#F39B7FFF", "#8491B4FF")
#     pic <- ggplot(df.sigsInputTrans, 
#                   aes(x=Mutational_Type, y=Mutation_Probability, group=Group, fill=Group)
#     ) + 
#         geom_bar(stat="identity")+ 
#         theme(panel.grid=element_blank(), 
#               panel.border=element_blank(), 
#               panel.background = element_blank(), 
#               legend.position='none', 
#               # axis.text.x=element_text(size=3, angle = 45, hjust = 1, vjust = 1), 
#               plot.title = element_text(size = 13, face = "bold", hjust = 0.5,vjust = 0),
#               axis.text.x=element_blank(), 
#               axis.ticks.x=element_blank(),
#               axis.line.y = element_blank(),
#               axis.ticks.length.y = unit(0.2, "cm"),
#               axis.text.y=element_text(size=6, color = "black")
#         ) +
#         annotate(geom = "segment", x=-1, xend = -1, y = 0, yend = 0.2, size = 0.6)+
#         ## background colors
#         geom_rect(aes(xmin=0, xmax=16.5, ymin=0, ymax=Inf),
#                   fill="#fce7e4", alpha=0.15) + 
#         geom_rect(aes(xmin=16.5, xmax=32.5, ymin=0, ymax=Inf),
#                   fill="#ecf8fa", alpha=0.25) + 
#         geom_rect(aes(xmin=32.5, xmax=48.5, ymin=0, ymax=Inf),
#                   fill="#dbfff9", alpha=0.05) + 
#         geom_rect(aes(xmin=48.5, xmax=64.5, ymin=0, ymax=Inf),
#                   fill="#e4e8f3", alpha=0.08) + 
#         geom_rect(aes(xmin=64.5, xmax=80.5, ymin=0, ymax=Inf),
#                   fill="#fdefeb", alpha=0.15) + 
#         geom_rect(aes(xmin=80.5, xmax=96.5, ymin=0, ymax=Inf),
#                   fill="#e5e8ef", alpha=0.1) + 
#         ## barplot
#         geom_bar(stat="identity") + 
#         ## combine different results
#         facet_grid(Alias ~ .) + 
#         ## color setting
#         scale_fill_manual(values=group.colors) + 
#         ## axis setting
#         xlab("Mutational type") + 
#         ylab("Mutation probability") + 
#         ggtitle(paste0("Mutational signatures of ",tree.mutSig$patientID,"'s phylogenetic tree ") )+
#         scale_y_continuous(limits=c(-0.03, 0.2), breaks=seq(0, 0.2, 0.1)) + 
#         ## signature notes and text parts
#         geom_text(data = df.sigsInputText, 
#                   aes(x=-Inf, y=Inf, label=paste(Signature, ": ", 
#                                                  round(as.numeric(levels(SigsWeight)[SigsWeight]), 3), 
#                                                  "    ", 
#                                                  "Aetiology: ", Aetiology, sep="")), 
#                   hjust = -0.02, vjust = 1.5, colour="#2B2B2B", fontface = "bold", size=3) + 
#         ## Mutational Type Labels
#         annotation_custom(grob = CA,  xmin = 10, xmax = 10, ymin = -0.065, ymax = -0) + 
#         annotation_custom(grob = CG,  xmin = 27, xmax = 27, ymin = -0.065, ymax = -0) + 
#         annotation_custom(grob = CT,  xmin = 43, xmax = 43, ymin = -0.065, ymax = -0) + 
#         annotation_custom(grob = TA,  xmin = 59, xmax = 59, ymin = -0.065, ymax = -0) + 
#         annotation_custom(grob = TC,  xmin = 75, xmax = 75, ymin = -0.065, ymax = -0) + 
#         annotation_custom(grob = TG,  xmin = 91, xmax = 91, ymin = -0.065, ymax = -0) + 
#         ## x axis bar
#         geom_rect(aes(xmin=0, xmax=16.405, ymin=-0.01, ymax=-0.005),
#                   fill=group.colors[1], alpha=1) + 
#         geom_rect(aes(xmin=16.595, xmax=32.405, ymin=-0.01, ymax=-0.005),
#                   fill=group.colors[2], alpha=0.25) + 
#         geom_rect(aes(xmin=32.595, xmax=48.405, ymin=-0.01, ymax=-0.005),
#                   fill=group.colors[3], alpha=0.05) + 
#         geom_rect(aes(xmin=48.595, xmax=64.405, ymin=-0.01, ymax=-0.005),
#                   fill=group.colors[4], alpha=0.08) + 
#         geom_rect(aes(xmin=64.595, xmax=80.405, ymin=-0.01, ymax=-0.005),
#                   fill=group.colors[5], alpha=0.15) + 
#         geom_rect(aes(xmin=80.595, xmax=96.5, ymin=-0.01, ymax=-0.005),
#                   fill=group.colors[6], alpha=0.1)
#     message("Mutational signature plot generation done!")
#     return(pic)
# }
# 
