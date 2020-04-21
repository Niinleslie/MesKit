Pathway_analysis <- function(genes = NULL, pathway.type = pathway.type, pval = pval, pAdjustMethod = "BH", qval = qval, 
                             patientID = patientID, Name = Name){
   
   trans = suppressMessages(bitr(genes, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db"))
   
   pathway.type <- toupper(pathway.type)
   pathway.type <- match.arg(pathway.type, c("KEGG", "REACTOME"))
   
   if (pathway.type == "KEGG"){
      pathway <- enrichKEGG(
         gene          = trans$ENTREZID,                      
         organism      = 'hsa',
         keyType       = 'kegg',                    
         pvalueCutoff  = pval,
         pAdjustMethod = pAdjustMethod,
         qvalueCutoff  = qval,
      )
   }
   else{
      pathway <- enrichPathway(
         gene          = trans$ENTREZID,                              
         organism      = 'human',                              
         pvalueCutoff  = pval,
         pAdjustMethod = pAdjustMethod,
         qvalueCutoff  = qval,
      )
      
   }
   
   
   if (!is.null(pathway) && nrow(pathway@result) > 0){     
      if(Name == "All"){
         pathway@result$Case <- patientID
      }
      else{
         pathway@result$branch <- Name
      }
   } 
   return(pathway)
}  

doTreePathway <- function(phyloTree = NULL,
                          selectedGenes = NULL,
                          pathway.type="KEGG",
                          pval=0.05,
                          pAdjustMethod="BH",
                          qval=0.2,
                          plotType="dot",
                          showCategory =  5,
                          withinType = FALSE){
   branches <- phyloTree@mut.branches
   patientID <- phyloTree@patientID
   
   Pathway.branch.result <- data.frame()
   all.genes <- c()
   pathPlot.list <- list()
   pathResult.list <- list()
   plot.branchNames <- c()
   result.branchNames <- c()
   
   mut.ref <- plyr::rbind.fill(phyloTree@mut.branches)
   
   if(withinType){
      branches <- unique(mut.ref$Branch_Tumor_Type)
   }
   else{
      branches <- unique(mut.ref$Branch_ID)
   }
   
   x <- 1
   y <- 1 
   for (i in 1:length(branches)){
      branchID <- branches[[i]]
      
      if(withinType){
         branch.ref <- mut.ref[mut.ref$Branch_Tumor_Type %in% branchID,]
         
      }
      else{
         branch.ref <- mut.ref[mut.ref$Branch_ID %in% branchID,]
      }
      
      ids <- unique(branch.ref$mut_id)
      if(length(ids) == 1){
          if(ids == "NoSigTag"){
              next
          }
      }
      
      #split the gene symbol by ","
      geneSymbol <- unique(unlist(strsplit(as.character(branch.ref$Hugo_Symbol), split = ",")))
      if(!is.null(selectedGenes)){
         geneSymbol <- geneSymbol[geneSymbol %in% selectedGenes]
      }
      all.genes <- unique(c(all.genes, geneSymbol))
      
      message(paste("Processing branch: ", branchID, sep = ""))
      Pathway.branch <- Pathway_analysis(geneSymbol, pathway.type, pval, pAdjustMethod,
                                         qval, patientID, branchID)
      if(!is.null(Pathway.branch)){
         pathResult.list[[x]] <- Pathway.branch@result
         result.branchNames <- c(result.branchNames, branchID)
         x <- x+1
         if(min(Pathway.branch@result$p.adjust) > pval | min(Pathway.branch@result$qvalue, na.rm = T) > qval){
            message(paste("0 enriched pathway found for branch ", branchID, sep = ""))
         }else{
            plot.branchNames <- c(plot.branchNames, branchID)
            if (plotType == "dot"){
               path.plot <- dotplot(Pathway.branch, showCategory = showCategory) + ggtitle(branchID)
            }else if (plotType == "bar"){
               path.plot <- barplot(Pathway.branch, showCategory = showCategory) + ggtitle(branchID)
            }
            pathPlot.list[[y]] <- path.plot
            y <- y+1
         }      
      }
      
   }  
   Pathway.all <- Pathway_analysis(all.genes, pathway.type, pval, pAdjustMethod,
                                   qval, patientID, Name = "All")
   Pathway.all.result <- Pathway.all@result
   pathResult.list[[x]] <- Pathway.all.result
   result.branchNames <- c(result.branchNames, "All")
   if(min(Pathway.all.result$p.adjust) > pval | min(Pathway.all.result$qvalue, na.rm = T) > qval){
      message(paste("0 enriched pathway found in ", patientID, sep = ""))
   }else{
      plot.branchNames <- c(plot.branchNames, "All")
      if (plotType == "dot"){
         path.plot <- dotplot(Pathway.all, showCategory = showCategory) + ggtitle(Pathway.all$Case)
      }else if (plotType == "bar"){
         path.plot <- barplot(Pathway.all, showCategory = showCategory) + ggtitle(Pathway.all$Case)
      }
      pathPlot.list[[y]] <- path.plot
   }
   
   names(pathResult.list) <- result.branchNames
   names(pathPlot.list) <- plot.branchNames
   
   
   pathway <- list(pathResult.list, pathPlot.list)
   names(pathway) <- c("pathway.category", "pathway.plot")
   
   return(pathway)
   
}