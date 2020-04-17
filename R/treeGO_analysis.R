
GO_analysis <- function(genes=NULL, GO.type=GO.type, pval=pval, pAdjustMethod=pAdjustMethod,
                        qval=qval, patientID=patientID, name=name){
   GO.type <- toupper(GO.type)
   GO.type <- match.arg(GO.type, c("BP", "MF", "CC", "ALL"))
   
   ego <- enrichGO(gene          = genes,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = GO.type,
                   pvalueCutoff  = pval,
                   pAdjustMethod = pAdjustMethod,
                   qvalueCutoff  = qval
   )
   
   if (!is.null(ego) && nrow(ego@result) > 0){
      if(name == "All"){
         ego@result$Case <- patientID
      }else{
         ego@result$branch <- name
      }
   }
   return(ego)
}

doTreeGO <- function(phyloTree = NULL,
                     selectedGenes = NULL,
                     GO.type="BP",
                     pval=0.05,
                     pAdjustMethod="BH", 
                     qval=0.2,
                     plotType="dot",
                     showCategory=5,
                     withinType = withinType){
   branches <- phyloTree@mut.branches
   patientID <- phyloTree@patientID
   
   GO.branch.result <- data.frame()
   all.genes <- c()
   egoPlot.list <- list()
   egoResult.list <- list() 
   plot.branchNames <- c()
   result.branchNames <- c()
   
   x <- 1
   y <- 1 
   
   mut.ref <- rbind.fill(phyloTree@mut.branches)
   
   if(withinType){
      branches <- unique(mut.ref$Branch_Tumor_Type)
   }
   else{
      branches <- unique(mut.ref$Branch_ID)
   }
   
   
   for (i in 1:length(branches)){
      branchID <- branches[[i]]
      
      if(withinType){
         branch.ref <- mut.ref[mut.ref$Branch_Tumor_Type %in% branchID,]
      }
      else{
         branch.ref <- mut.ref[mut.ref$Branch_ID %in% branchID,]
      }
      
      #split the gene symbol by ","
      geneSymbol <- unique(unlist(strsplit(as.character(branch.ref$Hugo_Symbol), split = ",")))
      
      if(!is.null(selectedGenes)){
         geneSymbol <- geneSymbol[geneSymbol %in% selectedGenes]
      }
      all.genes <- unique(c(all.genes, geneSymbol))
      
      message(paste("Processing branch: ", branchID, sep = ""))
      ego.branch <- GO_analysis(geneSymbol, GO.type, pval, pAdjustMethod,
                                qval, patientID, branchID)
      
      if(!is.null(ego.branch)){
         egoResult.list[[x]] <- ego.branch@result
         result.branchNames <- c(result.branchNames, branchID)
         x <- x+1
         if(min(ego.branch@result$p.adjust) > pval | min(ego.branch@result$qvalue, na.rm=T) > qval){
            message(paste("0 enriched term found for branch ", branchID, sep = ""))
         }
         else{      
            plot.branchNames <- c(plot.branchNames, branchID)
            if (plotType == "dot"){
               go.plot <- dotplot(ego.branch, showCategory=showCategory) + ggtitle(branchID)
            }else if (plotType == "bar"){
               go.plot <- barplot(ego.branch, showCategory=showCategory) + ggtitle(branchID)
            }
            egoPlot.list[[y]] <- go.plot
            y <- y+1
         }
      }
   }
   
   ego.all <- GO_analysis(unique(all.genes), GO.type, pval, pAdjustMethod,
                          qval, patientID, name="All")
   ego.all.result <- ego.all@result
   egoResult.list[[x]] <- ego.all.result
   result.branchNames <- c(result.branchNames, "All")
   if(min(ego.all.result$p.adjust) > pval | min(ego.all.result$qvalue, na.rm = T) > qval){
      message(paste("0 enriched terms found in ", patientID, sep = ""))
   }else{
      plot.branchNames <- c(plot.branchNames, "All")
      if (plotType == "dot"){
         go.plot <- dotplot(ego.all, showCategory=showCategory) + ggtitle(ego.all.result$Case)
      }else if (plotType == "bar"){
         go.plot <- barplot(ego.all, showCategory=showCategory)+ ggtitle(ego.all.result$Case)
      }
      egoPlot.list[[y]] <- go.plot  
   }
   
   names(egoResult.list) <- result.branchNames
   names(egoPlot.list) <- plot.branchNames
   ego <- list(egoResult.list, egoPlot.list)
   names(ego) <- c("GO.category", "GO.plot")
   
   return(ego)
}