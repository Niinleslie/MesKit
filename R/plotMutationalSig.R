#' plotMutSig
#' @description Visualize the mutational signatures of trunk/branches in a phylogenetic tree.
#' 
#' @import ggplot2 graphics utils cowplot
#' 
#' @param tree.mutSig The output of function treeMutationalSig.
#' @return A signature plot shows the signature for each branch and presents the distribution in each mutational categories.
#' 
#' @examples
#' plotMutSig(tree.mutSig)
#' 
#' @export plotMutSig

## plot.MutationalSigs
plotMutSig <- function(tree.mutSig){
    pms.list <- lapply(tree.mutSig, doPlotMutSig)
}

doPlotMutSig <- function(tree.mutSig) {
  sigsInput <- tree.mutSig$sigsInput
  mutSigsOutput <- tree.mutSig$mutSigsOutput
  df.aetiology <- tree.mutSig$df.aetiology
  
  ## calculate the Mutation Probability
  sigsInputSum <- as.data.frame(apply(sigsInput, 1, function(x) sum(x)))
  sigsInputTrans <- as.data.frame(t(sigsInput))
  ls.branchesName <- as.character(rownames(sigsInput))
  ls.mutationType <-as.character(rownames(sigsInputTrans))
  ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  sigsInputBranches <- data.frame()
  
  ## calculation process(maybe could be replaced by lapply)
  for (branch in ls.branchesName) {
    sigsInputTrans[branch] <- sigsInputTrans[branch]/sigsInputSum[branch, ]
    signature <- as.character(mutSigsOutput[which(mutSigsOutput$branch == branch), ]$sig)
    alias <- as.character(mutSigsOutput[which(mutSigsOutput$branch == branch), ]$alias)
    sigsWeight <- mutSigsOutput[which(mutSigsOutput$branch == branch), ]$sig.prob
    aetiology <- as.character(df.aetiology[which(df.aetiology$sig == signature), ]$aeti)
    sigsInputBranch <- data.frame(sigsInputTrans[branch], 
                                  rep(branch, length(sigsInputTrans[branch])), 
                                  rep(alias, length(sigsInputTrans[branch])), 
                                  rep(signature, length(sigsInputTrans[branch])), 
                                  rep(sigsWeight, length(sigsInputTrans[branch])), 
                                  rep(aetiology , length(sigsInputTrans[branch])), 
                                  stringsAsFactors = FALSE)
    colnames(sigsInputBranch) <- c("Mutation_Probability", "Branch", "Alias", "Signature", "SigsWeight", "Aetiology")
    sigsInputBranches <- rbind(sigsInputBranches, sigsInputBranch)
  }
  
  df.sigsInputTrans <- data.frame(Mutational_Type=ls.mutationType, 
                                  Group=ls.mutationType, 
                                  sigsInputBranches, 
                                  stringsAsFactors = FALSE)
  
  ## generate Mutation Type for every column
  for (mutationGroup in ls.mutationGroup) {
    df.sigsInputTrans$Group[which(grepl(mutationGroup, df.sigsInputTrans$Group))] <- mutationGroup
  }
  
  ## specific the label order of x axis
  orderlist <- c(ls.mutationType)
  df.sigsInputTrans <- transform(df.sigsInputTrans, Mutational_Type = factor(Mutational_Type, levels = orderlist))
  df.sigsInputTrans <- df.sigsInputTrans[which(df.sigsInputTrans$Signature != "No Signature"), ]
  df.sigsInputText <- dplyr::distinct(df.sigsInputTrans, Branch, .keep_all = TRUE)
  
  CA <- grid::textGrob(expression(bold("C > A")),
    gp=grid::gpar(fontsize=6, fontface="bold"), vjust=0,hjust=1)
  CG <- grid::textGrob(expression(bold("C > G")),
    gp=grid::gpar(fontsize=6, fontface="bold"), vjust=0,hjust=1)
  CT <- grid::textGrob(expression(bold("C > T")),
    gp=grid::gpar(fontsize=6, fontface="bold"), vjust=0,hjust=1)
  TA <- grid::textGrob(expression(bold("T > A")),
    gp=grid::gpar(fontsize=6, fontface="bold"), vjust=0,hjust=1)
  TC <- grid::textGrob(expression(bold("T > C")),
    gp=grid::gpar(fontsize=6, fontface="bold"), vjust=0,hjust=1)
  TG <- grid::textGrob(expression(bold("T > G")),
    gp=grid::gpar(fontsize=6, fontface="bold"), vjust=0,hjust=1)
  
  group.colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF",
                    "#3C5488FF", "#F39B7FFF", "#8491B4FF")
  pic <- ggplot(df.sigsInputTrans, 
    aes(x=Mutational_Type, y=Mutation_Probability, group=Group, fill=Group)
    ) + 
    geom_bar(stat="identity")+ 
    theme(panel.grid=element_blank(), 
          panel.border=element_blank(), 
          panel.background = element_blank(), 
          legend.position='none', 
          # axis.text.x=element_text(size=3, angle = 45, hjust = 1, vjust = 1), 
          plot.title = element_text(size = 13,face = "bold",hjust = 0.5,vjust = 0),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=5)) +
    ## background colors
    geom_rect(aes(xmin=0, xmax=16.5, ymin=0, ymax=Inf),
              fill="#fce7e4", alpha=0.15) + 
    geom_rect(aes(xmin=16.5, xmax=32.5, ymin=0, ymax=Inf),
              fill="#ecf8fa", alpha=0.25) + 
    geom_rect(aes(xmin=32.5, xmax=48.5, ymin=0, ymax=Inf),
              fill="#dbfff9", alpha=0.05) + 
    geom_rect(aes(xmin=48.5, xmax=64.5, ymin=0, ymax=Inf),
              fill="#e4e8f3", alpha=0.08) + 
    geom_rect(aes(xmin=64.5, xmax=80.5, ymin=0, ymax=Inf),
              fill="#fdefeb", alpha=0.15) + 
    geom_rect(aes(xmin=80.5, xmax=96.5, ymin=0, ymax=Inf),
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
    scale_y_continuous(limits=c(-0.03, 0.2), breaks=seq(0, 0.2, 0.1)) + 
    ## signature notes and text parts
    geom_text(data = df.sigsInputText, 
              aes(x=-Inf, y=Inf, label=paste(Signature, ": ", 
                                             round(as.numeric(levels(SigsWeight)[SigsWeight]), 3), 
                                             "    ", 
                                             "Aetiology: ", Aetiology, sep="")), 
              hjust = -0.02, vjust = 1.5, colour="#2B2B2B", fontface = "bold", size=2.75) + 
    ## Mutational Type Labels
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