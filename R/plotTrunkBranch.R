#' @title plotTrunkBranch
#' @description Prints the distribution of branch/trunk mutations based on mutational categories.
#' 
#' @param tree.mutSig the output of function treeMutationalSig.
#' @param conf.level confidence level of the interval for wilcox.test. Default: 0.95.
#' 
#' 
#' @examples
#' plotTrunkBranch(tree.mutSig, conf.level = 0.95)
#' @return Box plots based on mutational categories
#' @import ggplot2 grDevices graphics utils cowplot
#' @export plotTrunkBranch

plotTrunkBranch <- function(tree.mutSig, conf.level=0.95) {
    TB.plot <- lapply(tree.mutSig, doPlotTrunkBranch,
                      conf.level = conf.level)
    return(TB.plot)
}


doPlotTrunkBranch <- function(tree.mutSig, conf.level){
    ## input data from tree.mutSig
    ls.BT <- .dataProcessBT(tree.mutSig)
    df.pValue <- ls.BT$df.pValue
    sigsInputBoxplot <- ls.BT$sigsInputBoxplot
    
    ## p values of mutational list
    if (df.pValue[which(df.pValue$Group == "C>A"), ]$p.value < (1-conf.level)) {
        CApV <- grid::textGrob(paste("p = ", as.character(
            round(df.pValue[which(df.pValue$Group == "C>A"), ]$p.value, 
                  digits = 3)), "*", sep=""), 
            gp=grid::gpar(fontsize=12),vjust=0,hjust=1)    
    } else {
        CApV <- grid::textGrob(expression(""), 
                               gp=grid::gpar(fontsize=12),vjust=0,hjust=1) 
    }
    
    if (df.pValue[which(df.pValue$Group == "C>G"), ]$p.value < (1-conf.level)) {
        CGpV <- grid::textGrob(paste("p = ", as.character(
            round(df.pValue[which(df.pValue$Group == "C>G"), ]$p.value, 
                  digits = 3)), "*", sep=""), 
            gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    } else {
        CGpV <- grid::textGrob(expression(""), 
                               gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    }
    
    if (df.pValue[which(df.pValue$Group == "C>T"), ]$p.value < (1-conf.level)) {
        CTpV <- grid::textGrob(paste("p = ", as.character(
            round(df.pValue[which(df.pValue$Group == "C>T"), ]$p.value, 
                  digits = 3)), "*", sep=""), 
            gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    } else {
        CTpV <- grid::textGrob(expression(""), 
                               gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    }
    
    if (df.pValue[which(df.pValue$Group == "T>A"), ]$p.value < (1-conf.level)) {
        TApV <- grid::textGrob(paste("p = ", as.character(
            round(df.pValue[which(df.pValue$Group == "T>A"), ]$p.value, 
                  digits = 3)), "*", sep=""),
            gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    } else {
        TApV <- grid::textGrob(expression(""), 
                               gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    }
    
    if (df.pValue[which(df.pValue$Group == "T>C"), ]$p.value < (1-conf.level)) {
        TCpV <- grid::textGrob(paste("p = ", as.character(
            round(df.pValue[which(df.pValue$Group == "T>C"), ]$p.value, 
                  digits = 3)), "*", sep=""),
            gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    } else {
        TCpV <- grid::textGrob(expression(""), 
                               gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    }
    
    if (df.pValue[which(df.pValue$Group == "T>G"), ]$p.value < (1-conf.level)) {
        TGpV <- grid::textGrob(paste("p = ", as.character(
            round(df.pValue[which(df.pValue$Group == "T>G"), ]$p.value, 
                  digits = 3)), "*", sep=""),
            gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    } else {
        TGpV <- grid::textGrob(expression(""), 
                               gp=grid::gpar(fontsize=12),vjust=0,hjust=1)
    }
    
    ## names of mutational list
    CA <- grid::textGrob(expression(bold("C > A")),gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
    CG <- grid::textGrob(expression(bold("C > G")),gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
    CT <- grid::textGrob(expression(bold("C > T")),gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
    TA <- grid::textGrob(expression(bold("T > A")),gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
    TC <- grid::textGrob(expression(bold("T > C")),gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
    TG <- grid::textGrob(expression(bold("T > G")),gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
    
    group.colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF",
                      "#3C5488FF", "#F39B7FFF", "#8491B4FF")
    
    pic <- ggplot(sigsInputBoxplot, aes(x=GroupBT, y=mut.frac, fill=Group)) + 
        geom_boxplot(coef=100) + 
        theme(panel.grid=element_blank(), 
              panel.border=element_blank(), 
              panel.background = element_blank(), 
              legend.position='none', 
              axis.text.x=element_text(size=10, angle = 90, vjust = 0.5, hjust=1), 
              axis.ticks.x = element_blank(), 
              axis.text.y=element_text(size=10)) + 
        ## background colors
        geom_rect(aes(xmin=0.5, xmax=2.5, ymin=0, ymax=104),
                  fill="#fce7e4", alpha=0.15) + 
        geom_rect(aes(xmin=2.5, xmax=4.5, ymin=0, ymax=104),
                  fill="#ecf8fa", alpha=0.25) + 
        geom_rect(aes(xmin=4.5, xmax=6.5, ymin=0, ymax=104),
                  fill="#dbfff9", alpha=0.05) + 
        geom_rect(aes(xmin=6.5, xmax=8.5, ymin=0, ymax=104),
                  fill="#e4e8f3", alpha=0.08) + 
        geom_rect(aes(xmin=8.5, xmax=10.5, ymin=0, ymax=104),
                  fill="#fdefeb", alpha=0.15) + 
        geom_rect(aes(xmin=10.5, xmax=12.5, ymin=0, ymax=104),
                  fill="#e5e8ef", alpha=0.1) + 
        geom_boxplot(coef=100) + 
        ## color setting
        scale_fill_manual(values=group.colors) + 
        ## axis setting
        scale_x_discrete(name = "", labels=c("Branch", "Trunk", "Branch", "Trunk", 
                                             "Branch", "Trunk", "Branch", "Trunk",
                                             "Branch", "Trunk", "Branch", "Trunk")) + 
        scale_y_continuous(name = "Mutation fraction(%)", limits=c(-5, 104), breaks=seq(0, 100, 25)) + 
        coord_cartesian(ylim = c(0,100), expand = TRUE) + 
        ## x axis bar
        geom_rect(aes(xmin=0.5, xmax=2.5, ymin=-5, ymax=-0.5),
                  fill=group.colors[1], alpha=1) + 
        geom_rect(aes(xmin=2.5, xmax=4.5, ymin=-5, ymax=-0.5),
                  fill=group.colors[2], alpha=0.25) + 
        geom_rect(aes(xmin=4.5, xmax=6.5, ymin=-5, ymax=-0.5),
                  fill=group.colors[3], alpha=0.05) + 
        geom_rect(aes(xmin=6.5, xmax=8.5, ymin=-5, ymax=-0.5),
                  fill=group.colors[4], alpha=0.08) + 
        geom_rect(aes(xmin=8.5, xmax=10.5, ymin=-5, ymax=-0.5),
                  fill=group.colors[5], alpha=0.15) + 
        geom_rect(aes(xmin=10.5, xmax=12.5, ymin=-5, ymax=-0.5),
                  fill=group.colors[6], alpha=0.1) + 
        ## Mutational Type Labels
        annotation_custom(grob = CA,  xmin = 1, xmax = 3, ymin = -8.5, ymax = -0) + 
        annotation_custom(grob = CG,  xmin = 3, xmax = 5, ymin = -8.5, ymax = -0) + 
        annotation_custom(grob = CT,  xmin = 5, xmax = 7, ymin = -8.5, ymax = -0) + 
        annotation_custom(grob = TA,  xmin = 7, xmax = 9, ymin = -8.5, ymax = -0) + 
        annotation_custom(grob = TC,  xmin = 9, xmax = 11, ymin = -8.5, ymax = -0) + 
        annotation_custom(grob = TG,  xmin = 11, xmax = 13, ymin = -8.5, ymax = -0) + 
        ## Mutational Type p value of wilcox.test
        annotation_custom(grob = CApV,  xmin = 1.5, xmax = 3, ymin = 2 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>A"), ]$mut.frac), ymax = 4 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>A"), ]$mut.frac)) +
        annotation_custom(grob = CGpV,  xmin = 3.5, xmax = 5, ymin = 2 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>G"), ]$mut.frac), ymax = 4 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>G"), ]$mut.frac)) +
        annotation_custom(grob = CTpV,  xmin = 5.5, xmax = 7, ymin = 2 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>T"), ]$mut.frac), ymax = 4 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>T"), ]$mut.frac)) +
        annotation_custom(grob = TApV,  xmin = 7.5, xmax = 9, ymin = 2 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>A"), ]$mut.frac), ymax = 4 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>A"), ]$mut.frac)) +
        annotation_custom(grob = TCpV,  xmin = 9.5, xmax = 11, ymin = 2 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>C"), ]$mut.frac), ymax = 4 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>C"), ]$mut.frac)) +
        annotation_custom(grob = TGpV,  xmin = 11.5, xmax = 13, ymin = 2 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>G"), ]$mut.frac), ymax = 4 + max(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>G"), ]$mut.frac))
    message("Branch-trunk plot generation done!")
    return(pic)
}

