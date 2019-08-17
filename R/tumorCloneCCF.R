#' @description Plot Tumor clone and corresponding ccf distribution
#' #importFrom ggridges scale_discrete_manual
#
#' @param maf return from readMaf()
#' @return radar plot of clone ccf and dot plot of mutation ccf
#'
#' @examples
#' pyCloneCluster <- system.file("extdata/multi_lesion/ccf", "311252.cluster.tsv", package = "Meskit")
#' pyCloneLoci <- system.file("extdata/multi_lesion/ccf", "311252.loci.tsv", package = "Meskit")
#' tumorClonesPlot(maf)
#' @export tumorClonesPlot

#tumorClonesPlot <- function(ccfClusterFile = NULL, ccfLociFile = NULL, savePlot  = FALSE, out.dir = NULL, clone.min.mut = 5, clone.min.aveCCF = 0.1){
tumorClonesPlot <- function(maf, savePlot  = FALSE, out.dir = NULL, clone.min.mut = 5, clone.min.aveCCF = 0.1){
  ccfCluster <- maf@ccf.cluster
  if(is.null(ccfCluster)){
    stop("ccf data of cluster was not found when readMaf")
  }
  ccfCluster <- ccfCluster[which(ccfCluster$size>=5 & ccfCluster$mean>=0.1),c(1,2,4)]
  cluster.ids <- unique(ccfCluster$cluster_id)


  ccfLoci <- maf@ccf.loci
  if(is.null(ccfCluster)){
    stop("ccf data of loci was not found when readMaf")
  }
  ccfLoci <- ccfLoci[which(!is.na(ccfLoci$variant_allele_frequency)),]
  loci <- ccfLoci[which(ccfLoci$cluster_id %in% cluster.ids),]
  

  ##only for our test sample
  #ccf$sample_id <- unlist(lapply(ccf$sample_id, function(x) unlist(strsplit(x, "-"))[2]))
  #loci$sample_id <- unlist(lapply(loci$sample_id, function(x) unlist(strsplit(x, "-"))[2]))

  radar_plot <- ggplot(ccfCluster, aes(x=sample_id, y=mean, fill=factor(cluster_id))) + geom_bar(stat="identity") + coord_polar()+ 
          theme_bw() + scale_fill_npg(labels=as.character(seq(length(cluster.ids))))+
            theme(panel.grid = element_blank(),
            panel.border= element_blank(),
            axis.text.x = element_text(size = 9, hjust = 1, face = "bold"), 
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())+
          labs(fill = "Tumor clone")+
          theme(legend.position = "top")+
          theme(legend.title = element_text(size = 11, face = "bold"))+   
          theme(legend.text = element_text(size = 10, face = "bold"))+ 
          theme(plot.margin = unit(c(0.2, 0.05, 0.2, 0.05), "inches"))

  loci_dotplot <- ggplot(data=loci, aes(x=sample_id, y=cellular_prevalence, color=factor(cluster_id))) + 
          geom_point(size=1.5)+ theme_bw()+ scale_color_npg()+
          theme(panel.border = element_blank(),axis.line = element_line(colour = "black"))+
          theme(panel.grid = element_blank(),
          panel.border= element_blank(),
          axis.text.x = element_text(size = 10, angle = 45, vjust = 1,hjust = 1, face = "bold"), 
          axis.text.y = element_text(size = 10, hjust = 1, face = "bold"),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 11, face = "bold"))+
          labs(y = "Cellular prevalence")+
          theme(legend.position = "none")

  radar_dotplot <- grid.arrange(radar_plot, loci_dotplot, ncol = 2)
  radar_dotplot
  if (savePlot){
    ggsave(paste(out.dir, "radar_dotplot.pdf", sep = "/"), radar_dotplot, width = 9, height = 6.5) 
  }
   
}




