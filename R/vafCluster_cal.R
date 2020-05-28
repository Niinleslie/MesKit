
# .clusterGenerator <- function(data){
#    ## #Use Gaussian finite mixture model
#    cluster_result <- mclust::densityMclust(maf.dat$VAF, G=1:7, verbose=FALSE)
#    data$cluster <- as.character(cluster_result$classification)
#    return(subdata)
# }


## Functions for specific plotOption: "separate", "combine", "tsb"
## VAF painter
.drawVAF <- function(subdata, plotOption, id){
   ## A draft for density infomation(density_info) of ggplot
   picv <- ggplot(subdata, aes(x=VAF)) + 
      geom_line(size=1, colour="#00C0EB", stat="density")
   
   # build color vector for later use
   color_scale <- c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", 
                    "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", 
                    "#808180FF", "#1B1919FF")
   
  if (plotOption == "combine"){
     ## generate  plot and specific titles for minifigures
     p <- ggplot(subdata, aes(x = VAF)) + 
         theme_bw() +
         theme(legend.position='right', 
               plot.title=element_text(size=13,hjust = 1,vjust = 0.5),
               panel.grid=element_blank(), 
               panel.border=element_blank(), 
               axis.line=element_line(size=0.7), 
               axis.title=element_text(size=13),
               axis.text=element_text(size=12, colour = "black"))+
          geom_line(size=1, colour="#00C0EB", stat="density") + 
          geom_rug(aes(y=0, colour=cluster), sides="b") + 
          .vlineVAF(subdata, picv, tsb_ls, plotOption) + 
          geom_rug(aes(y=0, colour=cluster), sides="b") + 
          scale_colour_manual(values=color_scale) + 
          ggtitle(id)
          labs(y = "Density")
     
  } else {
     ## paint the picture
     p <-  ggplot(subdata, aes(x=VAF))+ 
                  theme_bw() + 
                  theme(plot.title=element_text(size=13),
                        panel.grid=element_blank(), 
                        panel.border=element_blank(),
                         axis.line=element_line(size=0.7),
                         axis.title=element_text(size=13), 
                         axis.text=element_text(size=12, colour = "black")
                        )+
                 geom_line(size=1, colour="#00C0EB",stat="density") +  
                 .vlineVAF(subdata, picv, tsb_ls, plotOption) +
                 geom_rug(aes(y=0, colour=cluster), sides="b") + 
                 scale_colour_manual(values=color_scale) + 
                 labs(y = "Density")
  } 
   return(p)
}


## Draw vlines for all plotOption
.vlineVAF <- function(subdata, picv, 
                      tsbLs, plotOption, 
                      tsb, ingredients=NULL)
{
   ## data prepare
   vafVlineCha <- ""
   cluster_list <- unique(subdata$cluster)
   ## density information of the curve for a tsb
   densityInfo <- data.frame(layer_data(picv))
   
   
   ## OFA specific: get scaling ratio
   # if (!is.null(ingredients)){
   #    iscale <- ingredients$iscale[1]
   #    scale <- ingredients$scale[1]
   # }
   vline <- NULL
   df_vline <- data.frame()
   ## Obtain vline Coordinate(x, xend, y, yend)
   for (cluster_name in cluster_list){
      x_end <- max(subdata[which(
         subdata$cluster == cluster_name), ]$VAF)
      x_end_alter <- densityInfo$x[which.min(
         abs(outer(densityInfo$x,x_end,FUN="-")))]
      y <- 0
      y_end <- densityInfo$y[which(
         densityInfo$x == x_end_alter)]
      sub <- data.frame(x = x_end_alter, xend = x_end_alter, y = 0, yend = y_end)
      df_vline <- rbind(df_vline, sub)
   }
   # if (plotOption == "compare"){
   #     ## Scale and draw lines
   #     density <- densityInfo$density[which(
   #         densityInfo$x == x_end_alter)]
   #     vline <- vline + 
   #         geom_segment(data=subdata,
   #                      aes(x= x_end_alter, xend= x_end_alter,
   #                          y= which(tsbLs == tsb), yend= which(tsbLs == tsb) + density*iscale*scale,  
   #                          size=0.5, colour= "grey ", linetype= "dashed "))
   # }
   vline <- geom_segment(aes(x= x, xend= xend, y= y, yend= yend),
                         data = df_vline,
                        size=0.5, colour= "grey", linetype= "dashed")
   return(vline)
}

## Functions for specific plotOption: "compare"
## VAF painter for OFA
.ofaVAF <- function(maf_data, tsbLs, 
                    plotOption, mathscore, patient, 
                    min.vaf, max.vaf){
    pic <- ggplot(maf_data,aes(x=VAF, y=ID)) +
           theme_bw() + 
            theme(plot.title=element_text(size=13, hjust=1,vjust=0.5, face='bold'), 
                  panel.grid=element_blank(), 
                  panel.border=element_blank(),
                  axis.title=element_text(size=13), 
                  axis.text=element_text(size=12, colour = "black"),
                  axis.line=element_line(size=0.7)) +  
            ggtitle("VAF clustering of ", patient) +  
            ggridges::geom_density_ridges(fill="whitesmoke",calc_ecdf=TRUE, alpha=0.5) + 
            geom_point(aes(x=VAF, y=ID, color=cluster),alpha=0.5, show.legend=T) +  
            ggridges::geom_density_ridges(color="#00C0EB",fill=NA, calc_ecdf=TRUE, alpha=0.5, size=1)+
            .ofaVlineVAF(clusterAll, tsbLs, plotOption)+ 
            scale_colour_manual(values=color_scale) +  
            labs(y = "Sample")+  
            geom_text(data=cbind(maf_data %>% 
                                     dplyr::group_by(ID) %>% 
                                     dplyr::summarise(), 
                                 VAF=  maf_data %>% 
                                           dplyr::group_by(ID) %>% 
                                           dplyr::summarise(VAF=max(VAF))$VAF),
                       aes(x=0.85*max(VAF), 
                           label=paste("MATH Score:", sprintf("%1.3f", MATH), sep="")), 
                       position=position_nudge(y=0.5), colour="black", size=3.5) +
           scale_x_continuous(limits = c(min.vaf,max.vaf))
  # pic <- paste("ggplot(clusterAll, ", 
  #                    "aes(x=VAF, y=Tumor_Sample_Barcode)) +
  #                  theme_bw() + 
  #                  theme(plot.title=element_text(size=13, hjust=1, ", 
  #                    "vjust=0.5, face='bold'), ", 
  #                    "panel.grid=element_blank(), ", 
  #                    "panel.border=element_blank(), ", 
  #                    "axis.title=element_text(size=13), ", 
  #                    "axis.text=element_text(size=12, colour = \"black\"), ", 
  #                    "axis.line=element_line(size=0.7)) + ", 
  #                    "ggtitle(\"VAF clustering of ", patient, "\") + ", 
  #                    "ggridges::geom_density_ridges(fill=\"whitesmoke\", ", 
  #                    "calc_ecdf=TRUE, alpha=0.5) + ",
  #                    "geom_point(aes(x=VAF, ", 
  #                    "y=Tumor_Sample_Barcode, ", 
  #                    "color=cluster), ", 
  #                    "alpha=0.5, show.legend=T) + ", 
  #                    "ggridges::geom_density_ridges(color=\"#00C0EB\", ", 
  #                    "fill=NA, calc_ecdf=TRUE, alpha=0.5, size=1) + ",
  #                    .ofaVlineVAF(clusterAll, tsbLs, plotOption), 
  #                    "scale_colour_manual(values=color_scale) + ", 
  #                    "labs(y = \"Sample\") + ", 
  #                    "geom_text(data=cbind(clusterAll %>% dplyr::group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(), 
  #                  MATH=unique(clusterAll$MATH), 
  #                  VAF=(clusterAll %>% dplyr::group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(VAF=max(VAF)))$VAF),
  #                  aes(x=0.85*max(VAF), label=paste(\"MATH Score: \", sprintf(\"%1.3f\", MATH), sep=\"\")), 
  #                  position=position_nudge(y=0.5), colour=\"black\", size=3.5) + ", 
  #                    "scale_x_continuous(limits = ", "c(", as.character(min.vaf), ",", as.character(max.vaf), "))", sep="")

    return(pic)
}


## VAF draw vlines for ofa
.ofaVlineVAF <- function(clusterAll, tsbLs, 
                         plotOption)
{
   ## data prepare
   .ofaVlineVAF <- ""
   ## density information of the curve for all tsbs
   gr <- ggplot(clusterAll, aes(x=VAF, y=Tumor_Sample_Barcode)) + 
      ggridges::geom_density_ridges() 
   ingredients <- ggplot_build(gr) %>% purrr::pluck("data", 1)
   
   for (tsb in tsbLs$samples)
   {
      ## renew tsb's data
      VAFVlineCha <- ""
      xEndLs <- data.frame()
      yEndLs <- data.frame()
      subdata <- clusterAll[which(
         clusterAll$Tumor_Sample_Barcode == tsb), ]
      cluster_list <- unique(subdata$cluster)
      ## A draft for density infomation(density_info) of ggplot
      picv <- ggplot(subdata, aes(x=VAF)) + 
         geom_line(size=1, colour="#00C0EB", stat="density")
      densityInfo <- data.frame(layer_data(picv))
      ## collect vlines for a single tsb
      VAFVlineCha <- .vlineVAF(subdata, picv, 
                               tsbLs, plotOption, 
                               tsb, ingredients)
      ## collect vlines for all tsbs
      .ofaVlineVAF <- paste(.ofaVlineVAF, VAFVlineCha)
   }
   .ofaVlineVAF
}


doVafCluster <- function(maf = NULL,
                         seg = NULL,
                         chrSilent = NULL,
                         mutType = "All",
                         use.indel = TRUE,
                         min.vaf=0.02,
                         max.vaf=1,
                         showMATH=TRUE, 
                         plotOption="combine",
                         use.adjVAF = FALSE){
    ## remove mutation in CNA regions
    if(!is.null(seg)){
        maf <- copyNumberFilter(maf,seg)
    }
   maf_data <- subsetMaf(maf,
                     chrSilent = chrSilent,
                     mutType = mutType,
                     use.indel = use.indel,
                     min.vaf = min.vaf,
                     max.vaf = max.vaf,
                     use.adjVAF = use.adjVAF)
   maf@data <- maf_data
   patient <- unique(maf_data$Patient_ID) 
   
   ## extract vaf info
   n <- length(maf_data$Hugo_Symbol)
   vafInputMt <- data.frame(maf_data$Hugo_Symbol, 
                            maf_data$VAF, 
                            maf_data$Tumor_Sample_Barcode)
   colnames(vafInputMt) <- c("Hugo_Symbol", "VAF", "Samples")
   clusterAll <- data.frame()
   ## extract all tumor sample barcode
   tsbLs <- data.frame(unique(vafInputMt$Samples))
   colnames(tsbLs) <- c("samples")
   
   # if(use.shiny){
   #     incProgress(amount=1)
   #     setProgress(message = paste('Generating ', "VAF density plot - ", patient, sep=""))
   # }
   
   
   # build color vector for later use
   color_scale <- c("#3B4992FF", "#EE0000FF", "#008B45FF", "#631879FF", 
                    "#008280FF", "#BB0021FF", "#5F559BFF", "#A20056FF", 
                    "#808180FF", "#1B1919FF")
   
   
   ## plot all samples' vaf distribution
   if (plotOption == "combine"){
      ## general data process for all samples 
      lsPicName <- c()
      lsSep <- list()
      lsSampleName <- c()
      for (counterMt in seq_along(tsbLs[,1])){
         sampleName <- as.character(tsbLs[,1][counterMt])
         ## calculate ScoreMATH
         mathscore <- .mathCal(maf, min.vaf, max.vaf, showMATH, plotOption, sampleName, use.adjVAF = use.adjVAF)
         sampleMt <- vafInputMt[which(
            vafInputMt$Samples %in% sampleName),]
         ## data cleaning
         sampleMt <- sampleMt[complete.cases(sampleMt), ]
         sampleMt <- sampleMt[which(sampleMt$VAF != 0),]
         if (length(sampleMt[,1]) < 3) {
            message(paste("Sample ", sampleName, " has too few mutaions",sep = ""))
            next()
         }
         
         ## infer possible cluster from maf_data
         subdata <- .clusterGenerator(maf_data, sampleName)
         subdata <- subdata[which(subdata$Tumor_Sample_Barcode == sampleName), ]
         
         # prepare separated pictures for later combination 
         pic_cha <- paste("separate", ".", counterMt, 
                          "<-.drawVAF(subdata, ", 
                          "sampleName, mathscore, ", 
                          "MIXOption=plotOption)", sep="")
         eval(parse(text=pic_cha))
         pic_name <- paste("separate", ".", counterMt, sep="")
         lsPicName <- c(lsPicName, pic_name)
      }
      
      ## combine: print VAF pictures for all samples in one document
      if (plotOption == "combine"){
         if (showMATH){
            mathtbscoreLs <- .mathCal(maf, min.vaf, max.vaf, showMATH, "compare", sampleName, use.adjVAF = use.adjVAF)
            
            ## set the columns of the picture and generate all single pictures
            combineTitle <- cowplot::ggdraw() + 
               cowplot::draw_label(
                  paste("VAF clustering of ", patient, sep=""),
                  fontface = 'bold',
                  x = 0.5,
                  hjust = 0.5,
                  size = 16
               ) +
               theme(
                  # add margin on the left of the drawing canvas,
                  # so title is aligned with left edge of first plot
                  plot.margin = margin(0, 0, 0, 7)
               )
            pic <- eval(parse(text=paste("cowplot::plot_grid(", 
                                         paste(lsPicName, collapse=","), 
                                         ", nrow=", 
                                         ceiling(length(lsPicName)/2), 
                                         ", ncol=2, align=\"v\", scale=0.85)" , 
                                         sep="")))
            pic <- cowplot::plot_grid(
               combineTitle, pic,
               ncol = 1,
               # rel_heights values control vertical title margins
               rel_heights = c(0.1, 1))
            
         } else {
            pic <- eval(parse(text=paste("cowplot::plot_grid(", 
                                         paste(lsPicName, collapse=","), 
                                         ", nrow=", 
                                         ceiling(length(lsPicName)/2), 
                                         ", ncol=2, align=\"v\", scale=0.85)" , 
                                         sep="")))
         }
         message(paste0("VAF density plot of ", patient, " has been generated!"))
         return(suppressWarnings(suppressMessages(pic)))
         # return(suppressWarnings(suppressMessages(pic)))
      }
   }
   
   ## plot all samples' vaf distribution with ggridges
   else if (plotOption == "compare"){
      ## calculate ScoreMATH
      mathtbscoreLs <- .mathCal(maf, min.vaf, max.vaf, showMATH, plotOption, use.adjVAF = use.adjVAF)
      mathscore <- mathtbscoreLs[mathtbscoreLs$Patient_ID == patient,]
      
      ## record sample with few mutations
      fs <- c()
      ## collect all samples' cluster results
      counterMt <- 1
      count <- 1
      for (i in seq_along(tsbLs[,1])){
         sampleName <- as.character(tsbLs[,1][counterMt])
         sampleDat <- maf_data[maf_data$Tumor_Sample_Barcode == sampleName,]
         sampleMt <- vafInputMt[which(
            vafInputMt$Samples %in% sampleName),]
         ## data cleaning
         sampleMt <- sampleMt[complete.cases(sampleMt), ]
         sampleMt <- sampleMt[which(sampleMt$VAF != 0),]
         if (nrow(sampleMt) < 3) {
            message(paste("Sample ", sampleName, " has too few mutaions",sep = ""))
            clusterAll <- clusterAll[clusterAll$Tumor_Sample_Barcode != sampleName,]
            mathscore <- mathscore[mathscore$Tumor_Sample_Barcode != sampleName,]
            # maf_data <- maf_data[maf_data$Tumor_Sample_Barcode != sampleName, ]
            fs <- append(fs,sampleName)
            counterMt <- counterMt + 1
            next()
         }
         
         ## generate data from different Tumor_Sample_Barcode
         clusterMtCha1 <- paste("clusterMt_", count, 
                                " <- .clusterGenerator(sampleDat, sampleName)", sep ="")
         eval(parse(text=clusterMtCha1))
         clusterMtCha2 <- paste("clusterMt_", count, "$MATH", 
                                " <- rep(mathscore[which(mathscore$Tumor_Sample_Barcode == sampleName), ]$MATH_Score,
                            nrow(clusterMt_", count, "))", sep ="")
         eval(parse(text=clusterMtCha2))
         clusterMtCha3 <- paste("clusterAll <- rbind(clusterAll, ", 
                                "clusterMt_", count, ")",sep ="")
         eval(parse(text=clusterMtCha3))
         counterMt <- counterMt + 1
         count <- count + 1
      }
      tsbLs <- data.frame(samples = tsbLs[!tsbLs$samples %in% fs,])
      # mathscore <- mathtbscoreLs$patientLevel$MATH_Score
      pic <- suppressMessages(eval(parse(text=.ofaVAF(clusterAll, 
                                                      tsbLs, plotOption, 
                                                      mathscore, patient, 
                                                      min.vaf, max.vaf))))
      pic <- pic +theme(plot.title = element_text(hjust = 0.5))
      message(paste0("VAF density plot of ", patient, " has been generated!"))
      return(suppressWarnings(suppressMessages(pic)))
      # return(suppressWarnings(suppressMessages(pic)))
   }
   
   ## plot specific sample's vaf plot
   # else if (plotOption %in% unique(vafInputMt$Samples))
   # {
   #    ## data preparation
   #    sampleName <- plotOption
   #    sampleMt <- vafInputMt[which(vafInputMt$Samples %in% plotOption),]
   #    ## data cleaning
   #    sampleMt <- sampleMt[complete.cases(sampleMt), ]
   #    sampleMt <- sampleMt[which(sampleMt$VAF != 0),]
   #    if (length(sampleMt[,1]) < 3) {
   #       stop(paste("Sample ", sampleName, " has too few mutaions",sep = ""))
   #    }
   #    subdata <- .clusterGenerator(maf_data, sampleName)
   #    ## calculate ScoreMATH
   #    mathscore <- .mathCal(maf, min.vaf, max.vaf, showMATH, plotOption, sampleName)
   #    ## VAF plot for specifc sample
   #    pic <- .drawVAF(subdata, plotOption, mathscore)
   #    message(paste0("VAF density plot of ", patient, " has been generated!"))
   #    return(suppressWarnings(suppressMessages(pic)))
   #    # return(suppressWarnings(suppressMessages(pic)))
   # }
   
}