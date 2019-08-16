#' A VAF plot painter
#' @description Read maf file as data.frame, cluster data with maftools
#'  and finally draw variant allele frequency(VAF) frequency distribution
#'  curve with ggplot2 as well as ggridges. We could use different parameters
#'  to control output images from different samples or conclude all samples' 
#'  VAF information in one image.
#' 
#' @import ggplot2 ggsci 
#' @importFrom maftools inferHeterogeneity
#' @importFrom ggridges geom_density_ridges
#' 
#' @param maf_file specify a maf document/directory as the input of the 
#' function
#' @param sampleOption specify single/all sample names 
#' (Tumor_Sample_Barcodes, tsb). Default "ridges".
#' @param themeOption select a coloring scheme from ggsci. 
#' Default "aaas".
#' @param fileFormat choose an output file format accessable for ggsave. 
#' Default "png".
#' @return Images of selected samples' VAF
#' 
#' @export plotVAF
#'
#' @examples
#' \dontrun{
#' VAF_plot(maf_file, sampleOption="ridges", themeOption="aaas") 
#' VAF_plot(maf_file, sampleOption="allSeparate") 
#' VAF_plot(maf_file, sampleOption="allCombined") 
#' VAF_plot(maf_file, sampleOption="tsb1", fileFormat="pdf") 
#' 
#'}

## Main function for VAF plot
vafCluster <-function(maf, plotOption="ridges", 
                      themeOption="aaas", fileFormat="png", 
                      showMATH=TRUE, outputDir=NULL){
    ## check the output directory
    if (is.null(outputDir)){
        warning("NOTE: It is recommended to provide proper output directory for pictures")
    }
    ## original data preparation
    ## read .maf file
    mafInput <- maf@data
    laml <- maf
    ## specify patienID
    patientID <- maf@patientID
    ## extract vaf info
    n <- length(mafInput$Hugo_Symbol)
    vafInputMt <- data.frame(mafInput$Hugo_Symbol, 
                             mafInput$VAF, 
                             mafInput$Tumor_Sample_Barcode)
    colnames(vafInputMt) <- c("Hugo_Symbol", "VAF", "Samples")
    clusterAll <- data.frame()
    ## extract all tumor sample barcode
    tsbLs <- data.frame(unique(vafInputMt$Samples))
    colnames(tsbLs) <- c("samples")
    
    ## print all samples respectively
    if ((plotOption == "allSeparate") | (plotOption == "allCombined")){
        lsPicName <- c()
        for (counterMt in seq_along(tsbLs[,1])){
            sampleName <- as.character(tsbLs[,1][counterMt])
            ## calculate ScoreMATH
            mathscore <- .mathCal(mafInput, showMATH, 
                                  plotOption, sampleName)
            sampleMt <- vafInputMt[which(
                vafInputMt$Samples %in% sampleName),]
            clusterMt <- inferHeterogeneity(maf=laml, 
                                            tsb=as.character(
                                                sampleMt$Samples[1]), 
                                            vafCol='VAF', 
                                            useSyn=TRUE)$clusterData
            colnames(clusterMt)[6]="VAF"
            ## print VAF pictures for all samples
            if (plotOption == "allSeparate"){
                pic <- .drawVAF(clusterMt, themeOption, 
                                sampleName, mathscore)
                ggsave(pic, filename=paste(sampleName, "_VAF_Cluster", 
                                           ".", fileFormat,sep=""), 
                       width=12, height=9, dpi=800, path=outputDir)
            }
            else {
                ## print VAF pictures for all samples
                pic_cha <- paste("allSeparate", patientID, ".", counterMt, 
                                 "<-.drawVAF(clusterMt, themeOption, ", 
                                 "sampleName, mathscore, ", 
                                 "MIXOption=plotOption)", sep="")
                eval(parse(text=pic_cha))
                pic_name <- paste("allSeparate", patientID, ".", counterMt, sep="")
                lsPicName <- c(lsPicName, pic_name)
            }
            
        }
        if (plotOption == "allCombined"){
            ## set the columns of the picture and generate all single pictures
            pic <- eval(parse(text=paste("plot_grid(", 
                                         paste(lsPicName, collapse=","), 
                                         ", nrow=", 
                                         ceiling(length(lsPicName)/2), 
                                         ", ncol=2, align=\"v\")" , 
                                         sep="")))
            ## save the cowplot picure
            ggsave(pic, filename=paste(sampleName, "_VAF_Cluster_MIX", ".", 
                                       fileFormat,sep=""), 
                   width=12, height=9, dpi=800, path=outputDir)
        }
    }
    
    ## one pic for all sample
    else if (plotOption == "ridges"){
        ## collect all samples' cluster results
        for (counterMt in seq_along(tsbLs[,1])){
            sampleName <- as.character(tsbLs[,1][counterMt])
            ## calculate ScoreMATH
            mathscore <- .mathCal(mafInput, showMATH, 
                                  plotOption, sampleName)
            sampleMt <- vafInputMt[which(
                vafInputMt$Samples %in% sampleName),]
            clusterMtCha <- paste("clusterMt_", counterMt, 
                                  " <- inferHeterogeneity(maf=laml, ", 
                                  "tsb=as.character(sampleMt[1,3]), ", 
                                  "vafCol=\'VAF\', ", 
                                  "useSyn=TRUE)$clusterData", sep ="")
            eval(parse(text=clusterMtCha))
            clusterMtCha <- paste("colnames(clusterMt_", counterMt, 
                                  ")[6]=\"VAF\"",sep ="")
            eval(parse(text=clusterMtCha))
            clusterMtCha <- paste("clusterAll <- rbind(clusterAll, ", 
                                  "clusterMt_", counterMt, ")",sep ="")
            eval(parse(text=clusterMtCha))
        }
        colnames(clusterAll)[6]="VAF"
        pic <- eval(parse(text=.ofaVAF(clusterAll, themeOption, tsbLs, 
                                       plotOption, mathscore)))
        ggsave(pic, filename= paste(patientID, "_VAF_Cluster", ".", 
                                    fileFormat, sep=""), 
               width=12, height=9, dpi=1200, path=outputDir)
    } 
    
    ## specific sample
    else 
    {
        ## calculate ScoreMATH
        mathscore <- .mathCal(mafInput, showMATH, plotOption)
        ## data preparation
        sampleMt <- vafInputMt[which(
            vafInputMt$Samples %in% plotOption),]
        clusterMt <- inferHeterogeneity(
            maf=laml, tsb=as.character(sampleMt[1,3]), 
            vafCol='VAF', useSyn=TRUE)$clusterData
        colnames(clusterMt)[colnames(clusterMt)=="t_vaf"] <- "VAF"
        ## VAF plot for specifc sample
        pic <- .drawVAF(clusterMt, themeOption, 
                        plotOption, mathscore)
        ggsave(pic, filename=paste(plotOption,"_VAF_Cluster",".", 
                                   fileFormat,sep=""), 
               width=12, height=9, path=outputDir)
    }
    message("VAF Plot Generation Done!")
}


## Functions for all plotOption
## Calculate ScoreMATH
.mathCal <- function(mafInput, showMATH,  
                     plotOption, sampleName = ""){
    if (showMATH){
        if ((plotOption == "allSeparate") | (plotOption == "allCombined")){
            mathscore <- mathScore(mafInput, c(sampleName))
            mathscore <- mathscore[which(
                mathscore$Tumor_Sample_Barcode == sampleName), 
                ]$MATH_score
        }
        else if (plotOption == "ridges"){
            mathscore <- mathScore(mafInput, c(sampleName))
            mathscore <- mathscore[which(
                mathscore$Tumor_Sample_Barcode == "ITH MATH score"), 
                ]$MATH_score
        }
        else {
            mathscore <- mathScore(mafInput, c(plotOption))
            mathscore <- mathscore[which(
                mathscore$Tumor_Sample_Barcode == plotOption), 
                ]$MATH_score
        }
    } else {
        mathscore <- NA
    }
    mathscore
}


## Draw vlines for all plotOption
.vlineVAF <- function(clusterMt, pic, 
                      tsbLs, plotOption, 
                      tsb, ingredients=NULL)
{
    ## data prepare
    vafVlineCha <- ""
    clusterLs <- unique(clusterMt$cluster)
    ## density information of the curve for a tsb
    densityInfo <- data.frame(layer_data(pic))
    
    ## OFA specific: get scaling ratio
    if (!is.null(ingredients)){
        iscale <- ingredients$iscale[1]
        scale <- ingredients$scale[1]
    }
    
    ## Obtain vline Coordinate(x, xend, y, yend)
    for (cluster_name in clusterLs){
        x_end <- max(clusterMt[which(
            clusterMt$cluster == cluster_name), ]$VAF)
        x_end_alter <- densityInfo$x[which.min(
            abs(outer(densityInfo$x,x_end,FUN="-")))]
        y_end <- densityInfo$y[which(
            densityInfo$x == x_end_alter)]
        if (plotOption == "ridges"){
            ## Scale and draw lines
            density <- densityInfo$density[which(
                densityInfo$x == x_end_alter)]
            vafVlineCha <- paste(
                vafVlineCha, 
                "geom_segment(data=clusterMt_", which(tsbLs == tsb), 
                " ,aes(x=", x_end_alter, 
                ", xend=", x_end_alter, 
                ", y=", which(tsbLs == tsb), 
                ", yend=", which(tsbLs == tsb) + density*iscale*scale, "), ", 
                "size=0.5, colour=\"grey\", linetype=\"dashed\") + ",sep="")
        }else{
            vafVlineCha <- paste(
                vafVlineCha, 
                "geom_segment(aes(x=", x_end_alter,", xend=", x_end_alter, 
                ", y=0, yend=", y_end,"), ", 
                "size=0.5, colour=\"grey\", linetype=\"dashed\") + ",sep="")
        }
    }
    vafVlineCha
}

## Functions for specific plotOption: "allSeparate", "allCombined", "tsb"
## VAF painter
.drawVAF <- function(clusterMt, themeOption, 
                     plotOption, mathscore, 
                     MIXOption=""){
    ## A draft for density infomation(density_info) of ggplot
    picv <- ggplot(clusterMt, aes(x=VAF)) + 
        geom_line(size=1, colour="#00C0EB", stat="density")
    if (is.na(mathscore)){
        if (MIXOption == "allCombined"){
            ## generate cha for ggplot and specific titles for minifigures
            vafDrawCha <- paste("ggplot(clusterMt, aes(x=VAF)) + 
                       theme_bw() + 
                       theme(legend.position=\'none\',", 
                                " title=element_text(size=10), ", 
                                "text=element_text(size=10), ", 
                                "panel.grid=element_blank(), ", 
                                "panel.border=element_blank(), ", 
                                "axis.line=element_line(size=0.25)) + ", 
                                "geom_line(size=1, colour=\"#00C0EB\", ", 
                                "stat=\"density\") + ", 
                                "geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
                                .vlineVAF(clusterMt, picv, 
                                          tsb_ls, plotOption),
                                "scale_color_", themeOption, "() + ", 
                                "scale_fill_", themeOption, "()", sep="")
            eval(parse(text=vafDrawCha))
        } else {
            ## generate character/string for ggplot and paint the picture
            vafDrawCha <- paste("ggplot(clusterMt, aes(x=VAF)) + 
                       theme_bw() + 
                       theme(title=element_text(size=18), ", 
                                "text=element_text(size=18), ", 
                                "panel.grid=element_blank(), ", 
                                "panel.border=element_blank(), ", 
                                "axis.line=element_line(size=0.25)) + ", 
                                "geom_line(size=1, colour=\"#00C0EB\", ", 
                                "stat=\"density\") + ", 
                                "geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
                                .vlineVAF(clusterMt, picv, tsb_ls, plotOption),
                                "scale_color_", themeOption, "() + ", 
                                "scale_fill_", themeOption, "()", sep="")
            eval(parse(text=vafDrawCha))
        } 
    }
    else {
        if (MIXOption == "allCombined"){
            ## generate cha for ggplot and specific titles for minifigures
            vafDrawCha <- paste("ggplot(clusterMt, aes(x=VAF)) + 
                       theme_bw() + 
                       theme(legend.position=\'none\', ", 
                                "plot.title=element_text(size=10, hjust=1, ", 
                                "vjust=0.5, face='bold'), ", 
                                "title=element_text(size=10), ", 
                                "text=element_text(size=10), ", 
                                "panel.grid=element_blank(), ", 
                                "panel.border=element_blank(), ", 
                                "axis.line=element_line(size=0.25)) + ", 
                                "ggtitle(\"", plotOption, 
                                "\'s MATH Score: ", as.character(mathscore), 
                                "\") + ", 
                                "geom_line(size=1, colour=\"#00C0EB\", ", 
                                "stat=\"density\") + ", 
                                "geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
                                .vlineVAF(clusterMt, picv, tsb_ls, plotOption),
                                "scale_color_", themeOption, "() + ", 
                                "scale_fill_", themeOption, "()", sep="")
            eval(parse(text=vafDrawCha))
        } else {
            ## generate character/string for ggplot and paint the picture
            vafDrawCha <- paste("ggplot(clusterMt, aes(x=VAF)) + 
                       theme_bw() + 
                       theme(plot.title=element_text(size=18, hjust=1, ", 
                                "vjust=0.5, face='bold'), ", 
                                "title=element_text(size=18), ",
                                "text=element_text(size=18), ", 
                                "panel.grid=element_blank(), ", 
                                "panel.border=element_blank(), ", 
                                "axis.line=element_line(size=0.25)) + ",  
                                "ggtitle(\"MATH Score: ", 
                                as.character(mathscore), "\") + ", 
                                "geom_line(size=1, colour=\"#00C0EB\", ", 
                                "stat=\"density\") + ", 
                                "geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
                                .vlineVAF(clusterMt, picv, tsb_ls, plotOption),
                                "scale_color_", themeOption, "() + ", 
                                "scale_fill_", themeOption, "()", sep="")
            eval(parse(text=vafDrawCha))
        }
    }
}

## Functions for specific plotOption: "ridges"
## VAF painter for OFA
.ofaVAF <- function(clusterAll, themeOption, 
                    tsbLs, plotOption, 
                    mathscore){
    if (is.na(mathscore)){
        vafOFACha <- paste("ggplot(clusterAll, ", 
                           "aes(x=VAF, y=Tumor_Sample_Barcode)) +
                      theme_bw() + 
                      theme(title=element_text(size=18), ", 
                           "text=element_text(size=18), ", 
                           "panel.grid=element_blank(), ", 
                           "panel.border=element_blank(), ", 
                           "axis.line=element_line(size=0.25)) + ",
                           "geom_density_ridges(fill=\"whitesmoke\", ",
                           "calc_ecdf=TRUE, alpha=0.8) + ",
                           "geom_point(aes(x=VAF, ", 
                           "y=Tumor_Sample_Barcode, ", 
                           "color=cluster), ", 
                           "alpha=0.5, show.legend=FALSE) + ", 
                           "geom_density_ridges(color=\"#00C0EB\", ", 
                           "fill=NA, calc_ecdf=TRUE, alpha=0.5) + ", 
                           .ofaVlineVAF(clusterAll, tsbLs, plotOption),  
                           "scale_color_", themeOption, "() + ", 
                           "scale_fill_", themeOption, "()", sep="")
    } else {
        vafOFACha <- paste("ggplot(clusterAll, ", 
                           "aes(x=VAF, y=Tumor_Sample_Barcode)) +
                      theme_bw() + 
                      theme(plot.title=element_text(size=18, hjust=1, ", 
                           "vjust=0.5, face='bold'), ", 
                           "title=element_text(size=18), ", 
                           "text=element_text(size=18), ", 
                           "panel.grid=element_blank(), ", 
                           "panel.border=element_blank(), ", 
                           "axis.line=element_line(size=0.25)) + ", 
                           "ggtitle(\"MATH Score: ", 
                           as.character(mathscore), "\") + ", 
                           "geom_density_ridges(fill=\"whitesmoke\", ", 
                           "calc_ecdf=TRUE, alpha=0.8) + ",
                           "geom_point(aes(x=VAF, ", 
                           "y=Tumor_Sample_Barcode, ", 
                           "color=cluster), ", 
                           "alpha=0.5, show.legend=FALSE) + ", 
                           "geom_density_ridges(color=\"#00C0EB\", ", 
                           "fill=NA, ", 
                           "calc_ecdf=TRUE, ", 
                           "alpha=0.5) + ",
                           .ofaVlineVAF(clusterAll, tsbLs, plotOption), 
                           "scale_color_", themeOption, "() + ", 
                           "scale_fill_", themeOption, "()", sep="")
    }
    vafOFACha
}

## VAF draw vlines for ofa
.ofaVlineVAF <- function(clusterAll, tsbLs, 
                         plotOption)
{
    ## data prepare
    .ofaVlineVAF <- ""
    ## density information of the curve for all tsbs
    gr <- ggplot(clusterAll, aes(x=VAF, y=Tumor_Sample_Barcode)) + 
        geom_density_ridges() 
    ingredients <- ggplot_build(gr) %>% purrr::pluck("data", 1)
    
    for (tsb in tsbLs$samples)
    {
        ## renew tsb's data
        VAFVlineCha <- ""
        xEndLs <- data.frame()
        yEndLs <- data.frame()
        clusterMt <- clusterAll[which(
            clusterAll$Tumor_Sample_Barcode == tsb), ]
        colnames(clusterMt)[6]="VAF"
        clusterLs <- unique(clusterMt$cluster)
        ## A draft for density infomation(density_info) of ggplot
        picv <- ggplot(clusterMt, aes(x=VAF)) + 
            geom_line(size=1, colour="#00C0EB", stat="density")
        densityInfo <- data.frame(layer_data(picv))
        ## collect vlines for a single tsb
        VAFVlineCha <- .vlineVAF(clusterMt, picv, 
                                 tsbLs, plotOption, 
                                 tsb, ingredients)
        ## collect vlines for all tsbs
        .ofaVlineVAF <- paste(.ofaVlineVAF, VAFVlineCha)
    }
    .ofaVlineVAF
}
