#' A VAF plot painter
#' @description Read maf file as data.frame, cluster data with maftools
#'  and finally draw variant allele frequency(VAF) frequency distribution
#'  curve with ggplot2 as well as ggridges. We could use different parameters
#'  to control output images from different samples or conclude all samples' 
#'  VAF information in one image.
#' 
#' @import ggplot2 ggridges ggsci 
#' @importFrom maftools inferHeterogeneity
#' 
#' @param maf_file specify a maf document/directory as the input of the 
#' function
#' @param sampleOption specify single/all sample names 
#' (Tumor_Sample_Barcodes, tsb). Default "OFA".
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
#' VAF_plot(maf_file, sampleOption="OFA", themeOption="aaas") 
#' VAF_plot(maf_file, sampleOption="All") 
#' VAF_plot(maf_file, sampleOption="MIX") 
#' VAF_plot(maf_file, sampleOption="tsb1", fileFormat="pdf") 
#' 
#'}

## Main function for VAF plot
plotVAF <-function(maf, sampleOption="OFA", 
                    themeOption="aaas", fileFormat="png", 
                    showMATH=TRUE)
{
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
    if ((sampleOption == "All") | (sampleOption == "MIX")){
        lsPicName <- c()
        for (counterMt in seq_along(tsbLs[,1])){
            sampleName <- as.character(tsbLs[,1][counterMt])
            ## calculate ScoreMATH
            mathScore <- .mathCal(mafInput, showMATH, 
                                   sampleOption, sampleName)
            sampleMt <- vafInputMt[which(
                vafInputMt$Samples %in% sampleName),]
            clusterMt <- inferHeterogeneity(maf=laml, 
                                             tsb=as.character(
                                                 sampleMt$Samples[1]), 
                                             vafCol='VAF', 
                                             useSyn=TRUE)$clusterData
            colnames(clusterMt)[6]="VAF"
            ## print VAF pictures for all samples
            if (sampleOption == "All"){
                pic <- .drawVAF(clusterMt, themeOption, 
                                sampleName, mathScore)
                ggsave(pic, filename=paste(sampleName, "_VAF_Cluster", 
                                           ".", fileFormat,sep=""), 
                       width=12, height=9, dpi=800, path="./output")
            }
            else {
                ## print VAF pictures for all samples
                pic_cha <- paste("all", patientID, ".", counterMt, 
                                 "<-.drawVAF(clusterMt, themeOption, ", 
                                 "sampleName, mathScore, ", 
                                 "MIXOption=sampleOption)", sep="")
                eval(parse(text=pic_cha))
                pic_name <- paste("all", patientID, ".", counterMt, sep="")
                lsPicName <- c(lsPicName, pic_name)
            }

        }
        if (sampleOption == "MIX"){
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
                   width=12, height=9, dpi=800, path="./output")
        }
    }
    
    ## one pic for all sample
    else if (sampleOption == "OFA"){
        ## collect all samples' cluster results
        for (counterMt in seq_along(tsbLs[,1])){
            sampleName <- as.character(tsbLs[,1][counterMt])
            ## calculate ScoreMATH
            mathScore <- .mathCal(mafInput, showMATH, 
                                   sampleOption, sampleName)
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
                                       sampleOption, mathScore)))
        ggsave(pic, filename= paste(patientID, "_VAF_Cluster", ".", 
                                    fileFormat, sep=""), 
               width=12, height=9, dpi=1200)
    } 
    
    ## specific sample
    else 
    {
        ## calculate ScoreMATH
        mathScore <- .mathCal(mafInput, showMATH, sampleOption)
        ## data preparation
        sampleMt <- vafInputMt[which(
            vafInputMt$Samples %in% sampleOption),]
        clusterMt <- inferHeterogeneity(
            maf=laml, tsb=as.character(sampleMt[1,3]), 
            vafCol='VAF', useSyn=TRUE)$clusterData
        colnames(clusterMt)[colnames(clusterMt)=="t_vaf"] <- "VAF"
        ## VAF plot for specifc sample
        pic <- .drawVAF(clusterMt, themeOption, 
                        sampleOption, mathScore)
        ggsave(pic, filename=paste(sampleOption,"_VAF_Cluster",".", 
                                   fileFormat,sep=""), 
               width=12, height=9)
    }
    message("VAF Plot Generation Done!")
}


## Functions for all sampleOption
## Calculate ScoreMATH
.mathCal <- function(mafInput, showMATH,  
                     sampleOption, sampleName = ""){
    if (showMATH){
        if ((sampleOption == "All") | (sampleOption == "MIX")){
            mathScore <- ScoreMATH(mafInput, c(sampleName))
            mathScore <- mathScore[which(
                mathScore$Tumor_Sample_Barcode == sampleName), 
                ]$MATH_score
        }
        else if (sampleOption == "OFA"){
            mathScore <- ScoreMATH(mafInput, c(sampleName))
            mathScore <- mathScore[which(
                mathScore$Tumor_Sample_Barcode == "ITH MATH score"), 
                ]$MATH_score
        }
        else {
            mathScore <- ScoreMATH(mafInput, c(sampleOption))
            mathScore <- mathScore[which(
                mathScore$Tumor_Sample_Barcode == sampleOption), 
                ]$MATH_score
        }
    } else {
        mathScore <- NA
    }
    mathScore
}


## Draw vlines for all sampleOption
.vlineVAF <- function(clusterMt, pic, 
                      tsbLs, sampleOption, 
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
        if (sampleOption == "OFA"){
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

## Functions for specific sampleOption: "All", "MIX", "tsb"
## VAF painter
.drawVAF <- function(clusterMt, themeOption, 
                     sampleOption, mathScore, 
                     MIXOption=""){
    ## A draft for density infomation(density_info) of ggplot
    picv <- ggplot(clusterMt, aes(x=VAF)) + 
        geom_line(size=1, colour="#00C0EB", stat="density")
    if (is.na(mathScore)){
        if (MIXOption == "MIX"){
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
                                 tsb_ls, sampleOption),
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
                       .vlineVAF(clusterMt, picv, tsb_ls, sampleOption),
                       "scale_color_", themeOption, "() + ", 
                       "scale_fill_", themeOption, "()", sep="")
            eval(parse(text=vafDrawCha))
        } 
    }
    else {
        if (MIXOption == "MIX"){
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
                       "ggtitle(\"", sampleOption, 
                            "\'s MATH Score: ", as.character(mathScore), 
                            "\") + ", 
                       "geom_line(size=1, colour=\"#00C0EB\", ", 
                                 "stat=\"density\") + ", 
                       "geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
                       .vlineVAF(clusterMt, picv, tsb_ls, sampleOption),
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
                                as.character(mathScore), "\") + ", 
                       "geom_line(size=1, colour=\"#00C0EB\", ", 
                                "stat=\"density\") + ", 
                       "geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
                       .vlineVAF(clusterMt, picv, tsb_ls, sampleOption),
                       "scale_color_", themeOption, "() + ", 
                       "scale_fill_", themeOption, "()", sep="")
            eval(parse(text=vafDrawCha))
        }
    }
}

## Functions for specific sampleOption: "OFA"
## VAF painter for OFA
.ofaVAF <- function(clusterAll, themeOption, 
                    tsbLs, sampleOption, 
                    mathScore){
    if (is.na(mathScore)){
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
                            .ofaVlineVAF(clusterAll, tsbLs, sampleOption),  
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
                                as.character(mathScore), "\") + ", 
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
                            .ofaVlineVAF(clusterAll, tsbLs, sampleOption), 
                            "scale_color_", themeOption, "() + ", 
                            "scale_fill_", themeOption, "()", sep="")
    }
    vafOFACha
}

## VAF draw vlines for ofa
.ofaVlineVAF <- function(clusterAll, tsbLs, 
                         sampleOption)
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
                                   tsbLs, sampleOption, 
                                   tsb, ingredients)
        ## collect vlines for all tsbs
        .ofaVlineVAF <- paste(.ofaVlineVAF, VAFVlineCha)
    }
    .ofaVlineVAF
}
