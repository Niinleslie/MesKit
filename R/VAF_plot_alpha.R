#' A VAF plot painter
#' @description Read maf file as data.frame, cluster data with maftools and 
#' finally draw variant allele frequency(VAF) frequency distribution curve 
#' with ggplot2 as well as ggridges. 
#' We could use different parameters to control output images from different 
#' samples or conclude all samples' VAF information in one image.
#' 
#' @import ggplot2 maftools ggridges ggsci dplyr
#' 
#' @param maf_file specify a maf document/directory as the input of the 
#' function.
#' @param sampleOption specify single/all sample names. 
#' Default "OFA".
#' @param themeOption select a coloring scheme from ggsci. Default "aaas".
#' @param fileFormat choose an output file format accessable for ggsave. 
#' Default "png".
#' @return Images of selected samples' VAF
#'
#' @examples
#' \dontrun{
#' vafPlot(maf_file, sampleOption="OFA", themeOption="aaas") # draw a 
#' VAF image that contains all samples' VAF distribution curves with 
#' different themes.
#' vafPlot(maf_file, sampleOption="All") # draw VAF images for every 
#' sample respectively.
#' vafPlot(maf_file, sampleOption="MIX") # draw VAF images for every 
#' sample in one picture.
#' vafPlot(maf_file, sampleOption="tsb1", fileFormat="pdf") # draw a 
#' VAF image for sample tsb1 and save as a pdf file.
#' 
#'}


## import pkgs
library(ggplot2)
library(maftools)
library(ggridges)
library(ggsci)
library(dplyr)
library(cowplot)

## source scoreMATH function
setwd("/home/ninomoriaty/R_Project/MesKit/R")
source("MATH_Score.R")

## vafPlot main function
vafPlot <-function(maf, sampleOption="OFA", themeOption="aaas", 
                    fileFormat="png", showMATH=T){
    ## original data preparation
    ## read .maf file
    mafInput <- maf@data[,-c("patient", "lesion", "time")]
    laml <- read.maf(mafInput)
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
    if (sampleOption == "All"){
        for (counterMt in 1:length(tsbLs[,1])){
            sampleName <- as.character(tsbLs[,1][counterMt])
            ## calculate scoreMATH
            if (showMATH){
                MATHScore <- scoreMATH(mafInput, c(sampleName))
                MATHScore <- MATHScore[which(
                    MATHScore$Tumor_Sample_Barcode == sampleName), 
                    ]$MATH_score
            } else {
                MATHScore <- NA
            }
            sampleMt <- vafInputMt[which(
                vafInputMt$Samples %in% sampleName),]
            clusterMt <- inferHeterogeneity(
                maf=laml, tsb=as.character(
                    sampleMt[1,3]), vafCol='VAF', useSyn=TRUE)$"clusterData"
            colnames(clusterMt)[6]="VAF"
            ## print VAF pictures for all samples
            pic <- .vafDraw(clusterMt, themeOption, sampleName, MATHScore)
            ggsave(pic, filename=paste(
                sampleName, "_VAF_Cluster", ".", fileFormat,sep=""), 
                width=12, height=9, dpi=800, path="./output")
        }
    } 
    
    ## draw all figures in one file
    else if (sampleOption == "MIX"){
        lsPicName <- c()
        ## draw each pictures and name them rescpectively
        for (counterMt in 1:length(tsbLs[,1])){
            sampleName <- as.character(tsbLs[,1][counterMt])
            ## calculate scoreMATH
            if (showMATH){
                MATHScore <- scoreMATH(mafInput, c(sampleName))
                MATHScore <- MATHScore[which(
                    MATHScore$Tumor_Sample_Barcode == sampleName), 
                    ]$MATH_score
            } else {
                MATHScore <- NA
            }
            sampleMt <- vafInputMt[which(
                vafInputMt$Samples %in% sampleName),]
            clusterMt <- inferHeterogeneity(
                maf=laml, 
                tsb=as.character(sampleMt[1,3]), 
                vafCol='VAF', useSyn=TRUE)$"clusterData"
            colnames(clusterMt)[6]="VAF"
            ## print VAF pictures for all samples
            picCha <- paste("all", 
                             patientID, 
                             ".", 
                             counterMt, 
                             "<-.vafDraw(
                             clusterMt, 
                             themeOption, 
                             sampleName, 
                             MATHScore, 
                             MIX_option=sampleOption)", 
                             sep="")
            eval(parse(text=picCha))
            picName <- paste("all", patientID, ".", counterMt, sep="")
            lsPicName <- c(lsPicName, picName)
        }
        ## set the columns of the picture and generate all single pictures above
        pic <- eval(parse(text=paste("plot_grid(", 
                                     paste(lsPicName, collapse=","), 
                                     ", nrow=", 
                                     ceiling(length(lsPicName)/2), 
                                     ", ncol=2, align=\"v\")" , sep="")))
        ## save the cowplot picure
        ggsave(pic, filename=paste(sampleName, 
                                   "_VAF_Cluster_MIX", 
                                   ".", 
                                   fileFormat,sep=""), 
               width=12, height=9, dpi=800, path="./output")
    } 
    
    ## one pic for all sample
    else if (sampleOption == "OFA"){
        ## collect all samples' cluster results
        for (counterMt in 1:length(tsbLs[,1])){
            sampleName <- as.character(tsbLs[,1][counterMt])
            ## calculate scoreMATH
            if (showMATH){
                MATHScore <- scoreMATH(mafInput, c(sampleName))
                MATHScore <- MATHScore[which(
                    MATHScore$Tumor_Sample_Barcode == "ITH MATH score"), 
                    ]$MATH_score
            } else {
                MATHScore <- NA
            }
            sampleMt <- vafInputMt[which(
                vafInputMt$Samples %in% sampleName),]
            clusterMtCha <- paste(
                "cluster_mt_", 
                counterMt, 
                " <- inferHeterogeneity(
                maf=laml, 
                tsb=as.character(sampleMt[1,3]), 
                vafCol=\'VAF\', 
                useSyn=TRUE)$\"clusterData\"", 
                sep ="")
            eval(parse(text=clusterMtCha))
            clusterMtCha <- paste("colnames(cluster_mt_", 
                                    counterMt, 
                                    ")[6]=\"VAF\"",sep ="")
            eval(parse(text=clusterMtCha))
            clusterMtCha <- paste("clusterAll <- rbind(
                                    clusterAll, cluster_mt_", 
                                    counterMt, ")",
                                    sep ="")
            eval(parse(text=clusterMtCha))
        }
        colnames(clusterAll)[6]="VAF"
        pic <- eval(parse(text=.vafOFA(clusterAll, 
                                       themeOption, 
                                       tsbLs, 
                                       sampleOption, 
                                       MATHScore)))
        ggsave(pic, filename= paste(patientID, 
                                    "_VAF_Cluster",
                                    ".", 
                                    fileFormat, 
                                    sep=""), 
               width=12, height=9, dpi=1200)
    } 
    
    ## specific sample
    else 
    {
        ## calculate scoreMATH
        if (showMATH){
            MATHScore <- scoreMATH(
                mafInput, c(sampleOption))
            MATHScore <- MATHScore[which(
                MATHScore$Tumor_Sample_Barcode == sampleOption), 
                                     ]$MATH_score
        } else {
            MATHScore <- NA
        }
        sampleMt <- vafInputMt[which(
            vafInputMt$Samples %in% sampleOption),]
        clusterMt <- inferHeterogeneity(maf=laml, tsb=as.character(
                sampleMt[1,3]), vafCol='VAF', useSyn=TRUE)$"clusterData"
        colnames(clusterMt)[6]="VAF"
        pic <- .vafDraw(clusterMt, themeOption, sampleOption, MATHScore)
        ggsave(pic, filename= paste(
            sampleOption,"_VAF_Cluster",".", fileFormat,sep=""), 
            width=12, height=9)
    }
}


############ General Toolbox ############
## VAF main draw vlines for all sampleOption
.vafVline <- function(clusterMt, pic, tsbLs, sampleOption, tsb, 
                      ingredients=NULL){
    ## data prepare
    VAF_vline_cha <- ""
    cluster_ls <- unique(clusterMt$cluster)
    ## density information of the curve for a tsb
    density_info <- data.frame(layer_data(pic))
    
    ## OFA specific: get scaling ratio
    if (typeof(ingredients) == "list"){
        iscale <- ingredients$iscale[1]
        scale <- ingredients$scale[1]
    }
    
    ## Obtain vline Coordinate(x, xend, y, yend)
    for (cluster_name in cluster_ls){
        x_end <- max(clusterMt[which(
            clusterMt$cluster == cluster_name)]$VAF)
        x_end_alter <- density_info$x[which.min(
            abs(outer(density_info$x,x_end,FUN="-")))]
        y_end <- density_info$y[which(
            density_info$x == x_end_alter)]
        if (sampleOption == "OFA"){
            ## Scale and draw lines
            density <- density_info$density[which(
                density_info$x == x_end_alter)]
            VAF_vline_cha <- paste(
                VAF_vline_cha, 
                "geom_segment(data=cluster_mt_", 
                which(tsbLs == tsb), 
                " ,aes(x=", x_end_alter, 
                ", xend=", x_end_alter, 
                ", y=", 
                which(tsbLs == tsb), 
                ", yend=", 
                which(tsbLs == tsb) + density*iscale*scale, 
                "), size=0.5, colour=\"grey\", linetype=\"dashed\") + "
                ,sep="")
        }else{
            VAF_vline_cha <- paste(
                VAF_vline_cha, 
                "geom_segment(aes(x=", 
                x_end_alter, 
                ", xend=", 
                x_end_alter, 
                ", y=0, yend=", 
                y_end, 
                "), size=0.5, 
                colour=\"grey\", linetype=\"dashed\") + ", 
                sep="")
        }
    }
    VAF_vline_cha
}

## Functions for specific sampleOption: "All","tsb"
## VAF painter
.vafDraw <- function(clusterMt, themeOption, sampleOption, MATHScore, 
                     MIX_option=""){
    ## A draft for density infomation(density_info) of ggplot
    picv <- ggplot(clusterMt, aes(x=VAF)) + geom_line(
        size=1, colour="#00C0EB", stat="density")
    if (is.na(MATHScore)){
        if (MIX_option == "MIX"){
            ## generate character/string for ggplot and specific titles for minifigures
            VAF_draw_cha <- paste(
            "ggplot(clusterMt, aes(x=VAF)) 
            + theme_bw() 
            + theme(legend.position=\'none\', 
               title=element_text(size=10), 
               text=element_text(size=10), 
               panel.grid=element_blank(), 
               panel.border=element_blank(), 
               axis.line=element_line(size=0.25)) 
            + geom_line(size=1, colour=\"#00C0EB\", stat=\"density\") 
            + geom_rug(aes(y=0, colour=cluster), sides=\"b\") 
            + ", 
            .vafVline(clusterMt, picv, tsbLs, sampleOption), 
            "scale_color_", 
            themeOption, "() + scale_fill_", 
            themeOption, "()", sep="")
            eval(parse(text=VAF_draw_cha))
        } else {
            ## generate character/string for ggplot and paint the picture
            VAF_draw_cha <- paste(
            "ggplot(clusterMt, aes(x=VAF)) 
            + theme_bw() 
            + theme(title=element_text(size=18), 
            text=element_text(size=18), 
            panel.grid=element_blank(), 
            panel.border=element_blank(), 
            axis.line=element_line(size=0.25)) 
            + geom_line(size=1, colour=\"#00C0EB\", stat=\"density\") 
            + geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
            .vafVline(clusterMt, picv, tsbLs, sampleOption),
            "scale_color_", themeOption, "() + scale_fill_", 
            themeOption, "()", sep="")
            eval(parse(text=VAF_draw_cha))
        } 
    }
    else {
        if (MIX_option == "MIX"){
            ## generate character/string for ggplot and specific titles for minifigures
            VAF_draw_cha <- paste(
            "ggplot(clusterMt, aes(x=VAF)) 
            + theme_bw() 
            + theme(legend.position=\'none\', 
                plot.title=element_text(size=10, 
                hjust=1, vjust=0.5, face='bold'), 
                title=element_text(size=10), 
                text=element_text(size=10), 
                panel.grid=element_blank(), 
                panel.border=element_blank(), 
                axis.line=element_line(size=0.25)) 
            + ggtitle(\"", sampleOption, "\'s MATH Score: ", 
            as.character(MATHScore), "\") 
            + geom_line(size=1, colour=\"#00C0EB\", stat=\"density\")
            + geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
            .vafVline(clusterMt, picv, tsbLs, sampleOption),
            "scale_color_", themeOption, "() + scale_fill_", 
            themeOption, "()", sep="")
            eval(parse(text=VAF_draw_cha))
        } else {
            ## generate character/string for ggplot and paint the picture
            VAF_draw_cha <- paste(
            "ggplot(clusterMt, aes(x=VAF)) 
            + theme_bw() 
            + theme(plot.title=element_text(size=18, 
                hjust=1, vjust=0.5, face='bold'), 
                title=element_text(size=18), text=element_text(size=18), 
                panel.grid=element_blank(),panel.border=element_blank(), 
                axis.line=element_line(size=0.25)) 
            + ggtitle(\"MATH Score: ", as.character(MATHScore), "\") 
            + geom_line(size=1, colour=\"#00C0EB\", stat=\"density\") 
            + geom_rug(aes(y=0, colour=cluster), sides=\"b\") + ", 
            .vafVline(clusterMt, picv, tsbLs, sampleOption),
            "scale_color_", themeOption, "() + scale_fill_", 
            themeOption, "()", sep="")
            eval(parse(text=VAF_draw_cha))
        }
    }
}

## Functions for specific sampleOption: "OFA"
## VAF drawer for OFA: generate character/string for follow-up painting with ggplot. 
.vafOFA <- function(clusterAll, themeOption, tsbLs, sampleOption, MATHScore){
    if (is.na(MATHScore)){
        VAF_ofa_cha <- paste(
        "ggplot(clusterAll, aes(x=VAF, y=Tumor_Sample_Barcode))
        + theme_bw() 
        + theme(title=element_text(size=18), 
            text=element_text(size=18), 
            panel.grid=element_blank(), 
            panel.border=element_blank(), 
            axis.line=element_line(size=0.25)) ",
        "+ geom_density_ridges(
            fill=\"whitesmoke\", 
            calc_ecdf=TRUE, 
            alpha=0.8) ",
        "+ geom_point(aes(x=VAF, y=Tumor_Sample_Barcode, color=cluster), 
            alpha=0.5, show.legend=F) ",
        "+ geom_density_ridges(
            color=\"#00C0EB\", fill=NA, calc_ecdf=TRUE, alpha=0.5) + ",
        .vafVlineOFA(clusterAll, tsbLs, sampleOption), 
        "scale_color_", themeOption, "() + scale_fill_", 
        themeOption, "()", sep="")
    } else {
        VAF_ofa_cha <- paste(
        "ggplot(clusterAll, 
        aes(x=VAF, y=Tumor_Sample_Barcode))
        + theme_bw() 
        + theme(plot.title=element_text(
                size=18, hjust=1, vjust=0.5, face='bold'), 
            title=element_text(size=18), text=element_text(size=18), 
            panel.grid=element_blank(), 
            panel.border=element_blank(), 
            axis.line=element_line(size=0.25)) ", 
        "+ ggtitle(\"MATH Score: ", as.character(MATHScore), "\") ", 
        "+ geom_density_ridges(fill=\"whitesmoke\", 
            calc_ecdf=TRUE, alpha=0.8) ",
        "+ geom_point(aes(x=VAF, y=Tumor_Sample_Barcode, 
            color=cluster), alpha=0.5, show.legend=F) ",
        "+ geom_density_ridges(color=\"#00C0EB\", fill=NA, 
            calc_ecdf=TRUE, alpha=0.5) ",
        .vafVlineOFA(clusterAll, tsbLs, sampleOption), 
        "scale_color_", themeOption, "() + scale_fill_", 
        themeOption, "()", sep="")
    }
    VAF_ofa_cha
}

## VAF draw vlines for ofa
.vafVlineOFA <- function(clusterAll, tsbLs, sampleOption)
{
    ## data prepare
    .vafVlineOFA <- ""
    ## density information of the curve for all tsbs
    gr <- ggplot(clusterAll, aes(
        x=VAF, y=Tumor_Sample_Barcode)) + geom_density_ridges() 
    ingredients <- ggplot_build(gr) %>% purrr::pluck("data", 1)
    
    for (tsb in tsbLs$samples)
    {
        ## renew tsb's data
        VAF_vline_cha <- ""
        x_end_ls <- data.frame()
        y_end_ls <- data.frame()
        clusterMt <- clusterAll[which(
            clusterAll$Tumor_Sample_Barcode == tsb)]
        colnames(clusterMt)[6]="VAF"
        cluster_ls <- unique(clusterMt$cluster)
        ## A draft for density infomation(density_info) of ggplot
        picv <- ggplot(clusterMt, aes(x=VAF)) 
            + geom_line(size=1, colour="#00C0EB", stat="density")
        density_info <- data.frame(layer_data(picv))
        ## collect vlines for a single tsb
        VAF_vline_cha <- .vafVline(clusterMt, picv, tsbLs, 
                                   sampleOption, tsb, ingredients)
        ## collect vlines for all tsbs
        .vafVlineOFA <- paste(.vafVlineOFA, VAF_vline_cha)
    }
    .vafVlineOFA
}


## Directory
setwd("/home/ninomoriaty/R_Project/data/maf")
maf_file1="311252.maf"
maf_file2="313544.maf"
maf_file3="313935.maf"
maf_file4="313953.maf"
maf_file5="314007.maf"
maf_file6="314069.maf"
maf_file7="314155.maf"
sample_info.dir="/home/ninomoriaty/R_Project/data/sample_info.txt"
maf_file_ls=c(maf_file1, maf_file2, maf_file3, maf_file4, 
              maf_file5, maf_file6, maf_file7)
for (counter in maf_file_ls){
    patientID <- strsplit(counter, ".maf")[[1]][1]
    maf <- readMaf(patientID , counter, sample_info.dir)
    ## vafPlot(maf, "OFA", fileFormat="pdf")
    ## vafPlot(maf, "All", fileFormat="pdf")
    vafPlot(maf, "MIX", fileFormat="pdf")
}
