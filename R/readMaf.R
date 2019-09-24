#' MAF/CCF/Sample_info Reader
#' @description Add sample information and ccf data to the original maf file 
#' to get a new maf. The new maf data adds three pieces of information:lesion, 
#' patient and time.
#'
#' @import ggplot2
#'
#' @param mafFile MAF file directory. 
#' @param sampleInfoFile sample_info.txt file directory.
#' @param mut.type select proper variant classification you need
#' @param ccfClusterTsvFile CCF cluster.tsv file directory if ccf data provided. Default NULL.
#' @param ccfLociTsvFile CCF loci.tsv file directory if ccf data provided. Default NULL.
#' @param refBuild BSgenome.Hsapiens.UCSC reference. Default "hg19". Full genome sequences for Homo sapiens (Human) as provided by UCSC.
#' 
#' @return a classMaf object/class includes information of sample_info and 
#' mut.id and summary figure of it
#' 
#' @examples
#' ## data information
#' maf.File <- system.file("extdata/multi_lesion/maf", "311252.maf", package = "Meskit")
#' sampleInfo.File <- system.file("extdata/multi_lesion", "sample_info.txt", package = "Meskit")
#' pyCloneCluster <- system.file("extdata/multi_lesion/ccf", "311252.cluster.tsv", package = "Meskit")
#' pyCloneLoci <- system.file("extdata/multi_lesion/ccf", "311252.loci.tsv", package = "Meskit")
#' ## manually usage
#' maf <- readMaf(mafFile=maf.File, sampleInfoFile=sampleInfo.File, refBuild="hg19")
#' ## if ccf data provided
#' maf <- readMaf(mafFile=maf.File, sampleInfoFile=sampleInfo.File, ccfClusterTsvFile=pyCloneCluster, ccfLociTsvFile=pyCloneLoci, refBuild="hg19")
#' 
#' @exportClass classMaf
#' @export readMaf

## read.maf main function
readMaf <- function(mafFile, sampleInfoFile, 
                    mut.type="All", mut.nonSilent=NULL, 
                    ccfClusterTsvFile=NULL, ccfLociTsvFile=NULL, 
                    refBuild="hg19"){
    
    ## read maf file
    if (.substrRight(mafFile, 3) == ".gz"){
        mafInput <- read.table(mafGz <- gzfile(mafFile, "r"), quote="", 
                               header=TRUE, fill=TRUE, 
                               sep='\t')
        close(mafGz)
    } else {
        mafInput <- read.table(mafFile, quote="", 
                               header=TRUE, fill=TRUE, 
                               sep='\t')
    }
    
    ## get patientID
    fileName <- unlist(strsplit(mafFile, "/"))[length(unlist(strsplit(mafFile, "/")))]
    patientID <- strsplit(fileName, ".maf")[[1]][1]
    
    ## read sample_info file
    sampleInfoInput <-  read.table(sampleInfoFile, quote="", 
                                   header=TRUE, fill=TRUE, 
                                   sep='', stringsAsFactors=FALSE)
    ## read ccf file
    if (!is.null(ccfClusterTsvFile) & !is.null(ccfLociTsvFile)) {
        ccfClusterInput <- read.table(ccfClusterTsvFile, quote="", 
                                      header=TRUE, fill=TRUE, 
                                      sep='\t', stringsAsFactors=FALSE)
        ccfLociTsvInput <- read.table(ccfLociTsvFile, quote="", 
                                      header=TRUE, fill=TRUE, 
                                      sep='\t', stringsAsFactors=FALSE)
    } else {
        ccfClusterInput <- NULL
        ccfLociTsvInput <- NULL
    }
    ## Generate patient,lesion and time information
    mafInput$patient <- ""
    mafInput$lesion <- ""
    mafInput$time <- ""
    
    ## combine sample_info_input with maf_input
    tsbSampleInfo <- unique(sampleInfoInput$sample)
    tsbLs <- tsbSampleInfo[which(
        tsbSampleInfo %in% mafInput$Tumor_Sample_Barcode)]
    for (tsb in tsbLs) {
        patient <- as.character(
            sampleInfoInput[which(
                sampleInfoInput$sample == tsb),]$patient) 
        lesion <- as.character(
            sampleInfoInput[which(
                sampleInfoInput$sample == tsb),]$lesion)
        time <- as.character(
            sampleInfoInput[which(
                sampleInfoInput$sample == tsb),]$time)
        mafInput[which(
            mafInput$Tumor_Sample_Barcode == tsb),]$patient <- patient
        mafInput[which(
            mafInput$Tumor_Sample_Barcode == tsb),]$lesion <- lesion
        mafInput[which(
            mafInput$Tumor_Sample_Barcode == tsb),]$time <- time
    }
    ## fix: Error in setattr(x, "row.names", rn)
    mafInput$Hugo_Symbol <- as.character(mafInput$Hugo_Symbol)
    if (mut.type == "nonSilent"){
        if (is.null(mut.nonSilent)){
            nonSilent <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                          "Translation_Start_Site", "Nonsense_Mutation", 
                          "Nonstop_Mutation", "In_Frame_Del",
                          "In_Frame_Ins", "Missense_Mutation")
        } else {
            nonSilent <- mut.nonSilent 
        }
        mafInput = mafInput[which(mafInput$Variant_Classification %in% nonSilent), ]
    } else if (mut.type == "All"){
        # message("All variant classification submitted")
    } else {
        error("mut.type setting error")
    }
    
    ## transform data.frame to data.table
    mafData <- data.table::setDT(mafInput)
    ccfClusterTsv <- data.table::setDT(ccfClusterInput)
    ccfLociTsv <- data.table::setDT(ccfLociTsvInput)
    
    ## generate classMaf
    maf <- classMaf(data=mafData, 
                    ccf.cluster=ccfClusterTsv, 
                    ccf.loci=ccfLociTsv, 
                    patientID=patientID, 
                    ref.build=refBuild)
    
    return(maf)
}

.substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
}

## classMaf class
classMaf <- setClass(Class="classMaf", 
                     slots= c(data='data.table', ccf.cluster='data.table', 
                              ccf.loci='data.table', patientID='character', 
                              ref.build='character'))
