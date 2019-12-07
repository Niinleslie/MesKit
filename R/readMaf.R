#' @description Read tab delimited MAF (can be plain text or gz compressed) file along with sample information file.
#'
#' @import ggplot2
#'
#' @param mafFile MAF file. 
#' @param sampleInfoFile sample information file.
#' @param mutType Default "All". "nonSilent" And you can select proper variant classification you need. 
#' @param mutNonSilent Default NULL. Option: "Default".And you can list variant classifications that you do not want them to be silent.
#' @param chrSilent Default NULL. select the chromosomes you want to dismiss.
#' @param use.indel Select SNP in Variant type
#' @param ccfClusterTsvFile CCF cluster.tsv file if ccf data provided. Default NULL.
#' @param ccfLociTsvFile CCF loci.tsv file if ccf data provided. Default NULL.
#' @param refBuild Default "hg19". You could choose human reference genome versions of hg19 or hg38 by UCSC.
#' @param name Default NULL. Correct name error of shiny fileinput
#' 
#' @return a Maf object/class 
#' 
#' @examples
#' ## data information
#' maf.File <- system.file("extdata/maf", "HCC6046.maf", package = "Meskit")
#' sampleInfo.File <- system.file("extdata", "HCC6046.sampleInfo.txt", package = "Meskit")
#' pyCloneCluster <- system.file("extdata/ccf", "HCC6046.cluster.tsv", package = "Meskit")
#' pyCloneLoci <- system.file("extdata/ccf", "HCC6046.loci.tsv", package = "Meskit")
#' ## manually usage
#' maf <- readMaf(mafFile=maf.File, sampleInfoFile=sampleInfo.File, refBuild="hg19")
#' ## if ccf data provided
#' maf <- readMaf(mafFile=maf.File, sampleInfoFile=sampleInfo.File, ccfClusterTsvFile=pyCloneCluster, ccfLociTsvFile=pyCloneLoci, refBuild="hg19")
#' 
#' @exportClass classMaf
#' @export readMaf


## read.maf main function
readMaf <- function(
    ## maf parameters
    mafFile, sampleInfoFile,
    ## filter selection
    mutType="All", mutNonSilent=NULL, chrSilent=NULL, use.indel=FALSE, 
    ## ccf parameters             
    ccfClusterTsvFile=NULL, ccfLociTsvFile=NULL, 
    ## supplyment
    refBuild="hg19",
    ## supply in shiny
    name=NULL){
    
    ## read maf file
    if (.substrRight(mafFile, 3) == ".gz"){
        mafInput <- read.table(mafGz <- gzfile(mafFile, "r"), quote="", 
                               header=TRUE, fill=TRUE, 
                               sep='\t', stringsAsFactors=FALSE)
        
        close(mafGz)
    } else {
        mafInput <- read.table(mafFile, quote="", 
                               header=TRUE, fill=TRUE, 
                               sep='\t', stringsAsFactors=FALSE)
    }
    
    ## get patientID
    fileName <- unlist(strsplit(mafFile, "/"))[length(unlist(strsplit(mafFile, "/")))]
    patientID <- strsplit(as.character(fileName), ".maf")[[1]][1]
    ## correct error name of fileinput on shiny app
    if(patientID == 0){
      patientID <- strsplit(name,"\\.")[[1]][1]
    }
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
    #mafInput$time <- ""
    
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
        #time <- as.character(
        #    sampleInfoInput[which(
        #        sampleInfoInput$sample == tsb),]$time)
        mafInput[which(
            mafInput$Tumor_Sample_Barcode == tsb),]$patient <- patient
        mafInput[which(
            mafInput$Tumor_Sample_Barcode == tsb),]$lesion <- lesion
        #mafInput[which(
        #    mafInput$Tumor_Sample_Barcode == tsb),]$time <- time
    }
    
    # ## fix: Error in setattr(x, "row.names", rn)
    # mafInput$Hugo_Symbol <- as.character(mafInput$Hugo_Symbol)
    
    ## filter variant classification
    if (mutType == "nonSilent"){
        if (mutNonSilent == "Default"){
            nonSilent <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                          "Translation_Start_Site", "Nonsense_Mutation", 
                          "Nonstop_Mutation", "In_Frame_Del",
                          "In_Frame_Ins", "Missense_Mutation")
        } else {
            nonSilent <- mutNonSilent 
        }
        mafInput <- mafInput[which(mafInput$Variant_Classification %in% nonSilent), ]
    } else if (mutType == "All"){
        # message("All variant classification submitted")
    } else {
        error("mut.type setting error")
    }
    
    ## use.indel filter
    if(!use.indel){
        mafInput <- mafInput[which(mafInput$Variant_Type == "SNP"),]
    }
    
    ## filter chromosome
    if (!is.null(chrSilent)){
        mafInput <- mafInput[which(!mafInput$Chromosome %in% chrSilent), ]
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
    
    # ## , vafColumn="VAF", 
    # colnames(maf@data)[colnames(maf@data) == vafColumn] <- "VAF"
    

    
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
