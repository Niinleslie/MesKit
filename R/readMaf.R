#' MAF/CCF/Sample_info Reader
#' @description Add sample information and ccf data to the original maf file 
#' to get a new maf. The new maf data adds three pieces of information:lesion, 
#' patient and time.
#'
#' @import ggplot2
#' @importClassesFrom maftools MAF
#' @importFrom maftools plotmafSummary
#' @importFrom maftools read.maf
#'
#' @param patientID patientID for all samples in mafFile.
#' @param mafFile MAF file directory. 
#' @param sampleInfoFile sample_info.txt file directory.
#' @param ccfClusterFile CCF cluster.tsv file directory if ccf data provided. Default NULL.
#' @param ccfLociTsvFile CCF loci.tsv file directory if ccf data provided. Default NULL.
#' @param refBuild BSgenome.Hsapiens.UCSC reference. Default "hg19". Full genome sequences for Homo sapiens (Human) as provided by UCSC.
#' @param MafSummary Option for whether printing a MafSummary plot or not. Default TRUE.
#' @param outputDir Directory for ouput files. Default NULL.
#' @param gzOption Option for whether generating a .maf.gz compressed file. Default "". If a .maf.gz file is needed, the gzOption could be set as "gz"
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
#' maf <- readMaf(patientID="311252", mafFile=maf.File, sampleInfoFile=sampleInfo.File, refBuild="hg19")
#' ## if ccf data provided
#' maf <- readMaf(patientID="311252", mafFile=maf.File, sampleInfoFile=sampleInfo.File, ccfClusterFile=pyCloneCluster, ccfLociTsvFile=pyCloneLoci, refBuild="hg19")
#' 
#' @exportClass classMaf
#' @export readMaf
#'

## classMaf class
classMaf <- setClass(Class="classMaf", contains="MAF", 
                     slots= c(ccf.cluster='data.table', ccf.loci='data.table', 
                              patientID='character', ref.build='character'))

## read.maf main function
readMaf <- function(patientID, mafFile, 
                    sampleInfoFile, ccfClusterFile=NULL, 
                    ccfLociTsvFile=NULL, refBuild="hg19", 
                    MafSummary=TRUE, outputDir=NULL, gzOption=""){
    
    ## read maf file
    mafInput <- read.table(mafFile, quote="", 
                           header=TRUE, fill=TRUE, 
                           sep='\t')
    
    ## save .maf file as gzip file (need to be confirmed)
    if (gzOption == "gz"){
        oriDir <- getwd()
        setwd(outputDir)
        mafGzName <- paste(patientID, ".maf.gz", sep="")
        mafGz <- gzfile(mafGzName, "w")
        write.table(mafInput, mafGz, quote=FALSE, sep='\t', row.names=FALSE)
        close(mafGz)
        message(paste(mafGzName," generation Done!", sep=""))
        setwd(oriDir)
    }
    
    ## read sample_info file
    sampleInfoInput <-  read.table(sampleInfoFile, quote="", 
                                   header=TRUE, fill=TRUE, 
                                   sep='', stringsAsFactors=FALSE)
    ## read ccf file
    if (!is.null(ccfClusterFile) & !is.null(ccfLociTsvFile)) {
        ccfClusterInput <- read.table(ccfClusterFile, quote="", 
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
    
    ## transform data.frame to data.table
    mafData <- data.table::setDT(mafInput)
    ccfClusterTsv <- data.table::setDT(ccfClusterInput)
    ccfLociTsv <- data.table::setDT(ccfLociTsvInput)
    
    ## summarize sample_info and mut.id with summarizeMaf
    mafSum <- suppressMessages(read.maf(mafData))
    # mafSum2 <- suppressMessages(.summarizeMaf(mafData))
    
    ## generate classMaf
    maf <- classMaf(data=mafData, 
                    variants.per.sample=mafSum@variants.per.sample, 
                    variant.type.summary=mafSum@variant.type.summary,
                    variant.classification.summary=mafSum@
                        variant.classification.summary, 
                    gene.summary=mafSum@gene.summary,
                    summary=mafSum@summary, 
                    maf.silent=mafSum@maf.silent, 
                    clinical.data=mafSum@clinical.data, 
                    ccf.cluster=ccfClusterTsv, 
                    ccf.loci=ccfLociTsv, 
                    patientID=patientID, 
                    ref.build=refBuild)
    
    ## print the summary plot
    if (MafSummary) {
        pic <- plotmafSummary(maf=maf, rmOutlier=TRUE, 
                              addStat='median', dashboard=TRUE, 
                              titvRaw=FALSE)
        
        ## check the output directory
        if (is.null(outputDir)){
            warning("NOTE: Missing output directory for pictures")
        } else {
            ## save the output of MafSummary plot.
            ggsave(pic, 
                   filename=paste(patientID, ".VariantSummary.pdf", sep=""), 
                   width=12, height=9, dpi=1200, path=outputDir)
            message("VariantSummary Plot Saved!")
        }
    }
    message("Class Maf Generation Done!")
    return(maf)
}
