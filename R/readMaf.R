#' Output an new maf document with information of sample_info 
#' @description Add sample_info to the original maf file to get a new maf.The
#'  new maf file adds three pieces of information:lesion,patient and time
#'
#' @import ggplot2 
#' @importClassesFrom maftools MAF
#' @importFrom maftools plotmafSummary
#' @importFrom maftools read.maf
#'
#' @param patientID patient/sample name
#' @param dat.dir specify a data directory as the input of the function
#' @param use.ccf=Default FALSE. TRUE if ccf data is provided
#' @param plot.mafSummary Default TRUE. FALSE if the summary figure is 
#' unnecessary
#' @param ref.build Default hg19. specify a referential genome including 
#' BSgenome.Hsapiens.UCSC.hg19 and BSgenome.Hsapiens.UCSC.hg38 
#' @return a classMaf object/class includes information of sample_info and 
#' mut.id and summary figure of it
#' 
#' @exportClass classMaf
#' @export readMaf
#'
#' @examples
#' maf.File <- system.file("extdata/multi_lesion/maf", "311252.maf", package = "Meskit")
#' sampleInfo.File <- system.file("extdata/multi_lesion", "sample_info.txt", package = "Meskit")
#' pyCloneCluster <- system.file("extdata/multi_lesion/ccf", "311252.cluster.tsv", package = "Meskit")
#' pyCloneLoci <- system.file("extdata/multi_lesion/ccf", "311252.loci.tsv", package = "Meskit")
#' maf <- readMaf(patientID = "311252", mafFile = maf.File, sampleInfo = sampleInfo.File, refBuild = "hg19")
#' # if use ccf 
#' maf <- readMaf(patientID = "311252", mafFile  = maf.File, sampleInfo = sampleInfo.File,
#'                 ccfClusterTsvDir = pyCloneCluster, ccfLociTsvInput = pyCloneLoci, refBuild = "hg19")

## classMaf class
classMaf <- setClass(Class="classMaf", contains="MAF", 
                     slots= c(ccf.cluster='data.table', ccf.loci='data.table', 
                              patientID='character', ref.build='character'))

## read.maf main function
readMaf <- function(patientID, mafFile, 
                    sampleInfo, ccfClusterTsvDir=NULL, ccfLociTsvInput=NULL, 
                    refBuild="hg19", plotMafSummary=TRUE, savePlot = FALSE, summaryPlot.dir = NULL){
    ## read maf file
    mafInput <- read.table(mafFile, quote="", 
                           header=TRUE, fill=TRUE, 
                           sep='\t')
    ## read info file
    sampleInfoInput <-  read.table(sampleInfo, quote="", 
                                   header=TRUE, fill=TRUE, 
                                   sep='', stringsAsFactors=FALSE)
    ## read ccf file
    if (!is.null(ccfClusterTsvDir) & !is.null(ccfLociTsvInput)) {
        ccfClusterTsvInput <- read.table(ccfClusterTsvDir, quote="", 
                                         header=TRUE, fill=TRUE, 
                                         sep='\t', stringsAsFactors=FALSE)
        ccfLociTsvInput <- read.table(ccfLociTsvInput, quote="", 
                                      header=TRUE, fill=TRUE, 
                                      sep='\t', stringsAsFactors=FALSE)
    } else {
        ccfClusterTsvInput <- NULL
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
    ccfClusterTsv <- data.table::setDT(ccfClusterTsvInput)
    ccfLociTsv <- data.table::setDT(ccfLociTsvInput)
    
    ## summarize sample_info and mut.id with summarizeMaf
    mafSum <- maftools::read.maf(mafData)
    
    # mafSum2 <- suppressMessages(.summarizeMaf(mafData))
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
    if (plotMafSummary) {
        pic <- maftools::plotmafSummary(maf=maf, rmOutlier=TRUE, 
                              addStat='median', dashboard=TRUE, 
                              titvRaw=FALSE)
        if (savePlot){
            ggsave(maftools::plotmafSummary(maf=maf, rmOutlier=TRUE, addStat='median', 
                              dashboard=T, titvRaw=FALSE), 
               filename=paste(patientID, ".VariantSummary.png", sep=""), 
               width=12, height=9, dpi=800, path=summaryPlot.dir)

        }
    }
    message("Class Maf Generation Done!")
    return(maf)
}
