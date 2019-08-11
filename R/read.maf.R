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
#' \dontrun{
#' read.maf(patientID, dat.dir = "./data", use.ccf = FALSE, plot.mafSummary = TRUE) 
#' }

## classMaf class
classMaf <- setClass(Class="classMaf", contains="MAF", 
                slots= c(ccf.cluster='data.table', ccf.loci='data.table', 
                         patientID='character', ref.build='character'))

## read.maf main function
readMaf <- function(patientID, mafDir, 
                    sampleInfoDir, ccfDir=NULL, 
                    plotMafSummary=TRUE, refBuild="hg19"){
    ## read maf file
    mafInput <- read.table(mafDir, quot="", 
                            header=TRUE, fill=TRUE, 
                            sep='\t')
    ## read info file
    sampleInfoInput <-  read.table(sampleInfoDir, quot="", 
                                     header=TRUE, fill=TRUE, 
                                     sep='', stringsAsFactors=FALSE)
    ## read ccf file
    if (!is.null(ccfDir)) {
        ccfClusterTsvInput <- read.table(ccfDir, quote="", 
                                            header=TRUE, fill=TRUE, 
                                            sep='\t', stringsAsFactors=FALSE)
        ccfLociTsvInput <- read.table(ccfDir, quote="", 
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
    mafSum <- read.maf(mafData)
    
    # mafSum2 <- suppressMessages(.summarizeMaf(mafData))
    maf <- classMaf(data=mafData, 
               variants.per.sample=mafSum@variants.per.sample, 
               variant.type.summary=mafSum@variant.type.summary,
               variant.classification.summary=mafSum@
                   variant.classification.summary, 
               gene.summary=mafSum@gene.summary,
               summary=mafSum@summary, 
               maf.silent=mafSilent, 
               clinical.data=mafSum@clinical.data, 
               ccf.cluster=ccfClusterTsv, 
               ccf.loci=ccfLociTsv, 
               patientID=patientID, 
               ref.build=refBuild)
    
    ## print the summary plot
    if (plotMafSummary) {
        pic <- plotmafSummary(maf=maf, rmOutlier=TRUE, 
                              addStat='median', dashboard=TRUE, 
                              titvRaw=FALSE)
        ggsave(plotmafSummary(maf=maf, rmOutlier=TRUE, addStat='median', 
                              dashboard=T, titvRaw=FALSE), 
               filename=paste(patientID, ".VariantSummary.png", sep=""), 
               width=12, height=9, dpi=800, path="./output")
    }
    message("Class Maf Generation Done!")
    return(maf)
}
