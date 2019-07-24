#' A Multifunctional MATH Calculator
#' @description read maf file as data.frame, calculate MATH score and
#' present results in different ways determined by parameters. This function 
#' requires VAF for MATH score calculation and Tumor_sample_Barcode for 
#' selection.
#' VAF can be on the scale 0-1 and you could filter VAF with minvaf and maxvaf.
#' 
#' @param maf_file specify a maf document/directory as the input of the function
#' @param tsb specify single/all sample names (Tumor_Sample_Barcodes)
#' @param minvaf filter low frequency variants caused by sequencing error. 
#' Default 0. (on the scale of 0 to 1)
#' @param maxvaf filter high frequency variants due to copy number alterations 
#' or impure tumor. Default 1. (on the scale of 0 to 1)
#' @return MATH score list for all/selected samples or a MATH score value.

#' @examples
#' \dontrun{
#' maf_file="/home/ninomoriaty/R_Project/311252_snv_indel.imputed.maf"
#' (maf_file)
#' (maf_file, tsb=c("tsb1"))
#' (maf_file, tsb=c("tsb1", "tsb2", "tsb3"))
#'}


## MATH Score main function
ScoreMATH <- function(mafData, tsb=c("OFA"), minvaf=0, maxvaf=1){
    ## get vaf-related infomation
    dataHugoSymbol <- mafData$Hugo_Symbol
    dataVaf <- mafData$VAF
    dataTsb <- data.frame(mafData$Tumor_Sample_Barcode)
    vafInputMt <- data.frame(dataHugoSymbol, dataVaf, dataTsb)
    colnames(vafInputMt) <- c("Hugo_Symbol", "VAF", "Tumor_Sample_Barcode")
    ## get all sample names
    tsbLs <- data.frame(unique(vafInputMt$Tumor_Sample_Barcode))
    
    ## MATH results for one/all sample
    if (any(tsb == c("OFA"))){
        ## list all samples' MATH scores
        mathOFA <- .multiSampleMATH(vafInputMt, tsbLs, minvaf, maxvaf)
        mathAll <- .patientMATH(vafInputMt, minvaf, maxvaf)
        mathAll <- data.frame(Tumor_Sample_Barcode=c("ITH MATH score"), 
                               MATH_score=c(mathAll))
        return(rbind(mathOFA, mathAll))
    } else{
        ## calculate specific samples' MATH score
        tsbLs <- data.frame(tsb)
        mathSp <- .multiSampleMATH(vafInputMt, tsbLs, minvaf, maxvaf)
        mathAll <- .patientMATH(vafInputMt, minvaf, maxvaf)
        mathAll <- data.frame(Tumor_Sample_Barcode=c("ITH MATH score"), 
                               MATH_score=c(mathAll))
        return(rbind(mathSp, mathAll))
    }
}


## Data cleaning
.dataClean <- function(vafInputMt, tsb, minvaf, maxvaf){
    vafColumn <- vafInputMt[which(
        vafInputMt$Tumor_Sample_Barcode == tsb), ]$VAF
    vafColumn <- vafColumn[which(
        !is.na(vafColumn))][which(
            vafColumn > minvaf & vafColumn < maxvaf)]
    vafColumn <- as.numeric(
        as.character(vafColumn))[which(!is.na(vafColumn))]
    vafColumn
}

## MATH Caculation
.calMATH <- function(vafColumn){
    MAD <- 1.4826*median(abs(vafColumn - median(vafColumn)))
    MATH <- 100 * MAD / median(vafColumn)
    MATH
}

## MATH multi-sample process
.multiSampleMATH <- function(vafInputMt, tsb_ls, minvaf, maxvaf){
    samplesMATH <- data.frame()
    for (counter in seq_along(tsb_ls[,1])){
        for (sampleNameMt in tsb_ls){
            vafColumn <- .dataClean(
                vafInputMt, 
                as.character(sampleNameMt)[counter], 
                minvaf, maxvaf)
            sampleMATH <- data.frame(
                as.character(sampleNameMt)[counter], 
                .calMATH(vafColumn))
            samplesMATH <- rbind(samplesMATH, sampleMATH)
        }
    }
    colnames(samplesMATH) <- c("Tumor_Sample_Barcode", "MATH_score")
    samplesMATH
}

## MATH patient calcualtion
.patientMATH <- function(vafInputMt, minvaf, maxvaf){
    vafColumn <- vafInputMt$VAF
    vafColumn <- vafColumn[which(
        !is.na(vafColumn))][which(
            vafColumn > minvaf & vafColumn < maxvaf)]
    vafColumn <- as.numeric(
        as.character(vafColumn))[which(
            !is.na(vafColumn))]
    .calMATH(vafColumn)
}
