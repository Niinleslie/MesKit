#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile tab delimited MAF file (plain text or *.gz compressed). The prefix to the filename of mafFile is considered as "patientID", which will be used for downstream analysis and visualization. Required.
#' @param mutType Select proper variant classification you need. Default "All".Option: "nonSilent".
#' @param mutNonSilent variant classifications which are considered as non-silent. Default NULL. Option: "Default".
#' @param chrSilent Select chromosomes needed to be dismissed. Default NULL.
#' @param use.indel logic. whether to use INDELs besides somatic SNVs. Default FALSE.
#' @param ccfFile CCF file of SNVs. Default NULL.
#' @param refBuild human reference genome versions of hg18, hg19 or hg38 by UCSC. Default "hg19".
#'
#'
#' @examples
#' maf.File <- system.file("extdata/maf", "HCC6046.maf", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC6046.CCF.txt", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, refBuild="hg19")
#' maf <- readMaf(mafFile=maf.File, ccfFile=ccf.File, refBuild="hg19")
#' @return an object of class Maf.
#'
#' @exportClass classMaf
#' @export readMaf


## read.maf main function
readMaf <- function(## maf parameters
    mafFile,
    ## filter selection
    mutType = "All",
    mutNonSilent = "Default",
    chrSilent = NULL,
    use.indel = FALSE,
    ## ccf parameters
    ccfFile = NULL,
    refBuild = "hg19") {

    ref.options = c('hg18', 'hg19', 'hg38')
    if(!refBuild %in% ref.options){
        stop("refBuild can only be either 'hg18', 'hg19' or 'hg38'")
    }


    mafData <- data.table::fread(
            file = mafFile,
            quote = "",
            header = TRUE,
            fill = TRUE,
            sep = '\t',
            skip = "Hugo_Symbol",
            stringsAsFactors = FALSE
        )
    
    ## get patientID
    patientID <- unlist(strsplit(basename(mafFile), split = "[.]"))[1]
    
    ## read ccf files
    if (!is.null(ccfFile)) {
        ccfInput <- read.table(
            ccfFile,
            quote = "",
            header = TRUE,
            fill = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
        )
        
        mafData <- mafData %>%
            uniteCCF(ccfInput) %>%
            getMutStatus() %>%
            dplyr::mutate(VAF_adj = CCF/2) ## calculate adjusted VAF based on CCF
    }
    
    ## filter variant classification
    if (mutType == "nonSilent") {
        if (mutNonSilent == "Default") {
            nonSilent <- c(
                "Frame_Shift_Del",
                "Frame_Shift_Ins",
                "Splice_Site",
                "Translation_Start_Site",
                "Nonsense_Mutation",
                "Nonstop_Mutation",
                "In_Frame_Del",
                "In_Frame_Ins",
                "Missense_Mutation"
            )
        } else {
            nonSilent <- mutNonSilent
        }
        mafData <-
            mafData[which(mafData$Variant_Classification %in% nonSilent),]
    } else if (mutType == "All") {
        # message("All variant classification submitted")
    } else {
        error(
            "parameter `mut.type` error.
              The mut.type should be either 'All' or 'nonSilent'.
              You could further settle the filter by parameter 'mutNonSilent'."
        )
    }
    
    ## use.indel filter
    if (!use.indel) {
        mafData <- mafData[which(mafData$Variant_Type == "SNP"), ]
    }
    
    ## chromosome filter
    if (!is.null(chrSilent)) {
        mafData <- mafData[which(!mafData$Chromosome %in% chrSilent),]
    }
    
    
    ## generate classMaf
    maf <- classMaf(
        data = data.table::setDT(mafData),
        patientID = patientID,
        ref.build = refBuild
    )
    
    # ## for parameter vafColumn="VAF", select particular VAF column
    # colnames(maf@data)[colnames(maf@data) == vafColumn] <- "VAF"
    
    return(maf)
}


##--- combine CCF into maf object
uniteCCF <- function(mafdata, ccf) {
    mafdata <- tidyr::unite(
        mafdata,
        "mutID",
        c(
            "Tumor_Sample_Barcode",
            "Chromosome",
            "Start_Position",
            "Variant_Type"
        ),
        sep = ":",
        remove = FALSE
    )
    ccf <- ccf %>%
        dplyr::mutate(Variant_Type = "SNP") %>%
        tidyr::unite(
            "mutID",
            c(
                "Sample_ID",
                "Chromosome",
                "Start_Position",
                "Variant_Type"
            ),
            sep = ":",
            remove = FALSE
        ) %>%
        dplyr::select(mutID, CCF, CCF_std)
    
    mafdata_merge_ccf <-
        merge(mafdata, ccf, by = "mutID", all.x = TRUE) %>%
        dplyr::select(-mutID)
}




getMutStatus <- function(mafdata, ccf.conf.level = 0.95) {
    mafdata <-
        mafdata %>%
        # 95% confidence interval
        # normal distribution
        dplyr::mutate(CCF_max = CCF +
                          qnorm((1 - ccf.conf.level) / 2, lower.tail = FALSE) * CCF_std) %>%
        dplyr::mutate(Status =
                          dplyr::case_when(CCF_max >= 1 ~ "Clonal",
                                           CCF_max < 1 ~ "Subclonal")) %>%
        dplyr::select(-CCF_max)
}

##--- classMaf class
classMaf <- setClass(
    Class = "classMaf",
    slots = c(
        data = 'data.table',
        patientID = 'character',
        ref.build = 'character'
    )
)
