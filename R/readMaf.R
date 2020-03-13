#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile tab delimited MAF file (plain text or *.gz compressed). Required.
#' @param ccfFile CCF file of SNVs. Default NULL.
#' @param mutType select proper variant classification you need. Default "All".Option: "nonSilent".
#' @param mutNonSilent variant classifications which are considered as non-silent. Default NULL.
#' @param chrSilent Select chromosomes needed to be dismissed. Default NULL.
#' @param use.indel logical value. whether to use INDELs besides somatic SNVs. Default FALSE.
#' @param refBuild human reference genome versions of "hg18", "hg19" or "hg38" by UCSC. Default "hg19".
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
    ccfFile = NULL,
    ## filter selection
    mutType = "All",
    mutNonSilent = NULL,
    chrSilent = NULL,
    use.indel = FALSE,
    refBuild = "hg19") {

    ref.options = c('hg18', 'hg19', 'hg38')
    if(!refBuild %in% ref.options){
        stop("refBuild can only be either 'hg18', 'hg19' or 'hg38'")
    }

    mutType.options = c("All", "nonSilent")
    if(!mutType %in% mutType.options){
        stop("mutType should be either 'All' or 'nonSilent.")
    }

    ## get patientID
    # patientID <- unlist(strsplit(basename(mafFile), split = "[.]"))[1]

    mafData <- data.table::fread(
            file = mafFile,
            quote = "",
            header = TRUE,
            fill = TRUE,
            sep = '\t',
            skip = "Hugo_Symbol",
            stringsAsFactors = FALSE
        )
    

    
    ## read ccf files
    if (!is.null(ccfFile)) {
        ccfInput <- suppressWarnings(read.table(
            ccfFile,
            quote = "",
            header = TRUE,
            fill = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
        ))
        
        mafData <- mafData %>%
            uniteCCF(ccfInput) %>%
            getMutStatus() %>%
            dplyr::mutate(VAF_adj = CCF/2) ## calculate adjusted VAF based on CCF
    }
    
    ## filter variant classification
    if (mutType == "nonSilent") {
        if (is.null(mutNonSilent)) {
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
    }
    
    ## use.indel filter
    if (!use.indel) {
        mafData <- mafData[which(mafData$Variant_Type == "SNP"), ]
    }
    
    ## chromosome filter
    if (!is.null(chrSilent)) {
        mafData <- mafData[which(!mafData$Chromosome %in% chrSilent),]
    }
    
    ## Rescale vaf coloum 0-1
    if(max(mafData$VAF, na.rm = TRUE) > 1){
        mafData$VAF <- as.numeric(as.character(mafData$VAFdat))/100
    }    
    
    ## generate classMaf
    maf <- classMaf(
        data = data.table::setDT(mafData),
        #patientID = patientID,
        ref.build = refBuild
    )
    
    # ## for parameter vafColumn="VAF", select particular VAF column
    # colnames(maf@data)[colnames(maf@data) == vafColumn] <- "VAF"
    
    return(maf)
}


##--- combine CCF into maf object
uniteCCF <- function(mafData, ccf) {
    mafData <- tidyr::unite(
        mafData,
        "mutID",
        c(
            "Patient_ID",
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
                "Patient_ID",
                "Tumor_Sample_Barcode",
                "Chromosome",
                "Start_Position",
                "Variant_Type"
            ),
            sep = ":",
            remove = FALSE
        ) %>%
        dplyr::select(mutID, CCF, CCF_std)
    
    mafData_merge_ccf <-
        merge(mafData, ccf, by = "mutID", all.x = TRUE) %>%
        dplyr::select(-mutID)
}




getMutStatus <- function(mafData, ccf.conf.level = 0.95) {
    mafData <-
        mafData %>%
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
        ref.build = 'character'
    )
)
