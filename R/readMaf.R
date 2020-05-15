#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile Tab delimited MAF file (plain text or *.gz compressed). Required.
#' @param ccfFile CCF file of somatic mutations. Default NULL.
#' @param mutType select Proper variant classification you need. Default "All". Option: "nonSilent".
#' @param mutNonSilent Variant classifications which are considered as non-silent. Default NULL.
#' @param ccf.conf.level The confidence level of CCF to identify clonal or subclonal. Only works when "CCF_std" or "CCF_CI_high" is provided in ccfFile. Default: 0.95
#' @param refBuild Human reference genome versions of "hg18", "hg19" or "hg38" by UCSC. Default: "hg19".
#'
#'
#' @examples
#' maf.File <- system.file("extdata/", "HCC6046.maf", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC6046.CCF.txt", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, refBuild="hg19")
#' maf <- readMaf(mafFile=maf.File, ccfFile=ccf.File, refBuild="hg19")
#' @return an object of class Maf.
#'
#' @exportClass classMaf
#' @export readMaf


## read.maf main function
readMaf <- function(
    mafFile,
    ccfFile = NULL,
    adjusted.VAF = FALSE,
    mutNonSilent = NULL,
    ccf.conf.level = 0.95,
    refBuild = "hg19") {

    refBuild <- match.arg(refBuild, choices =  c('hg18', 'hg19', 'hg38'), several.ok = FALSE)
    
    ## get non-silent muation types
    if (is.null(mutNonSilent)) {
        mutNonSilent <- c(
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
    }

    mafData <- data.table::fread(
            file = mafFile,
            quote = "",
            header = TRUE,
            data.table = TRUE,
            fill = TRUE,
            sep = '\t',
            skip = "Hugo_Symbol",
            stringsAsFactors = FALSE
        )
    
    ## check maf data
    mafData <- validMaf(mafData)
    
    ## get mutation id:
    mafData <- mafData %>% 
        tidyr::unite(
            "Mut_ID",
            c(
                "Patient_ID",
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2"
            ),
            sep = ":",
            remove = FALSE
        )
    
    if(adjusted.VAF){
        mafData$VAF_adj <- mafData$VAF
    } 

    # Add sampleinfo
    patients.dat <- split(mafData, mafData$Patient_ID)
    sample.info <- lapply(patients.dat,
                          function(x){
                              tsb.info <- x %>% 
                                  dplyr::select(Tumor_Sample_Barcode,Tumor_ID) %>%
                                  dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
                              if(nrow(tsb.info) < 2){
                                  stop("Errors: each patient should have at least two tumor samples.")
                              }
                              return(tsb.info)
                          })

    
    
    ## read ccf files
    if (!is.null(ccfFile)) {
        ccf <- suppressWarnings(data.table::fread(
            ccfFile,
            quote = "",
            header = TRUE,
            fill = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
        ))
        ## check ccf data
        ccf <- validCCF(ccf)
        ## merge ccf to maf
        mafData <- readCCF(mafData, ccf, ccf.conf.level, sample.info, adjusted.VAF)
    }
    

        
        
    ## generate classMaf
    maf <- classMaf(
        data = data.table::setDT(mafData),
        sample.info = sample.info,
        mutNonSilent = mutNonSilent,
        ref.build = refBuild
    )
    
    # ## for parameter vafColumn="VAF", select particular VAF column
    # colnames(maf@data)[colnames(maf@data) == vafColumn] <- "VAF"
    
    return(maf)
}

##--- classMaf class
classMaf <- setClass(
    Class = "classMaf",
    slots = c(
        data = 'data.table',
        sample.info = 'list',
        mutNonSilent = 'character',
        ref.build = 'character'
        
    )
)



