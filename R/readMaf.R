#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile Tab delimited MAF file (plain text or *.gz compressed). Required.
#' @param ccfFile CCF file of somatic mutations. Default NULL.
#' @param min.vaf The minimum VAF for filtering variants. Default: 0.02.
#' @param max.vaf The maximum VAF for filtering variants. Default: 1.
#' @param min.ref.depth The minimum reference allele depth for filtering variants. Default: 4.
#' @param max.ref.depth The maximum reference allele depth for filtering variants. Default: 4.
#' @param mutType select Proper variant classification you need. Default "All". Option: "nonSilent".
#' @param mutNonSilent Variant classifications which are considered as non-silent. Default NULL.
#' @param chrSilent Select chromosomes needed to be dismissed. Default NULL.
#' @param use.indel Logical value. Whether to use INDELs besides somatic SNVs. Default FALSE.
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
    min.vaf = 0.02,
    max.vaf = 1,
    min.average.vaf = 0,
    min.average.adj.vaf = 0,
    min.ccf = 0,
    min.ref.depth = 4,
    min.alt.depth = 4,
    mutType = "All",
    mutNonSilent = NULL,
    chrSilent = NULL,
    use.indel = FALSE,
    ccf.conf.level = 0.95,
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
            data.table = TRUE,
            fill = TRUE,
            sep = '\t',
            skip = "Hugo_Symbol",
            stringsAsFactors = FALSE
        )
    
    ## check required columns
    maf.standardCol <- c("Hugo_Symbol","Chromosome","Start_Position","End_Position",
                     "Variant_Classification", "Variant_Type", "Reference_Allele",
                     "Tumor_Seq_Allele2","Ref_allele_depth","Alt_allele_depth",
                     "VAF", "Tumor_Sample_Barcode","Patient_ID","Tumor_ID")
    
    if(!all(maf.standardCol %in% colnames(mafData))){
        missing_fileds_maf <- maf.standardCol[!maf.standardCol %in% colnames(mafData)]
        info <- paste(missing_fileds_maf, collapse = ",")
        stop(paste0("missing fileds from MAF :", info))
    }

    # pre-process with gene symbols
    mafData$Tumor_Sample_Barcode <- as.character(mafData$Tumor_Sample_Barcode)
    mafData$Patient_ID <- as.character(mafData$Patient_ID)
    mafData$Tumor_ID <- as.character(mafData$Tumor_ID)

    mafData <- preprocess_HugoSymbol(mafData)

    ## Rescale vaf coloum 0-1
    if(max(mafData$VAF, na.rm = TRUE) > 1){
        mafData$VAF <- as.numeric(as.character(mafData$VAFdat))/100
    } 


    ## VAF and allele depth filter
    mafData$Ref_allele_depth[is.na(mafData$Ref_allele_depth)] <- 0
    mafData$Alt_allele_depth[is.na(mafData$Alt_allele_depth)] <- 0

    mafData <- mafData %>%
        dplyr::filter(VAF > min.vaf,
                      VAF < max.vaf,
                      Ref_allele_depth > min.ref.depth,
                      Alt_allele_depth > min.alt.depth)

    
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
    
    ## calculate type average vaf
    mafData <- mafData %>% 
        tidyr::unite(
            "Mut_ID_Type",
            c(
                "Patient_ID",
                "Tumor_ID",
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2"
            ),
            sep = ":",
            remove = FALSE
        ) %>% 
        dplyr::group_by(Mut_ID_Type) %>% 
        dplyr::mutate(Total_allele_depth = Ref_allele_depth + Alt_allele_depth) %>% 
        dplyr::mutate(Tumor_Average_VAF = round(sum(VAF * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
        dplyr::filter(Tumor_Average_VAF > min.average.vaf)
    
    if(adjusted.VAF){
        mafData <- mafData %>% 
            dplyr::mutate(VAF_adj = VAF) %>% 
            dplyr::group_by(Patient_ID,Mut_ID_Type) %>% 
            dplyr::mutate(Tumor_Average_VAF_adj = round(sum(VAF_adj * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
            dplyr::filter(Tumor_Average_VAF_adj > min.average.adj.vaf)
    } 
    
    mafData <- mafData %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
    
  
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
        
        ccf.standardCol <- c("Patient_ID", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "CCF")
        if(!all(ccf.standardCol %in% colnames(ccf))){
            missing_fileds_ccf <- ccf.standardCol[!ccf.standardCol %in% colnames(ccf)]
            info <- paste(missing_fileds_ccf, collapse = ",")
            stop(paste0("Missing fields from CCF data :",info) )
        }
        

        mafData <- readCCF(mafData, ccf, ccf.conf.level, sample.info, adjusted.VAF, min.average.adj.vaf)
        
            #getMutStatus() %>%
            #dplyr::mutate(VAF_adj = CCF/2) ## calculate adjusted VAF based on CCF
    }
    
    mafData <- mafData %>% 
        dplyr::select(-Total_allele_depth, -Mut_ID_Type)

        
        
    ## generate classMaf
    maf <- classMaf(
        data = data.table::setDT(mafData),
        #patientID = patientID,
        ref.build = refBuild,
        sample.info = sample.info
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
        ref.build = 'character'
        
    )
)



