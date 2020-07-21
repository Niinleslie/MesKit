#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile Tab delimited MAF file (plain text or *.gz compressed). Required.
#' @param ccfFile CCF file of somatic mutations. Default NULL.
#' @param adjusted.VAF Let VAF = VAF_adj.Default FALSE.
#' @param nonSyn.vc List of Variant classifications which are considered as non-silent. Default NULL, use Variant Classifications with "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Missense_Mutation"
#' @param ccf.conf.level The confidence level of CCF to identify clonal or subclonal. Only works when "CCF_std" or "CCF_CI_high" is provided in ccfFile. Default: 0.95
#' @param refBuild Human reference genome version. Default: 'hg19'. Optional: 'hg18' or 'hg38'.
#'
#'
#' @examples
#' maf.File <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, refBuild="hg19")
#' maf <- readMaf(mafFile=maf.File, ccfFile=ccf.File, refBuild="hg19")
#' @return an object of Maf or MafList.
#' @import methods
#' @importFrom data.table fread setkey 
#' @importFrom stats qnorm
#' @export readMaf


## read.maf main function
readMaf <- function(
    mafFile,
    ccfFile = NULL,
    adjusted.VAF = FALSE,
    nonSyn.vc = NULL,
    ccf.conf.level = 0.95,
    refBuild = "hg19") {

    refBuild <- match.arg(refBuild, choices =  c('hg18', 'hg19', 'hg38'), several.ok = FALSE)
    
    ## get non-silent muation types
    if (is.null(nonSyn.vc)) {
        nonSyn.vc <- c(
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

    maf_data <- data.table::fread(
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
    maf_data <- validMaf(maf_data)
    
    ## calculate average VAF
    maf_data <- maf_data %>% 
        dplyr::group_by(.data$Patient_ID, .data$Tumor_ID, .data$Chromosome, .data$Start_Position, .data$Reference_Allele, .data$Tumor_Seq_Allele2) %>%
        dplyr::mutate(Total_allele_depth = .data$Ref_allele_depth + .data$Alt_allele_depth) %>% 
        dplyr::mutate(Tumor_Average_VAF = round(sum(.data$VAF * .data$Total_allele_depth)/sum(.data$Total_allele_depth),3)) %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
    
    if(adjusted.VAF){
        maf_data$VAF_adj <- maf_data$VAF
    } 


    ## read ccf files
    if (!is.null(ccfFile)) {
        ccf_data <- suppressWarnings(data.table::fread(
            ccfFile,
            quote = "",
            header = TRUE,
            fill = TRUE,
            sep = '\t',
            stringsAsFactors = FALSE
        ))
        ## check ccf_data
        ccf_data <- validCCF(ccf_data)
        ## merge ccf_data to maf_data
        maf_data <- readCCF(maf_data, ccf_data, ccf.conf.level, sample.info, adjusted.VAF)
    }
    
    ## calculate average adjust VAF
    if("VAF_adj" %in% colnames(maf_data)){
        maf_data <- maf_data %>%
            dplyr::group_by(.data$Patient_ID, .data$Tumor_ID, .data$Chromosome, 
                            .data$Start_Position, .data$Reference_Allele,.data$Tumor_Seq_Allele2) %>%
            dplyr::mutate(Tumor_Average_VAF_adj = round(sum(.data$VAF_adj * .data$Total_allele_depth)/sum(.data$Total_allele_depth),3))
    }
    
    maf_data <- maf_data %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-"Total_allele_depth") %>% 
        as.data.frame()
    
    data_list <- split(maf_data, maf_data$Patient_ID)
    maf_patient_list <- list()
    for(data in data_list){
        patient <- unique(data$Patient_ID)
        sample.info <- data %>% 
            dplyr::select("Tumor_Sample_Barcode","Tumor_ID") %>%
            dplyr::distinct(.data$Tumor_Sample_Barcode, .keep_all = TRUE)
        if(nrow(sample.info) < 2){
            stop("Error: each patient should have at least two tumor samples.")
        }
        ## set Maf
        maf <- Maf(
            data = data.table::setDT(data),
            sample.info = as.data.frame(sample.info),
            nonSyn.vc = nonSyn.vc,
            ref.build = refBuild
        )
        maf_patient_list[[patient]] <- maf
    }
    
    if(length(data_list) > 1){
        ## set MafList
        maf_list <-  MafList(maf_patient_list)
        return(maf_list)
    }else{
        return(maf_patient_list[[1]])
    }
}





