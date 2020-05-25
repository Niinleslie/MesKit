#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile Tab delimited MAF file (plain text or *.gz compressed). Required.
#' @param ccfFile CCF file of somatic mutations. Default NULL.
#' @param nonSyn.vc List of Variant classifications which are considered as non-silent. Default NULL, use Variant Classifications with "Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site","Nonsense_Mutation","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Missense_Mutation"
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
#' @exportClass classMaf_list
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
        dplyr::group_by(Patient_ID,Tumor_ID,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2) %>%
        dplyr::mutate(Total_allele_depth = Ref_allele_depth + Alt_allele_depth) %>% 
        dplyr::mutate(Tumor_Average_VAF = round(sum(VAF * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
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
            dplyr::group_by(Patient_ID,Tumor_ID,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2) %>%
            dplyr::mutate(Tumor_Average_VAF_adj = round(sum(VAF_adj * Total_allele_depth)/sum(Total_allele_depth),3))
    }
    
    maf_data <- maf_data %>% 
        dplyr::ungroup() %>% 
        dplyr::select(-Total_allele_depth) %>% 
        as.data.frame()
    
    data_list <- split(maf_data, maf_data$Patient_ID)
    maf_patient_list <- list()
    for(data in data_list){
        patient <- unique(data$Patient_ID)
        info <- data %>% 
            dplyr::select(Tumor_Sample_Barcode,Tumor_ID) %>%
            dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
        if(nrow(info) < 2){
            stop("Errors: each patient should have at least two tumor samples.")
        }
        maf <- classMaf(
            data = data.table::setDT(data),
            sample.info = as.data.frame(info),
            nonSyn.vc = nonSyn.vc,
            ref.build = refBuild
        )
        maf_patient_list[[patient]] <- maf
    }
    
    if(length(data_list) > 1){
        ## set claassMaf_list
        maf_list <-  classMaf_list(patient.list = maf_patient_list)
        return(maf_list)
    }else{
        return(maf_patient_list[[1]])
    }
}

##--- classMaf class
classMaf <- setClass(
    Class = "classMaf",
    slots = c(
        data = 'data.table',
        sample.info = 'data.frame',
        nonSyn.vc = 'character',
        ref.build = 'character'
        
    )
)

##--- classMaf_list class
classMaf_list <- setClass(
    Class = "classMaf_list",
    slots = c(
        patient.list = 'list'
    )
)



