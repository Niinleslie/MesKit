#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile tab delimited MAF file (plain text or *.gz compressed). Required.
#' @param ccfFile CCF file of SNVs. Default NULL.
#' @param min.vaf the minimum VAF for filtering variants. Default: 0.02.
#' @param max.vaf the maximum VAF for filtering variants. Default: 1.
#' @param min.ref.depth the minimum reference allele depth for filtering variants. Default: 8.
#' @param max.ref.depth the maximum reference allele depth for filtering variants. Default: 8.
#' @param mutType select proper variant classification you need. Default "All".Option: "nonSilent".
#' @param mutNonSilent variant classifications which are considered as non-silent. Default NULL.
#' @param chrSilent Select chromosomes needed to be dismissed. Default NULL.
#' @param use.indel logical value. whether to use INDELs besides somatic SNVs. Default FALSE.
#' @param ccf.conf.level the confidence level of CCF to identify clonal or subclonal. Only works when "CCF_std" or "CCF_CI_high" is provided in ccfFile. Default: 0.95
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
    adjusted.VAF = FALSE,
    ## filter selection
    min.vaf = 0.02,
    max.vaf = 1,
    min.average.vaf = 0,
    min.average.adj.vaf = 0,
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
                     "VAF", "Tumor_Sample_Barcode","Patient_ID","Tumor_Type")
    
    if(!all(maf.standardCol %in% colnames(mafData))){
        stop(paste("MAF file should contain Hugo_Symbol, Chromosome, Start_Position, End_Position, ",
                    "Variant_Classification, Variant_Type, Reference_Allele, ", 
                    "Tumor_Seq_Allele2, Ref_allele_depth, Alt_allele_depth, VAF, Tumor_Sample_Barcode, Patient_ID, Tumor_Type")
        )
    }

    # pre-process with gene symbols
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
                "Tumor_Type",
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2"
            ),
            sep = ":",
            remove = FALSE
        ) %>% 
        dplyr::group_by(Patient_ID,Mut_ID_Type) %>% 
        dplyr::mutate(Total_allele_depth = Ref_allele_depth + Alt_allele_depth) %>% 
        dplyr::mutate(Type_Average_VAF = round(sum(VAF * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
        dplyr::filter(Type_Average_VAF > min.average.vaf)
    
    if(adjusted.VAF){
        mafData <- mafData %>% 
            dplyr::mutate(VAF_adj = VAF) %>% 
            dplyr::group_by(Patient_ID,Mut_ID_Type) %>% 
            dplyr::mutate(Type_Average_VAF_adj = round(sum(VAF_adj * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
            dplyr::filter(Type_Average_VAF_adj > min.average.adj.vaf)
    } 
    
    mafData <- mafData %>% 
        dplyr::ungroup() %>% 
        as.data.frame()
    
  
    # Add sampleinfo
    patients.dat <- split(mafData, mafData$Patient_ID)
    sample.info <- lapply(patients.dat,
                          function(x){
                              tsb.info <- x %>% 
                                  dplyr::select(Tumor_Sample_Barcode,Tumor_Type) %>%
                                  dplyr::distinct(Tumor_Sample_Barcode, .keep_all = TRUE)
                              if(nrow(tsb.info) < 2){
                                  stop("Errors: each patient should have least two tumor samples.")
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
            stop("CCF file should contain Patient_ID,Tumor_Sample_Barcode,Chromosome,Start_Position and CCF")
        }
        

        mafData <- uniteCCF(mafData, ccf, ccf.conf.level, sample.info, adjusted.VAF, min.average.adj.vaf)
        
            #getMutStatus() %>%
            #dplyr::mutate(VAF_adj = CCF/2) ## calculate adjusted VAF based on CCF
    }
    
    mafData <- mafData %>% 
        dplyr::select(-Total_allele_depth, -Mut_ID_Type)
        ## split gene
        # tidyr::separate_rows(Hugo_Symbol,sep = ",|;") %>% 
        # dplyr::filter(Hugo_Symbol != "NONE") %>% 
        # tidyr::unite("Mut_ID",
        #              c("Patient_ID",
        #               "Chromosome",
        #               "Start_Position",
        #               "Reference_Allele",
        #               "Tumor_Seq_Allele2",
        #                "Hugo_Symbol"
        #     ),
        #     sep = ":",
        #     remove = FALSE
        # ) %>% 
        # dplyr::distinct()
        
        
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


##--- combine CCF into maf object
uniteCCF <- function(mafData, ccf, ccf.conf.level, sample.info, adjusted.VAF,  min.average.adj.vaf) {
    mafData <- tidyr::unite(
        mafData,
        "Mut_ID_CCF",
        c(
            "Patient_ID",
            "Tumor_Sample_Barcode",
            "Chromosome",
            "Start_Position"
        ),
        sep = ":",
        remove = FALSE
    )
    ccf <- ccf %>%
        tidyr::unite(
            "Mut_ID_CCF",
            c(
                "Patient_ID",
                "Tumor_Sample_Barcode",
                "Chromosome",
                "Start_Position"
            ),
            sep = ":",
            remove = TRUE
        ) 
    
    mafData_merge_ccf <-  
        merge(mafData, ccf, by = "Mut_ID_CCF", all.x = TRUE) %>% 
        dplyr::mutate(
            CCF = dplyr::if_else(
                VAF == 0,
                0,
                CCF
            )
        )%>%
        dplyr::mutate(
            CCF = dplyr::if_else(
                CCF > 1,
                1,
                CCF
            ) 
        )
    if (!"CCF_CI_High" %in% colnames(ccf) & !"CCF_Std" %in% colnames(ccf)){
        mafData_merge_ccf <-
            mafData_merge_ccf %>%
            dplyr::select(-Mut_ID_CCF)
    }else{
        if("CCF_CI_High" %in% colnames(ccf)){
            mafData_merge_ccf <- mafData_merge_ccf
        }else if("CCF_Std" %in% colnames(ccf)){
            mafData_merge_ccf <- mafData_merge_ccf %>%
                dplyr::mutate(
                    CCF_CI_High = CCF + qnorm((1 - ccf.conf.level) / 2, lower.tail = FALSE) * CCF_Std
                    ) 
        }

        mafData_merge_ccf <-  mafData_merge_ccf %>%
            dplyr::mutate(Clonal_Status =
                    dplyr::case_when(CCF_CI_High >= 1 ~ "Clonal",
                                    CCF_CI_High < 1 ~ "Subclonal"))
        
        ## classify clonal status by tumor type
        mafData_merge_ccf <- 
            tidyr::unite(
                mafData_merge_ccf,
                "Mut_ID_Sample",
                c(
                    "Patient_ID",
                    "Tumor_Sample_Barcode",
                    "Chromosome",
                    "Start_Position",
                    "Reference_Allele",
                    "Tumor_Seq_Allele2"
                ),
                sep = ":",
                remove = FALSE
            )
        
        ## if any region CCFm < 0.5
        t1 <- mafData_merge_ccf %>%
            # dplyr::filter(!is.na(CCF)) %>%
            dplyr::select(Mut_ID_Type,Tumor_Sample_Barcode, CCF) %>% 
            dplyr::group_by(Mut_ID_Type) %>%
            dplyr::summarise(condition2 = dplyr::if_else(
                any(CCF< 0.5)  |length(CCF) == 1,
                "yes",
                "no"
            )) %>% 
            as.data.frame()
        
        
        t2 <- mafData_merge_ccf %>%
            # dplyr::filter(!is.na(CCF)) %>%
            dplyr::select(Mut_ID_Type, Mut_ID_Sample, Clonal_Status) %>%
            merge(t1, by = "Mut_ID_Type", all = T) %>% 
            dplyr::mutate(Clonal_Status_2 = dplyr::if_else(
                Clonal_Status ==  "Subclonal" & condition2 == "yes",
                "Subclonal",
                "Clonal"
            )) %>% 
            dplyr::select(Mut_ID_Sample, Clonal_Status_2) %>%
            dplyr::rename(Clonal_Status = Clonal_Status_2)
        
        ## calculation Type_Average_CCF and Type_Average_VAF
        mafData_merge_ccf <- mafData_merge_ccf %>% 
            dplyr::group_by(Patient_ID,Mut_ID_Type) %>% 
            dplyr::mutate(Type_Average_CCF = round(sum(CCF * Total_allele_depth)/sum(Total_allele_depth),3))
        
        if(!adjusted.VAF){
            mafData_merge_ccf <- mafData_merge_ccf %>%
                dplyr::mutate(VAF_adj = CCF/2) %>% 
                dplyr::mutate(Type_Average_VAF_adj = round(sum(VAF_adj * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
                dplyr::filter(Type_Average_VAF_adj >  min.average.adj.vaf)
        }
        
        mafData_merge_ccf <- mafData_merge_ccf %>%
            dplyr::ungroup() %>% 
            data.table::as.data.table()
            
        ## merge clonal status
        mafData_merge_ccf <- mafData_merge_ccf %>%
            dplyr::select(-Clonal_Status) %>% 
            merge(t2, by = "Mut_ID_Sample")%>%
            dplyr::select(-Mut_ID_CCF, -Mut_ID_Sample,, -CCF_CI_High)
    }

    
    return(mafData_merge_ccf)
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
