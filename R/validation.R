validMaf <- function(maf_data){
   ## check required columns
   maf_standardcol <- c("Hugo_Symbol","Chromosome","Start_Position","End_Position",
                        "Variant_Classification", "Variant_Type", "Reference_Allele",
                        "Tumor_Seq_Allele2","Ref_allele_depth","Alt_allele_depth",
                        "VAF", "Tumor_Sample_Barcode","Patient_ID","Tumor_ID")
   
   if(!all(maf_standardcol %in% colnames(maf_data))){
       missing_fileds_maf <- maf_standardcol[!maf_standardcol %in% colnames(maf_data)]
       info <- paste(missing_fileds_maf, collapse = ",")
       stop(paste0("Missing ", info, " from mafFile"))
   }
   
   maf_data$Chromosome <- as.character(maf_data$Chromosome)
   maf_data$Tumor_Sample_Barcode <- as.character(maf_data$Tumor_Sample_Barcode)
   maf_data$Patient_ID <- as.character(maf_data$Patient_ID)
   maf_data$Tumor_ID <- as.character(maf_data$Tumor_ID)
   
   ## remove VAF = 0 
   maf_data <- maf_data[maf_data$VAF!=0]
   
   ## remove mutation in chromosome M and chromosome MT
   maf_data <- maf_data[!maf_data$Chromosome %in% c("M", "MT")]
   
   ## sort HugoSymbol
   # maf_data <- preprocess_HugoSymbol(maf_data)
   
   ## Rescale vaf coloum 0-1
   if(max(maf_data$VAF, na.rm = TRUE) > 1){
      maf_data$VAF <- as.numeric(as.character(maf_data$VAF))/100
   }
   
   return(maf_data)
}


validCCF <- function(ccf_data, maf_data){
  
    patients_in_maf <- sort(unique(maf_data$Patient_ID))
    patients_in_ccf <- sort(unique(ccf_data$Patient_ID))  
    if(!identical(patients_in_maf, patients_in_ccf)){
      patient_setdiff <- setdiff(patients_in_maf, patients_in_ccf)
      warning("Patient: ",paste0(paste(patient_setdiff, collapse = ", "), " are not in ccf data"))
    }
    
    tsb_in_maf <- sort(unique(maf_data$Tumor_Sample_Barcode))
    tsb_in_ccf <- sort(unique(ccf_data$Tumor_Sample_Barcode)) 
    if(!identical(tsb_in_maf, tsb_in_ccf)){
      tsb_setdiff <- setdiff(tsb_in_maf, tsb_in_ccf)
      warning("Tumor sample barcodes: ",paste0(paste(tsb_setdiff, collapse = ", "), " are not in ccf data"))
    }  
    
    
    ccf_standardcol <- c("Patient_ID", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "CCF")
    if(!all(ccf_standardcol %in% colnames(ccf_data))){
        missing_fileds_ccf <- ccf_standardcol[!ccf_standardcol %in% colnames(ccf_data)]
        info <- paste(missing_fileds_ccf, collapse = ",")
        stop(paste0("Missing ", info, " from ccfFile"))
    }
    
    ccf_data$Patient_ID <- as.character(ccf_data$Patient_ID)
    ccf_data$Tumor_Sample_Barcode <- as.character(ccf_data$Tumor_Sample_Barcode)
    ccf_data$Chromosome <- as.character(ccf_data$Chromosome)
    ccf_data$CCF <- as.numeric(ccf_data$CCF)
    # if("CCF_Std" %in% colnames(ccf_data)){
    #     ccf_data$CCF_Std <- as.numeric(ccf_data$CCF_Std)
    # }
    
    if("CCF_Std" %in% colnames(ccf_data)){
      ccf_data$CCF_Std <- as.numeric(ccf_data$CCF_Std)
      ccf_data <- dplyr::select(ccf_data, "Patient_ID", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "CCF", "CCF_Std")
    }else{
      ccf_data <- dplyr::select(ccf_data, "Patient_ID", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "CCF")
    }
    
    
    return(ccf_data)
}

validSeg <- function(seg){
    seg_standardcol <- c("Patient_ID","Tumor_Sample_Barcode",
                         "Chromosome","Start_Position",
                          "End_Position")
    if(!all(seg_standardcol %in% colnames(seg))){
        missing_fileds_seg <- seg_standardcol[!seg_standardcol %in% colnames(seg)]
        info <- paste(missing_fileds_seg, collapse = ",")
        stop(paste0("Missing ", info, " from segFile"))
    }
    seg$Chromosome = gsub(pattern = 'chr', replacement = '', x = seg$Chromosome, fixed = TRUE)
    seg$Chromosome = gsub(pattern = 'X', replacement = '23', x = seg$Chromosome, fixed = TRUE)
    seg$Chromosome = gsub(pattern = 'Y', replacement = '24', x = seg$Chromosome, fixed = TRUE) 
    
    seg$Patient_ID <- as.character(seg$Patient_ID)
    seg$Tumor_Sample_Barcode <- as.character(seg$Tumor_Sample_Barcode)
    seg$Start_Position <- as.numeric(seg$Start_Position)
    seg$End_Position <- as.numeric(seg$End_Position)
    
    return(seg)
}


validClinicalData <- function(clin_data, maf_data){
  ## check Tumor_Sample_Barcode of maf data and clinical data
  clin_tb_count <- table(clin_data$Tumor_Sample_Barcode)
  if(length(which(clin_tb_count > 1)) > 0){
    rep_tb <- names(clin_tb_count)[which(clin_tb_count > 1)]
    stop(paste0("There are more than one ", paste(rep_tb, collapse = ", "), " in clinical data!"))
  }
  
  
  maf_tb <- unique(maf_data$Tumor_Sample_Barcode)
  clin_tb <- unique(clin_data$Tumor_Sample_Barcode)
  tb_setdiff <- setdiff(maf_tb, clin_tb)
  if(length(tb_setdiff) > 0){
    stop(paste0("Information about Tumor_Sample_Barcode ", paste(tb_setdiff, collapse = ", "), " cannot be found in clinical data!"))
  }
}


