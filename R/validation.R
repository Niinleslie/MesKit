validMaf <- function(maf_data){
   ## check required columns
   maf_standardcol <- c("Hugo_Symbol","Chromosome","Start_Position","End_Position",
                        "Variant_Classification", "Variant_Type", "Reference_Allele",
                        "Tumor_Seq_Allele2","Ref_allele_depth","Alt_allele_depth",
                        "VAF", "Tumor_Sample_Barcode","Patient_ID","Tumor_ID")
   
   if(!all(maf_standardcol %in% colnames(maf_data))){
      missing_fileds_maf <- maf_standardcol[!maf_standardcol %in% colnames(maf_data)]
      info <- paste(missing_fileds_maf, collapse = ",")
      stop(paste0("missing fileds from MAF :", info))
   }
   
   
   maf_data$Tumor_Sample_Barcode <- as.character(maf_data$Tumor_Sample_Barcode)
   maf_data$Patient_ID <- as.character(maf_data$Patient_ID)
   maf_data$Tumor_ID <- as.character(maf_data$Tumor_ID)
   
   ## remove VAF = 0 
   maf_data <- maf_data[VAF!=0]
   
   ## remove mutation in chromosome M and chromosome MT
   maf_data <- maf_data[!Chromosome %in% c("M", "MT")]
   
   ## sort HugoSymbol
   # maf_data <- preprocess_HugoSymbol(maf_data)
   
   ## Rescale vaf coloum 0-1
   if(max(maf_data$VAF, na.rm = TRUE) > 1){
      maf_data$VAF <- as.numeric(as.character(maf_data$VAF))/100
   }
   
   return(maf_data)
}


validCCF <- function(ccf_data){
    
    ccf_standardcol <- c("Patient_ID", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "CCF")
    if(!all(ccf_standardcol %in% colnames(ccf_data))){
        missing_fileds_ccf <- ccf_standardcol[!ccf_standardcol %in% colnames(ccf_data)]
        info <- paste(missing_fileds_ccf, collapse = ",")
        stop(paste0("Missing fields from CCF data : ",info))
    }
    
    ccf_data$Patient_ID <- as.character(ccf_data$Patient_ID)
    ccf_data$Tumor_Sample_Barcode <- as.character(ccf_data$Tumor_Sample_Barcode)
    ccf_data$Chromosome <- as.character(ccf_data$Chromosome)
    ccf_data$CCF <- as.numeric(ccf_data$CCF)
    if("CCF_Std" %in% colnames(ccf_data)){
        ccf_data$CCF_Std <- as.numeric(ccf_data$CCF_Std)
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
        stop(paste0("Missing fields from CCF data : ",info) )
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



