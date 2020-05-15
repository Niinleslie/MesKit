validMaf <- function(mafData){
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
   
   
   mafData$Tumor_Sample_Barcode <- as.character(mafData$Tumor_Sample_Barcode)
   mafData$Patient_ID <- as.character(mafData$Patient_ID)
   mafData$Tumor_ID <- as.character(mafData$Tumor_ID)
   
   mafData <- preprocess_HugoSymbol(mafData)
   
   ## Rescale vaf coloum 0-1
   if(max(mafData$VAF, na.rm = TRUE) > 1){
      mafData$VAF <- as.numeric(as.character(mafData$VAF))/100
   }
   
   return(mafData)
}