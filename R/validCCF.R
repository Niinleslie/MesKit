validCCF <- function(ccf){
   
   ccf.standardCol <- c("Patient_ID", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "CCF")
   if(!all(ccf.standardCol %in% colnames(ccf))){
      missing_fileds_ccf <- ccf.standardCol[!ccf.standardCol %in% colnames(ccf)]
      info <- paste(missing_fileds_ccf, collapse = ",")
      stop(paste0("Missing fields from CCF data :",info) )
   }
   
   ccf$Chromosome <- as.character(ccf$Chromosome)
   ccf$CCF <- as.numeric(ccf$CCF)
   if("CCF_Std" %in% colnames(ccf)){
      ccf$CCF_Std <- as.numeric(ccf$CCF_Std)
   }

   return(ccf)
}