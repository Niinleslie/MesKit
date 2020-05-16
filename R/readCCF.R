##--- combine CCF into maf object
readCCF <- function(mafData, ccf, ccf.conf.level, sample.info, adjusted.VAF,  min.average.adj.vaf) {
   mafData_merge_ccf <-  
      dplyr::left_join(mafData, ccf, 
                       by = c("Patient_ID",
                              "Tumor_Sample_Barcode",
                               "Chromosome",
                               "Start_Position")) %>% 
      dplyr::mutate(CCF = dplyr::if_else(VAF == 0, 0,CCF))%>%
      dplyr::mutate(CCF = dplyr::if_else(CCF > 1,1,CCF))
   
   if(!adjusted.VAF){
      mafData_merge_ccf$VAF_adj <- mafData_merge_ccf$CCF/2
   }
   
   ## get clonal status of each mutation
   if (!"CCF_CI_High" %in% colnames(ccf) & !"CCF_Std" %in% colnames(ccf)){
      mafData_merge_ccf <- mafData_merge_ccf 
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
         dplyr::mutate(Clonal_Status = dplyr::case_when(CCF_CI_High >= 1 ~ "Clonal",
                                                        CCF_CI_High < 1 ~ "Subclonal"))
      
      
      ## classify clonal status by tumor type
      ## condition1:if any region CCFm < 0.5
      mafData_merge_ccf <- mafData_merge_ccf %>%
         dplyr::group_by(Patient_ID,Tumor_ID,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2)%>%
         dplyr::mutate(condition1 = dplyr::if_else(
            any(CCF< 0.5)  |length(CCF) == 1,
            "yes",
            "no"
         )) %>% 
         dplyr::ungroup() %>% 
         dplyr::mutate(Clonal_Status = dplyr::if_else(
            Clonal_Status ==  "Subclonal" & condition1 == "yes",
            "Subclonal",
            "Clonal"
         ))%>%
         dplyr::select(-CCF_CI_High, -condition1)

   }
   
   
   return(mafData_merge_ccf)
}

