##--- combine CCF into maf object
readCCF <- function(maf_data, ccf_data, ccf.conf.level, sample.info, adjusted.VAF,  min.average.adj.vaf) {
   
   mafData_merge_ccf <-  
      dplyr::left_join(maf_data, ccf_data, 
                       by = c("Patient_ID",
                              "Tumor_Sample_Barcode",
                               "Chromosome",
                               "Start_Position")) %>% 
       # dplyr::mutate(CCF = dplyr::if_else(is.na(CCF), 0,CCF))%>%
      dplyr::mutate(CCF = dplyr::if_else(CCF > 1,1,CCF))
   
   if(!adjusted.VAF){
      mafData_merge_ccf$VAF_adj <- mafData_merge_ccf$CCF/2
   }
   
   ## calculate tumor average ccf
   mafData_merge_ccf <- mafData_merge_ccf %>% 
       dplyr::group_by(Patient_ID,Tumor_ID,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2) %>%
       dplyr::mutate(Tumor_Average_CCF = round(sum(CCF * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
       dplyr::ungroup() %>% 
       as.data.frame()
   
   # print(mafData_merge_ccf[is.na(mafData_merge_ccf$Tumor_Average_CCF),]$Chromosome)

   ## get clonal status of each mutation
   if (!"CCF_CI_High" %in% colnames(ccf_data) & !"CCF_Std" %in% colnames(ccf_data)){
      mafData_merge_ccf <- mafData_merge_ccf 
   }else{
      if("CCF_CI_High" %in% colnames(ccf_data)){
         mafData_merge_ccf <- mafData_merge_ccf
      }else if("CCF_Std" %in% colnames(ccf_data)){
         mafData_merge_ccf <- mafData_merge_ccf %>%
            # dplyr::mutate(CCF_Std = dplyr::if_else(is.na(CCF_Std), 0,CCF_Std))%>%
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

