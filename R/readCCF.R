##--- combine CCF into maf object
readCCF <- function(mafData, ccf, ccf.conf.level, sample.info, adjusted.VAF,  min.average.adj.vaf) {
   ccf$Chromosome <- as.character(ccf$Chromosome)
   
   mafData_merge_ccf <-  
      dplyr::left_join(mafData, ccf, 
                       by = c("Patient_ID",
                              "Tumor_Sample_Barcode",
                               "Chromosome",
                               "Start_Position")) %>% 
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
   
   ## calculation Tumor_Average_CCF and Tumor_Average_VAF
   mafData_merge_ccf <- mafData_merge_ccf %>% 
      dplyr::group_by(Mut_ID_Type) %>% 
      dplyr::mutate(Tumor_Average_CCF = round(sum(CCF * Total_allele_depth)/sum(Total_allele_depth),3))
   
   if(!adjusted.VAF){
      mafData_merge_ccf <- mafData_merge_ccf %>%
         dplyr::mutate(VAF_adj = CCF/2) %>% 
         dplyr::mutate(Tumor_Average_VAF_adj = round(sum(VAF_adj * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
         dplyr::filter(Tumor_Average_VAF_adj >  min.average.adj.vaf)
   }
   
   mafData_merge_ccf <- dplyr::ungroup(mafData_merge_ccf)
   
   
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
         dplyr::mutate(Clonal_Status =
                          dplyr::case_when(CCF_CI_High >= 1 ~ "Clonal",
                                           CCF_CI_High < 1 ~ "Subclonal"))
      
      
      ## classify clonal status by tumor type
      ## condition1:if any region CCFm < 0.5
      c1 <- mafData_merge_ccf %>%
         dplyr::group_by(Mut_ID_Type) %>%
         dplyr::mutate(condition1 = dplyr::if_else(
            any(CCF< 0.5)  |length(CCF) == 1,
            "yes",
            "no"
         )) %>% 
         ## ungroup Mut_ID_Type
         dplyr::ungroup() %>% 
         dplyr::mutate(Clonal_Status = dplyr::if_else(
            Clonal_Status ==  "Subclonal" & condition1 == "yes",
            "Subclonal",
            "Clonal"
         ))%>%
         dplyr::select(-CCF_CI_High)

   }
   
   
   return(mafData_merge_ccf)
}

