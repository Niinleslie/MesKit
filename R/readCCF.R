##--- combine CCF into maf object
readCCF <- function(maf_data, ccf_data, ccf.conf.level, sample.info, adjusted.VAF,  min.average.adj.vaf) {
   
   mafData_merge_ccf <-  
      dplyr::left_join(maf_data, ccf_data, 
                       by = c("Patient_ID",
                              "Tumor_Sample_Barcode",
                               "Chromosome",
                               "Start_Position")) %>% 
       # dplyr::mutate(CCF = dplyr::if_else(is.na(CCF), 0,CCF))%>%
      dplyr::mutate(CCF = dplyr::if_else(.data$CCF > 1,1,.data$CCF))
   
   if(!adjusted.VAF){
      mafData_merge_ccf$VAF_adj <- mafData_merge_ccf$CCF/2
   }
   ## calculate tumor average ccf
   mafData_merge_ccf <- mafData_merge_ccf %>% 
       dplyr::group_by(.data$Patient_ID, .data$Tumor_ID, .data$Chromosome, .data$Start_Position, .data$Reference_Allele, .data$Tumor_Seq_Allele2) %>%
       dplyr::mutate(Tumor_Average_CCF = round(sum(.data$CCF * .data$Total_allele_depth)/sum(.data$Total_allele_depth),3)) %>% 
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
               CCF_CI_High = .data$CCF + stats::qnorm((1 - ccf.conf.level) / 2, lower.tail = FALSE) * .data$CCF_Std
            ) 
      }
      
      mafData_merge_ccf <-  mafData_merge_ccf %>%
         dplyr::mutate("Tumor_Mut_ID" = paste(.data$Patient_ID,
                                              .data$Tumor_ID,
                                              .data$Chromosome,
                                              .data$Start_Position,
                                              .data$Reference_Allele,
                                              .data$Tumor_Seq_Allele2,
                                              sep = ":"))
      
      mut_in_tumor_num <- mafData_merge_ccf %>%
         dplyr::group_by(.data$Tumor_Mut_ID) %>%
         dplyr::summarise(mut_in_tumor_num = length(.data$Tumor_Mut_ID))

      any_ccf_lower_05 <- mafData_merge_ccf %>%
         dplyr::group_by(.data$Tumor_Mut_ID) %>%
         dplyr::summarise(judge_any_ccf_lower_05 = dplyr::if_else(any(.data$CCF< 0.5),"yes","no"))
      
      any_ccf_ci_lower_1 <- mafData_merge_ccf %>%
         dplyr::group_by(.data$Tumor_Mut_ID) %>%
         dplyr::summarise(judge_any_ccf_ci_lower_1 = dplyr::if_else(any(.data$CCF_CI_High< 1),"yes","no"))
      
      # print(any_ccf_ci_lower_1[any_ccf_ci_lower_1$judge_any_ccf_ci_lower_1 == "yes",])
      # print(mafData_merge_ccf[mafData_merge_ccf$Tumor_Mut_ID == "V974:P:5:148745675:G:A",])
      
      mafData_merge_ccf <- mafData_merge_ccf %>%
         dplyr::left_join(mut_in_tumor_num, by = "Tumor_Mut_ID") %>%
         dplyr::left_join(any_ccf_lower_05, by = "Tumor_Mut_ID") %>% 
         dplyr::left_join(any_ccf_ci_lower_1, by = "Tumor_Mut_ID")
      
      ## classify clonal status within tumors
      ## condition1: if any region CCFm < 0.5
      mafData_merge_ccf <- mafData_merge_ccf %>%
         dplyr::group_by(.data$Tumor_Mut_ID) %>% 
         dplyr::mutate(
            Clonal_Status = dplyr::case_when(
               (.data$mut_in_tumor_num == 1 & .data$CCF_CI_High < 1) ~ "Subclonal",
               (.data$mut_in_tumor_num > 1 &
                   .data$judge_any_ccf_ci_lower_1 == "yes" & 
                   .data$judge_any_ccf_lower_05 == "yes") ~ "Subclonal",
               TRUE ~ "Clonal"
            )
         ) %>% 
         dplyr::ungroup() %>% 
         dplyr::select(-"CCF_CI_High",
                       -"Tumor_Mut_ID",
                       -"mut_in_tumor_num",
                       -"judge_any_ccf_lower_05",
                       -"judge_any_ccf_ci_lower_1")
 
   }
   
   return(mafData_merge_ccf)
}

