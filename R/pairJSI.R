JSI_dist <- function(df, pairByType){
   
   patientid <- unique(df$Patient_ID) 

   if(pairByType){
      if(length(unique(df$Tumor_Type))  < 2 ){
         message(paste0("Warnings: Only one tumor type was found of ",patientid, ". It you want to compare CCF between regions, pairByType should be set as FALSE"))
         return(NA)
      }
   }else{
      if(length(unique(df$Tumor_Sample_Barcode))  < 2 ){
         message(paste0("Warnings: Only one sample was found in ", patientid, "."))
         return(NA)
      } 
   }
   
   ## pairwise heterogeneity
   if(pairByType){
      types <- as.character(unique(df$Tumor_Type))
      pairs <- combn(length(types), 2, simplify = FALSE)
      dist <- diag(1, nrow = length(types), ncol = length(types))
      rownames(dist) <- types
      colnames(dist) <- types
   }else{
      samples <- as.character(unique(df$Tumor_Sample_Barcode))
      pairs <- combn(length(samples), 2, simplify = FALSE)
      dist <- diag(1, nrow = length(samples), ncol = length(samples))
      rownames(dist) <- samples
      colnames(dist) <- samples
   }
   
   PC_1.list <- c()
   PC_2.list <- c()
   SS_12.list <- c()
   JSI.df <- data.frame()
   for (pair in pairs){
      
      if(pairByType){
         name <- paste(types[pair[1]],types[pair[2]], sep = "_")
         vaf.pair <- subset(df, Tumor_Type %in% c(types[pair[1]],types[pair[2]])) %>%
            tidyr::unite("mutation_id2",
                         c("mutation_id",
                           "Tumor_Type"),
                         sep = ":",
                         remove = FALSE
            ) %>%
            dplyr::distinct(mutation_id2, .keep_all = T) %>%
            dplyr::select(mutation_id, Tumor_Type, Clonal_Status, VAF_adj) %>% 
            tidyr::pivot_wider(
               names_from = Tumor_Type,       
               values_from = c(VAF_adj, Clonal_Status),
               values_fill = c(VAF_adj = 0, Clonal_Status = 'NA')
            ) %>%
            dplyr::ungroup()
         colnames(vaf.pair) <- c("mutation_id", "vaf1", "vaf2", "status1", "status2")
         
         
      }
      else{
         name <- paste(samples[pair[1]],samples[pair[2]], sep = "_")
         
         vaf.pair <- subset(df, Tumor_Sample_Barcode %in% c(samples[pair[1]],samples[pair[2]])) %>%
            dplyr::select(mutation_id, Tumor_Sample_Barcode, Clonal_Status, VAF_adj) %>% 
            tidyr::pivot_wider(
               names_from = Tumor_Sample_Barcode,       
               values_from = c(VAF_adj, Clonal_Status),
               values_fill = c(VAF_adj = 0, Clonal_Status = 'NA')
            ) %>%
            dplyr::ungroup()
         colnames(vaf.pair) <- c("mutation_id", "vaf1", "vaf2", "status1", "status2")
      }
      
      vaf.pair <- vaf.pair %>% 
         dplyr::filter(vaf1 + vaf2 !=0) %>% 
         data.table::setDT()
      
      
      PC_1 <- nrow(vaf.pair[status1 == "Clonal" & vaf1>0 & vaf2==0])
      PC_1.list <- c(PC_1.list, PC_1)
      PC_2 <- nrow(vaf.pair[status2 == "Clonal" & vaf1==0 & vaf2>0])
      PC_2.list <- c(PC_2.list, PC_2)
      SS_12 = nrow(vaf.pair[status2=="Subclonal" & status1 == "Subclonal" & vaf1>0 & vaf2>0 ])
      SS_12.list <- c(SS_12.list, SS_12)
      
      
      dist[pair[1],pair[2]] <- dist[pair[2],pair[1]] <- SS_12/(PC_1+PC_2+SS_12)
      
      sub <- data.frame(patientid,name, SS_12/(PC_1+PC_2+SS_12))
      JSI.df <- rbind(JSI.df, sub)
   }
   
   multi <- mean(SS_12.list)/(mean(PC_1.list) + mean(PC_2.list) + mean(SS_12.list))
   JSI.multi <- data.frame(Patient_ID = patientid, JSI.multi = multi)
   
   #return(JSI.dist)
   return(list(JSI.multi = JSI.multi, JSI.pair =  dist, JSI.df = JSI.df))
}    