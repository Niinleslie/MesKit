fst.hudson.pair <- function(maf.pair) {
   
   vafCol <- which(grepl("vaf", colnames(maf.pair)))
   
   maf.pair <- mutate(
      maf.pair, 
      covariance = (vaf1-vaf2)^2-(vaf1*(1-vaf1))/(depth1-1)-(vaf2*(1-vaf2))/(depth2-1), 
      sd = vaf1*(1-vaf2)+vaf2*(1-vaf1)
   )
   
   Fst.h = mean(maf.pair$covariance)/mean(maf.pair$sd)
   return(Fst.h)
}

fst.hudson.patient <- function(df, min.vaf, plot = TRUE, use.circle = TRUE, title = NULL, 
                               withinType = FALSE, number.cex = 8, number.col = "#C77960") {
   
   df <- as.data.frame(df)
   patient <- unique(df$Patient_ID)
   ## filter by min.vaf
   if(withinType){
      
      type <- unique(df$Tumor_Type)
      
      if(length(unique(df$Tumor_Sample_Barcode))  < 2 ){
         message(paste0("Warnings: Only one sample was found of ", type,
                        " in ", patient, ". It you want to compare CCF between regions, withinType should be set as FALSE"))
         return(NA)
      }
      
      df <- df %>%
         dplyr::select(-Tumor_Type)
      # ## filter by min.vaf
      # idx <- (df %>% 
      #             dplyr::group_by(Mut_ID) %>% 
      #             dplyr::summarise(max = max(VAF_adj)) %>% 
      #             dplyr::filter(max > min.vaf) %>% 
      #             as.data.frame())$Mut_ID  
      # 
      # df <- df[df$Mut_ID %in% idx,] %>% 
      #     dplyr::select(-Tumor_Type)
      # 
      # if(length(unique(df$Tumor_Sample_Barcode))  < 2 ){
      #     message(paste0("Warnings: Only one sample was found of ", type,
      #                    " in ", patient, " after filter by min.vaf."))
      #     return(NA)
      # }
      
   }
   else{
      if(length(unique(df$Tumor_Sample_Barcode))  < 2 ){
         message(paste0("Warnings: Only one sample was found in ", patient, "."))
         return(NA)
      } 
   }
   
   
   ## pairwise heterogeneity
   patientID <- as.character(unique(df$Patient_ID))
   samples <- as.character(unique(df$Tumor_Sample_Barcode))
   pairs <- combn(length(samples), 2, simplify = FALSE)
   
   Fst.dist <- diag(1, nrow = length(samples), ncol = length(samples))
   rownames(Fst.dist) <- samples
   colnames(Fst.dist) <- samples
   
   Fst.df <- list()
   
   for (pair in pairs){
      maf.pair <- subset(df, Tumor_Sample_Barcode %in% c(samples[pair[1]],samples[pair[2]])) %>%
         tidyr::pivot_wider(
            names_from = Tumor_Sample_Barcode,       
            values_from = c(VAF_adj, totalDepth),
            values_fill = c(VAF_adj = 0, totalDepth = 0)
         ) %>%
         dplyr::ungroup() %>%
         dplyr::select(-Patient_ID)
      colnames(maf.pair) <- c("Mut_ID", "vaf1", "vaf2", "depth1", "depth2")
      
      name <- paste(samples[pair[1]],samples[pair[2]],sep = "_")
      # print(maf.pair)
      Fst.dist[pair[1],pair[2]] <- Fst.dist[pair[2],pair[1]] <- fst.hudson.pair(maf.pair)
      sub <- data.frame(patientID, name, fst.hudson.pair(maf.pair))
      Fst.df <- rbind(Fst.df,sub)
      
   }
   colnames(Fst.df) <- c("Patient_ID","Pair", "Fst")
   Fst.avg <- mean(Fst.dist)
   if(is.null(title)){
      if(withinType){
         title <- paste0("Fst of ",type," in ", patientID, ": ",round(Fst.avg,2))
      }
      else{
         title <- paste0("Fst of patient ", patientID, ": ",round(Fst.avg,2))
      }
   }
   
   return(list(
      Fst.df = Fst.df,
      Fst.avg = Fst.avg, 
      Fst.pair = Fst.dist,
      Fst.plot = if(plot) plotCorr(
         Fst.dist, 
         use.circle, 
         number.cex = number.cex,
         number.col = number.col,
         title = if(!is.null(title)) title else{paste0("Fst of patient ", patientID, ": ",round(Fst.avg,2))}) else{NA} 
   )
   )
}


groupByType <- function(df, min.vaf, plot = TRUE, use.circle = TRUE, 
                        title = NULL, withinType = FALSE, number.cex = 8, number.col = "#C77960"){
   
   types <- unique(df$Tumor_Type)
   
   Fst.type.out <- df %>% 
      group_by(Tumor_Type) %>% 
      dplyr::group_map(
         ~fst.hudson.patient(.,
                             min.vaf = min.vaf,
                             plot = plot,
                             use.circle = use.circle,
                             title = title,
                             number.cex = number.cex, 
                             number.col = number.col,
                             withinType = withinType),keep = TRUE) %>% 
      rlang::set_names(types)
   
   Fst.type.out <- Fst.type.out[!is.na(Fst.type.out)]
   
   return(Fst.type.out)
   
}