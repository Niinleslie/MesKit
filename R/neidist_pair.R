nei_dist <- function(df, withinType = FALSE) {
   
   patient <- unique(df$Patient_ID)
   
   if(withinType){
      type <- unique(df$Tumor_Type)
      if(length(unique(df$Tumor_Sample_Barcode)) < 2){
         message(paste0("Warnings: Only one sample was found of ", type,
                        " in ", patient, ". It you want to compare CCF between regions, withinType should be set as FALSE\n"))
         return(NA)
      }
      
   }
   
   if(length(unique(df$Tumor_Sample_Barcode)) < 2){
      message(paste0("Warnings: Only one sample was found of ",patient,"."))
      return(NA)
   }
   
   df <- tidyr::pivot_wider(df,
                            names_from = Tumor_Sample_Barcode,
                            values_from = CCF,
                            values_fill = list(CCF = 0)
   ) %>%
      dplyr::select(-Patient_ID, -Mut_ID, -Tumor_Type)
   
   dist <- diag(0, nrow = ncol(df), ncol = ncol(df))
   
   dist.df <- data.frame() 
   
   for (i in 1:(ncol(df) - 1)) {
      s1 <- colnames(df)[i]
      for (j in (i + 1):ncol(df)) {
         s2 <- colnames(df)[j]
         x <- as.vector(df[, i])
         y <- as.vector(df[, j])
         
         x_ <- sum(x ^ 2 + (1 - x) ^ 2)
         y_ <- sum(y ^ 2 + (1 - y) ^ 2)
         
         xy <- sum(x * y + (1 - x) * (1 - y))
         
         
         name <- paste(s1,s2,sep = "_")
         
         sub <- data.frame(patient, name, -log(xy / sqrt(x_ * y_)))
         dist.df <-  rbind(dist.df, sub)
         
         dist[i, j] <- dist[j, i] <- -log(xy / sqrt(x_ * y_))
      }
   }
   rownames(dist) <- colnames(df)
   colnames(dist) <- colnames(df)
   return(list(dist.mat = dist, dist.df = dist.df))
   
}

groupByTypeND <- function(df){
   
   types <- unique(df$Tumor_Type)
   
   Nei.dist  <- df %>% 
      dplyr::group_by(Tumor_Type) %>% 
      dplyr::group_map(~nei_dist(., withinType = TRUE), keep = TRUE) %>% 
      rlang::set_names(types)
   
   Nei.dist <- Nei.dist[!is.na(Nei.dist)]
   
}
