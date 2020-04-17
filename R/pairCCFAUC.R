ccfAUC.sample <- function(ccf){
   ccf.sort <- data.frame(CCF = as.vector(sort(ccf)), prop = c(1:length(ccf))/length(ccf))
   area <- suppressWarnings(stats::integrate(approxfun(ccf.sort$CCF,ccf.sort$prop),
                                             min(ccf.sort$CCF),
                                             max(ccf.sort$CCF),
                                             subdivisions = length(ccf),
                                             stop.on.error = F)$value) 
   return(area)
}

sortCCF <- function(df, withinType){
   
   if(withinType){
      CCF.sort <- df %>%
         dplyr::arrange(Type_Average_CCF) %>%
         dplyr::mutate(prop = c(1:nrow(.)/nrow(.)))
   }else{
      CCF.sort <- df %>%
         dplyr::arrange(CCF) %>%
         dplyr::mutate(prop = c(1:nrow(.)/nrow(.)))  
   }
}

plotDensity <- function(df,withinType){
   if(withinType){
      CCF.sort <- df %>%
         dplyr::select(Patient_ID, Tumor_Type, Type_Average_CCF) %>%
         dplyr::group_by(Tumor_Type) %>%
         dplyr::group_map(~sortCCF(.x, withinType), keep = TRUE) %>%
         do.call("rbind",.)%>%
         as.data.frame()   
   }else{
      CCF.sort <- df %>%
         dplyr::select(Patient_ID, Tumor_Sample_Barcode, CCF) %>%
         dplyr::group_by(Tumor_Sample_Barcode) %>%
         dplyr::group_map(~sortCCF(.x, withinType), keep = TRUE) %>%
         do.call("rbind",.)%>%
         as.data.frame()  
   }
   
   
   if(withinType){
      p <- ggplot2::ggplot(CCF.sort, 
                           aes(x=Type_Average_CCF, y=prop, group=Tumor_Type, color=Tumor_Type))
   }else{
      p <- ggplot2::ggplot(CCF.sort, 
                           aes(x=CCF, y=prop, group=Tumor_Sample_Barcode, color=Tumor_Sample_Barcode))
   }
   p <- p + 
      #geom_smooth(na.rm = TRUE, se = FALSE, size = 1.2, formula = y ~ s(x, bs = "cs"), method = "gam") +
      theme_bw() + 
      geom_line(size=1.2) +
      xlim(0,1) + ylim(0,1) +     
      coord_fixed() +
      theme(
         #legend.position='none', 
         legend.title = element_blank(),
         title=element_text(size=13), 
         panel.grid=element_blank(), 
         panel.border=element_blank(), 
         axis.line=element_line(size=0.7, colour = "black"),
         axis.text = element_text(size=11, colour = "black"),
         legend.text = element_text(size=11, colour = "black"),
         panel.grid.major = element_line(linetype = 2, color = "grey")
      )+
      labs(x = "CCF", y = "Proportion", 
           title = paste0("AUC plot of CCF : ", CCF.sort$Patient_ID))
   
   return(p)
   
   
}
