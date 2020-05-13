plotMutSigProfile_branch <- function(spectrum_df, mode){
   
   # products <- plyr::rbind.fill(products)
   # 
   # spectrum_df <- products %>%
   #     dplyr::filter(Branch_Tumor_ID == "Private_P", Patient_ID == "Uchi2")
   if(nrow(spectrum_df) == 0){
      return(NA)
   }
   spectrum_df <- as.data.frame(spectrum_df)
   ## calculate RSS and Cosine similarity
   reconstructed_spectrum <- spectrum_df[spectrum_df$spectrum.type == "Reconstructed",]
   original_spectrum <- spectrum_df[spectrum_df$spectrum.type == "Original",]
   m1 <- reconstructed_spectrum$Mutation_Probability
   m2 <- original_spectrum$Mutation_Probability
   RSS <- round(sum((m1 - m2)^2),3)
   cosine_sim <- round(sum(m1*m2)/(sqrt(sum(m1^2))*sqrt(sum(m2^2))),3)
   
   diff_spectrum <- original_spectrum
   diff_spectrum$Mutation_Probability <- original_spectrum$Mutation_Probability - reconstructed_spectrum$Mutation_Probability
   diff_spectrum$spectrum.type <- "Difference"
   spectrum_df <- rbind(spectrum_df,diff_spectrum)
   
   if(!is.null(mode)){
      spectrum_df <- dplyr::filter(spectrum_df, spectrum.type == mode)
   }
   
   ls.mutationType <-as.character(spectrum_df$Group)
   ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
   
   ## generate Mutation Type for every column
   spectrum_df$Group <- as.character(spectrum_df$Group)
   for (mutationGroup in ls.mutationGroup) {
      spectrum_df$Group[which(grepl(mutationGroup, spectrum_df$Group))] <- mutationGroup
   }
   
   ## specific the label order of x axis
   orderlist <- c(ls.mutationType)
   spectrum_df$Type <- factor(ls.mutationType, levels = unique(orderlist) )
   spectrum_df$Group <- factor(spectrum_df$Group, levels = ls.mutationGroup)
   
   ## specific the order of y axis
   if(is.null(mode)){
      spectrum_df$spectrum.type <- factor(spectrum_df$spectrum.type,
                                       levels = c("Original","Reconstructed","Difference"))
      spectrum_df <- dplyr::arrange(spectrum_df,spectrum.type)
   }
   else{
      spectrum_df$spectrum.type <- factor(spectrum_df$spectrum.type,
                                       levels = unique(spectrum_df$spectrum.type)) 
   }
   
   if(nrow(spectrum_df) == 0){
      warning(paste0("Patient ", patientID, ": There is no enough eligible mutations can be used."))
      return(NA)
   }
   
   
   ## set color
   group.colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF",
                     "#3C5488FF", "#F39B7FFF", "#8491B4FF")
   names(group.colors) <- ls.mutationGroup
   background.colors <- rep(c("#fce7e4","#ecf8fa","#dbfff9",
                              "#e4e8f3","#fdefeb","#e5e8ef"))
   names(background.colors) <- paste(ls.mutationGroup,1,sep = "")
   all.scales <- c(group.colors,background.colors)
   spectrum_df$Group.background <- paste(spectrum_df$Group,1,sep = "")
   
   if(!is.null(mode)){
       spectrum_df <- spectrum_df[spectrum_df$spectrum.type == mode, ]
   }

   name <- paste0(unique(spectrum_df$Patient_ID),"_",spectrum_df[1,1])

   pic <- ggplot(spectrum_df, aes(x=Type, y=Mutation_Probability, group=Group, fill=Group))
   if("Original" %in% spectrum_df$spectrum.type){
       original_spectrum <- spectrum_df[spectrum_df$spectrum.type == "Original",]
       pic <- pic + geom_rect(data = original_spectrum,
                              aes(fill = Group.background),
                              xmin = -Inf,
                              xmax = Inf,
                              ymin = min(original_spectrum$Mutation_Probability),
                              ymax = max(original_spectrum$Mutation_Probability))
   }
   if("Reconstructed" %in% spectrum_df$spectrum.type){
       reconstructed_spectrum <- spectrum_df[spectrum_df$spectrum.type == "Reconstructed",]
       pic <- pic + geom_rect(data = reconstructed_spectrum,
                              aes(fill = Group.background),
                              xmin = -Inf,
                              xmax = Inf,
                              ymin = min(reconstructed_spectrum$Mutation_Probability),
                              ymax = max(reconstructed_spectrum$Mutation_Probability))
   }
   if("Difference" %in% spectrum_df$spectrum.type){
       sig.Difference <- spectrum_df[spectrum_df$spectrum.type == "Difference",]
       pic <- pic + geom_rect(data = sig.Difference,
                              aes(fill = Group.background),
                              xmin = -Inf,
                              xmax = Inf,
                              ymin = min(sig.Difference$Mutation_Probability),
                              ymax = max(sig.Difference$Mutation_Probability))
   }
   pic <- pic + 
       geom_bar(stat="identity") +
       # geom_hline(yintercept = 0)+
       theme(
           legend.position='none', 
           strip.background = element_rect(colour = "black"),
           strip.text.x=element_text(size=14),
           strip.text.y=element_text(size=14),
           panel.spacing.y = unit(0.5, "lines"),
           panel.spacing.x = unit(0, "lines"),
           plot.title = element_text(size = 13, hjust = 0,vjust = 0),
           axis.text.x= element_text(angle = 90,vjust = 0.2,size = 5,colour = "black"),
           axis.ticks.x=element_line(color = "black"),
           axis.ticks.length.y = unit(0.2, "cm"),
           axis.text.y=element_text(size=6, color = "black"))+
       facet_grid(spectrum.type ~ Group,scales =  "free") + 
       ## color setting
       scale_fill_manual(values= all.scales) +
       scale_y_continuous(expand = c(0,0)) + 
       xlab("Mutational type") + 
       ylab("Mutation probability") + 
       ggtitle(paste(name,"\n",
                     "RSS = ", RSS, "; Cosine similarity = ", cosine_sim,"\n",
                    unique(spectrum_df$Signature),sep = ""))
   # message(paste0(patientID," mutational signature plot generation done!"))
   return(pic)
}