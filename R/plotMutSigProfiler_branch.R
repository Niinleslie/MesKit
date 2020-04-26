plotMutSigProfiler_branch <- function(sig.product, mode){
   
   # products <- plyr::rbind.fill(products)
   # 
   # sig.product <- products %>%
   #     dplyr::filter(Branch_Tumor_Type == "Private_P", Patient_ID == "Uchi2")
   if(nrow(sig.product) == 0){
      return(NA)
   }
   sig.product <- as.data.frame(sig.product)
   ## calculate RSS and Cosine similarity
   sig.Reconstructed <- sig.product[sig.product$value.type == "Reconstructed",]
   sig.Original <- sig.product[sig.product$value.type == "Original",]
   m1 <- sig.Reconstructed$Mutation_Probability
   m2 <- sig.Original$Mutation_Probability
   RSS <- round(sum((m1 - m2)^2),3)
   cosine_sim <- round(sum(m1*m2)/(sqrt(sum(m1^2))*sqrt(sum(m2^2))),3)
   
   sig.diff <- sig.Original
   sig.diff$Mutation_Probability <- sig.Original$Mutation_Probability - sig.Reconstructed$Mutation_Probability
   sig.diff$value.type <- "Difference"
   sig.product <- rbind(sig.product,sig.diff)
   
   if(!is.null(mode)){
      sig.product <- dplyr::filter(sig.product, value.type == mode)
   }
   
   ls.mutationType <-as.character(sig.product$Group)
   ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
   
   ## generate Mutation Type for every column
   sig.product$Group <- as.character(sig.product$Group)
   for (mutationGroup in ls.mutationGroup) {
      sig.product$Group[which(grepl(mutationGroup, sig.product$Group))] <- mutationGroup
   }
   
   ## specific the label order of x axis
   orderlist <- c(ls.mutationType)
   sig.product$Type <- factor(ls.mutationType, levels = unique(orderlist) )
   sig.product$Group <- factor(sig.product$Group, levels = ls.mutationGroup)
   
   ## specific the order of y axis
   if(is.null(mode)){
      sig.product$value.type <- factor(sig.product$value.type,
                                       levels = c("Original","Reconstructed","Difference"))
      sig.product <- dplyr::arrange(sig.product,value.type)
   }
   else{
      sig.product$value.type <- factor(sig.product$value.type,
                                       levels = unique(sig.product$value.type)) 
   }
   
   if(nrow(sig.product) == 0){
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
   sig.product$Group.background <- paste(sig.product$Group,1,sep = "")
   
   if(!is.null(mode)){
       sig.product <- sig.product[sig.product$value.type == mode, ]
   }

   name <- paste0(unique(sig.product$Patient_ID),"_",sig.product[1,1])

   pic <- ggplot(sig.product, aes(x=Type, y=Mutation_Probability, group=Group, fill=Group))
   if("Original" %in% sig.product$value.type){
       sig.Original <- sig.product[sig.product$value.type == "Original",]
       pic <- pic + geom_rect(data = sig.Original,
                              aes(fill = Group.background),
                              xmin = -Inf,
                              xmax = Inf,
                              ymin = min(sig.Original$Mutation_Probability),
                              ymax = max(sig.Original$Mutation_Probability))
   }
   if("Reconstructed" %in% sig.product$value.type){
       sig.Reconstructed <- sig.product[sig.product$value.type == "Reconstructed",]
       pic <- pic + geom_rect(data = sig.Reconstructed,
                              aes(fill = Group.background),
                              xmin = -Inf,
                              xmax = Inf,
                              ymin = min(sig.Reconstructed$Mutation_Probability),
                              ymax = max(sig.Reconstructed$Mutation_Probability))
   }
   if("Difference" %in% sig.product$value.type){
       sig.Difference <- sig.product[sig.product$value.type == "Difference",]
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
           plot.title = element_text(size = 13, face = "bold", hjust = 0,vjust = 0),
           axis.text.x= element_text(angle = 90,vjust = 0.2,size = 5,colour = "black"),
           axis.ticks.x=element_line(color = "black"),
           axis.ticks.length.y = unit(0.2, "cm"),
           axis.text.y=element_text(size=6, color = "black"))+
       facet_grid(value.type ~ Group,scales =  "free") + 
       ## color setting
       scale_fill_manual(values= all.scales) +
       scale_y_continuous(expand = c(0,0)) + 
       xlab("Mutational type") + 
       ylab("Mutation probability") + 
       ggtitle(paste(name,"\n",
                     "RSS = ", RSS, "; Cosine similarity = ", cosine_sim,"\n",
                    unique(sig.product$Signature),sep = ""))
   # message(paste0(patientID," mutational signature plot generation done!"))
   return(pic)
}