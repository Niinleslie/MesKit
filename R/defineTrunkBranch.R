doMutTrunkBranch <- function(mtb_input,CT){
    
   ls.BT <- NA
   ## input data
   ls.BT <- .dataProcessBT(mtb_input,CT)
   # print(mtb_input$patientID)
   if(class(ls.BT) != "list"){
      return(NA)
   }
   df.pValue <- ls.BT$df.pValue
   tri_matrixBoxplot <- ls.BT$tri_matrixBoxplot
   
   if(CT){
       ls.mutationGroup <- c("C>A","C>G","C>T at CpG","C>T other","T>A","T>C","T>G")
   }else{
       ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
   }

   ## generate output data.frame with quantity of mutations in different categories
   output <- data.frame(matrix(0,nrow=length(ls.mutationGroup), ncol=2))
   colnames(output) <- c("Trunk", "Branch")
   output <- cbind(Group=ls.mutationGroup, output)
   output <- merge(output, df.pValue, by=c("Group"))
   output <- cbind(output, Significance=rep("-", nrow(output)))
   output$Significance <- as.character(output$Significance)
   for (mutationGroup in ls.mutationGroup) {
      output$Branch[which(output$Group == mutationGroup)] <- sum(tri_matrixBoxplot[which(
         tri_matrixBoxplot$Group == mutationGroup & tri_matrixBoxplot$BT == "Branch"),]$mut.num)
      output$Trunk[which(output$Group == mutationGroup)] <- sum(tri_matrixBoxplot[which(
         tri_matrixBoxplot$Group == mutationGroup & tri_matrixBoxplot$BT == "Trunk"),]$mut.num)
      
      # significant level
      if (!is.null(output[which(output$p.value < 0.05), ]$Significance)) {
         output[which(output$p.value < 0.05), c("Significance")] <- "*"
      }
      else if(!is.null(output[which(output$p.value < 0.01), ]$Significance)) {
         output[which(output$p.value < 0.01), c("Significance")] <- "**"
      }
   }
   
   output <- output %>% 
      dplyr::mutate(Patient_ID = mtb_input$patientID) %>% 
      dplyr::rename(P_Value = p.value) %>% 
      dplyr::select(Patient_ID, Group, Trunk, Branch, P_Value)
   return(output)
}

.dataProcessBT <- function(mtb_input,CT) {
   ## input data from mtb_input
   tri_matrix <- mtb_input$tri_matrix
   ## label the Trunk
   if (length(mtb_input$trunk_name) != 0){
      trunk_name <- mtb_input$trunk_name
   } else {
      warning(paste0("Patient ",mtb_input$patientID,": no trunk mutations are detected!"))
      return(NA)
   } 
   ## separate trunk and branch data
   tri_matrix.trunk <- tri_matrix[which(rownames(tri_matrix) == trunk_name), ]
   tri_matrix.branch <- tri_matrix[which(rownames(tri_matrix) != trunk_name), ]
   tri_matrix.branch <- colSums(tri_matrix.branch)
   tri_matrixBT <- rbind(Trunk=tri_matrix.trunk, Branch=tri_matrix.branch)
   tri_matrixBTTrans <- data.frame(Mutational_Type=colnames(tri_matrixBT), t(tri_matrixBT))
   
   if(CT){
       ls.mutationGroup <- c("C>A","C>G","C>T at CpG","C>T other","T>A","T>C","T>G")
   }else{
       ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
   }

   ## generate Mutation Type for every column
   # print(mtb_input$patientID)
   # print(tri_matrixBTTrans)
   tri_matrixBTTrans$Group <- ''
   for (mutationGroup in ls.mutationGroup) {
       # print(mutationGroup)
      tri_matrixBTTrans$Group[which(grepl(mutationGroup, tri_matrixBTTrans$Mutational_Type))] <- mutationGroup
   }
   # print(1)
   tri_matrixBSum <- tri_matrixBTTrans %>% dplyr::group_by(Group) %>% dplyr::summarise(sum = sum(Branch))
   tri_matrixTSum <- tri_matrixBTTrans %>% dplyr::group_by(Group) %>% dplyr::summarise(sum = sum(Trunk))
   
   tri_matrixBTTrans <- cbind(tri_matrixBTTrans, 
                             BranchFrac=rep(0, nrow(tri_matrixBTTrans)), 
                             TrunkFrac=rep(0, nrow(tri_matrixBTTrans)))
   for (mutationGroup in ls.mutationGroup) {
      groupBSum <- tri_matrixBSum$sum[which(tri_matrixBSum$Group == mutationGroup)]
      if (groupBSum == 0) {
         tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$Branch <- 0
         # tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$BranchFrac <- 0
      }
      groupTSum <- tri_matrixTSum$sum[which(tri_matrixTSum$Group == mutationGroup)]
      if (groupBSum == 0) {
         tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$Trunk <- 0
         # tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$TrunkFrac <- 0
      }
      
      tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$BranchFrac <-
         100*tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$Branch/groupBSum
      tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$TrunkFrac <-
         100*tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]$Trunk/groupTSum
   }
   
   tri_matrixBoxplot <- data.frame(matrix(nrow=0, ncol=5))
   colnames(tri_matrixBoxplot) <- c("GroupBT", "Group", "BT", "mut.frac", "mut.num")
   for (mutationGroup in ls.mutationGroup) {
      dat.group <- tri_matrixBTTrans[which(tri_matrixBTTrans$Group == mutationGroup), ]
      df.groupT <- data.frame(rep(paste(mutationGroup, "Trunk", sep=" "), nrow(dat.group)), 
                              rep(mutationGroup, nrow(dat.group)), 
                              rep("Trunk", nrow(dat.group)),
                              dat.group$TrunkFrac, 
                              dat.group$Trunk)
      df.groupB <- data.frame(rep(paste(mutationGroup, "Branch", sep=" "), nrow(dat.group)), 
                              rep(mutationGroup, nrow(dat.group)), 
                              rep("Branch", nrow(dat.group)), 
                              dat.group$BranchFrac, 
                              dat.group$Branch)
      colnames(df.groupB) <- c("GroupBT", "Group", "BT", "mut.frac", "mut.num")
      colnames(df.groupT) <- c("GroupBT", "Group", "BT", "mut.frac", "mut.num")
      tri_matrixBoxplot <- rbind(tri_matrixBoxplot, df.groupT, df.groupB)
   }
   
   df.pValue <- data.frame(matrix(ncol = 2, nrow = 0))
   colnames(df.pValue) <- c("Group", "p.value")
   for (mutationGroup in ls.mutationGroup) {
      branch.mut.num <- sum(tri_matrixBoxplot[which(tri_matrixBoxplot$Group == mutationGroup &  tri_matrixBoxplot$BT == "Branch"), ]$mut.num) 
      branch.mut.num2 <- sum(tri_matrixBoxplot[which(tri_matrixBoxplot$Group != mutationGroup &  tri_matrixBoxplot$BT == "Branch"), ]$mut.num)
      trunk.mut.num <-  sum(tri_matrixBoxplot[which(tri_matrixBoxplot$Group == mutationGroup & tri_matrixBoxplot$BT == "Trunk"), ]$mut.num) 
      trunk.mut.num2 <-  sum(tri_matrixBoxplot[which(tri_matrixBoxplot$Group != mutationGroup & tri_matrixBoxplot$BT == "Trunk"), ]$mut.num) 
      # pValue <- wilcox.test(branch.mut.frac,
      #                       trunk.mut.frac,
      #                       paired=TRUE,
      #                       alternative = "two.sided",
      #                       exact=FALSE)$p.value
      if(all(!is.nan(branch.mut.num)) & all(!is.nan(trunk.mut.num))){
         # message("branch.mut.frac")
         # print(branch.mut.frac)
         # message("trunk.mut.frac")
         # print(sum(trunk.mut.frac) )
         
         m <- matrix(c(branch.mut.num, trunk.mut.num, branch.mut.num2, trunk.mut.num2),ncol = 2)
         pValue <- fisher.test(m,alternative = "two.sided")$p.value
      }
      else{
         # warning(paste0("Patient ", mtb_input$patientID, ": There is no enough eligible mutations can be used."))
         pValue <- NA
      }
      row.pValue <- data.frame(mutationGroup, pValue)
      colnames(row.pValue) <- c("Group", "p.value")
      df.pValue <- rbind(df.pValue, row.pValue)
   }
   output <- list(df.pValue=df.pValue, tri_matrixBoxplot=tri_matrixBoxplot)
   return(output)
}


doPlotTrunkBranch <- function(mtb_output, pvalue = 0.05, CT){
   if(sum(mtb_output$Trunk) == 0){
       mtb_output$trunk.frac <- 0
   }
   else{
       mtb_output$trunk.frac <- mtb_output$Trunk/sum(mtb_output$Trunk)  
   }
   if(sum(mtb_output$Branch) == 0){
       mtb_output$branch.frac <-  0
   }
   else{
       mtb_output$branch.frac <- mtb_output$Branch/sum(mtb_output$Branch)
   }
   dat <- mtb_output %>% 
       tidyr::pivot_longer(cols = c("trunk.frac","branch.frac"),names_to = "BT",values_to = "fraction") %>%
       # reshape2::melt(id.vars = c("Group","P_Value","Patient_ID","Branch","Trunk"),variable.name = "BT",value.name = "fraction") %>% 
       dplyr::rowwise() %>% 
       dplyr::mutate(group.name = paste0(Patient_ID,BT)) %>% 
       as.data.frame() 
   dat$group.name <- factor(dat$group.name, levels = unique(dat$group.name))
   dat <- dplyr::arrange(dat, group.name, plyr::desc(Group))
   ## set data table for bar plot
   dat <- as.data.table(dat)
   dat$rect.xmin <- 0
   dat$rect.xmax <- 0
   dat$rect.ymin <- 0
   dat$rect.ymax <- 0
   groups <- unique(dat$group.name)
   # patient1 <- unique(dat[group.name == groups[1],]$Patient_ID)
   i <- 0
   for(g in groups){
      ## get position for x axis
      # patient2 <- unique(dat[group.name == g,]$Patient_ID)
      # if(patient1 != patient2){
      #    ## The space between the patients
      #    i <- i + 0.1
      # }
      # patient1 <- patient2
      dat[group.name == g,]$rect.xmin <- i
      dat[group.name == g,]$rect.xmax <- i + 0.1
      i <- i + 0.15
      
      ## get position for y axis
      yy <- cumsum(dat[group.name == g]$fraction)
      len.yy <- length(yy)
      dat[group.name == g,]$rect.ymin <- c(0,yy[-len.yy])
      dat[group.name == g,]$rect.ymax <- yy
      
   }
   
   ## set table for lable of patientid
   patients <- unique(dat$Patient_ID)
   patient.text.table <- data.table(Patient_ID = patients)
   patient.text.table$p.x <- 0
   for(p in patients){
      p.xmin <- min(dat[Patient_ID == p]$rect.xmin)
      p.xmax <- max(dat[Patient_ID == p]$rect.xmax)
      patient.text.table[Patient_ID == p]$p.x <- p.xmin + (p.xmax - p.xmin)/2
   }
   
   ## set table for number of mutation
   num.table <- data.table(group.name = groups)
   num.table$g.x <- 0
   num.table$num <- 0
   for(g in groups){
      g.xmin <- min(dat[group.name == g]$rect.xmin)
      g.xmax <- max(dat[group.name == g]$rect.xmax)
      num.table[group.name == g]$g.x <- g.xmin + (g.xmax-g.xmin)/2
      if(grepl("trunk",g)){
         num.table[group.name == g]$num <- sum(dat[group.name == g]$Trunk) 
      }
      else{
         num.table[group.name == g]$num <- sum(dat[group.name == g]$Branch) 
      }
   }
   
   ## set table for pvalue
   dat[,PG:=paste0(Patient_ID,Group)]
   pvalue.table <- data.table()
   for(pg in unique(dat$PG)){
      v <- unique(dat[PG == pg]$P_Value) 
      if(is.na(v)){
         next
      }
      if(v < pvalue){
         x <- max(dat[PG == pg]$rect.xmax)
         b.ymin <- dat[PG == pg]$rect.ymin[2]
         b.ymax <- dat[PG == pg]$rect.ymax[2]
         y <- b.ymin + (b.ymax - b.ymin)/2
         sub <- data.table(g,x,y,"*")
         pvalue.table <- rbind(pvalue.table,sub)
      }
      # else if(v < 0.05){
      #     x <- max(dat[PG == pg]$rect.xmax)
      #     b.ymin <- dat[PG == pg]$rect.ymin[2]
      #     b.ymax <- dat[PG == pg]$rect.ymax[2]
      #     y <- b.ymin + (b.ymax - b.ymin)/2
      #     sub <- data.table(g,x,y,"*")
      #     pvalue.table <- rbind(pvalue.table,sub)
      # } 
   }
   
   ## set segment from trunk to branches
   segment.table <- data.table()
   for(pg in unique(dat$PG)){
      fractions <- dat[PG == pg]$fraction
      if(all(fractions == 0)){
         next
      }
      x1 <- min(dat[PG == pg]$rect.xmax)
      x2 <- max(dat[PG == pg]$rect.xmin)
      y1 <- dat[PG == pg]$rect.ymax[1]
      y2 <- dat[PG == pg]$rect.ymax[2]
      sub <- data.table(x1,x2,y1,y2)
      segment.table <- rbind(segment.table,sub)
   }
   
   ## set color 
   if(CT){
       all.colors <- c("C>A" = "#E64B35FF",
                       "C>G" = "#4DBBD5FF",
                       "T>A" = "#3C5488FF",
                       "T>C" = "#F39B7FFF",
                       "T>G" = "#8491B4FF",
                       "C>T at CpG" = "#00A087FF",
                       "C>T other" = "#91D1C2FF") 
   }else{
       all.colors <- c("C>A" = "#E64B35FF",
                       "C>G" = "#4DBBD5FF",
                       "C>T" = "#00A087FF",
                       "T>A" = "#3C5488FF",
                       "T>C" = "#F39B7FFF",
                       "T>G" = "#8491B4FF") 
   }

   group.colors <- all.colors[unique(dat$Group)]
   
   ## plot
   pic <- ggplot() + 
      geom_rect(data = dat,
                aes(xmin = rect.xmin, xmax = rect.xmax,
                    ymin = rect.ymin, ymax = rect.ymax,
                    fill = Group),color = "black")+ 
      
      ## lable patientid
      geom_text(data = patient.text.table,
                aes(x = p.x, y = 1.1, label = Patient_ID),size = 6)+

      # label number of mutation
      geom_text(data = num.table,
                aes(x = g.x, y = 1.03, label = num),size = 4)+
      
      
      geom_segment(aes(y = 0 ,yend = 1,x=-Inf,xend=-Inf), size = 1.5)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.48,size = 17,vjust = -2),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 11,colour = "black"),
            axis.text.x = element_text(angle = 60,size = 11,colour = "black",hjust = 1,margin = margin(t = -10)),
            axis.text.y = element_text(size = 10,colour = "black"),
            axis.ticks.length.y = unit(.25, "cm"),
            axis.ticks.y = element_line(size = 1)) +
      scale_fill_manual(values = group.colors) + 
      xlab("") + 
      ylab("Proportion")+
      # ggtitle(unique(dat$Patient_ID)) + 
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1))+
      scale_x_continuous(breaks = unique(dat$rect.xmin+(dat$rect.xmax - dat$rect.xmin)/2) ,
                         labels = rep(c("Trunk","Branches"),length(unique(dat$Patient_ID))))
   
   ## label significant pvalue
   if(nrow(pvalue.table) > 0){
      colnames(pvalue.table) <- c("group.name","pv.x","pv.y","pv.label")
      pic <- pic + 
         geom_text(data = pvalue.table,
                   aes(x = pv.x+0.01, y = pv.y, label = pv.label),
                   size = 7)
   }
   
   ## segment between trunk and branches
   if(nrow(segment.table) > 0){
      pic <- pic + 
         geom_segment(data = segment.table,
                      aes(x = x1, xend = x2, y = y1, yend = y2))
   }
   
   return(pic)
}

