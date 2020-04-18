doMutTrunkBranch <- function(tree.mutSig){
   ls.BT <- NA
   ## input data from tree.mutSig
   ls.BT <- .dataProcessBT(tree.mutSig)
   # print(tree.mutSig$patientID)
   if(is.na(ls.BT)){
      return(NA)
   }
   df.pValue <- ls.BT$df.pValue
   sigsInputBoxplot <- ls.BT$sigsInputBoxplot
   ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
   
   ## generate output data.frame with quantity of mutations in different categories
   output <- data.frame(matrix(nrow=6, ncol=2))
   colnames(output) <- c("Trunk", "Branch")
   output <- cbind(Group=ls.mutationGroup, output)
   output <- merge(output, df.pValue, by=c("Group"))
   output <- cbind(output, Significance=rep("-", nrow(output)))
   for (mutationGroup in ls.mutationGroup) {
      output$Branch[which(output$Group == mutationGroup)] <- sum(sigsInputBoxplot[which(
         sigsInputBoxplot$Group == mutationGroup & sigsInputBoxplot$BT == "Branch"),]$mut.num)
      output$Trunk[which(output$Group == mutationGroup)] <- sum(sigsInputBoxplot[which(
         sigsInputBoxplot$Group == mutationGroup & sigsInputBoxplot$BT == "Trunk"),]$mut.num)
      
      # significant level
      if (!is.null(output[which(output$p.value < 0.05), ]$Significance)) {
         output[which(output$p.value < 0.05), c("Significance")] <- "*"
      }
      else if(!is.null(output[which(output$p.value < 0.01), ]$Significance)) {
         output[which(output$p.value < 0.01), c("Significance")] <- "**"
      }
   }
   
   output <- output %>% 
      dplyr::mutate(Patient_ID = tree.mutSig$patientID) %>% 
      dplyr::rename(P_Value = p.value) %>% 
      dplyr::select(Patient_ID, Group, Trunk, Branch, P_Value)
   return(output)
}

.dataProcessBT <- function(tree.mutSig) {
   ## input data from tree.mutSig
   sigsInput <- tree.mutSig$sigsInput
   ## label the Trunk
   if (length(tree.mutSig$trunkName) != 0){
      trunkName <- tree.mutSig$trunkName
   } else {
      warning(paste0("Patient ",tree.mutSig$patientID,": no trunk mutations are detected!"))
      return(NA)
   } 
   ## separate trunk and branch data
   sigsInput.trunk <- sigsInput[which(rownames(sigsInput) == trunkName), ]
   sigsInput.branch <- sigsInput[which(rownames(sigsInput) != trunkName), ]
   sigsInput.branch <- colSums(sigsInput.branch)
   sigsInputBT <- rbind(Trunk=sigsInput.trunk, Branch=sigsInput.branch)
   sigsInputBTTrans <- data.frame(Mutational_Type=colnames(sigsInputBT), t(sigsInputBT))
   ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
   
   ## generate Mutation Type for every column
   for (mutationGroup in ls.mutationGroup) {
      sigsInputBTTrans$Group[which(grepl(mutationGroup, sigsInputBTTrans$Mutational_Type))] <- mutationGroup
   }
   
   sigsInputBSum <- sigsInputBTTrans %>% dplyr::group_by(Group) %>% dplyr::summarise(sum = sum(Branch))
   sigsInputTSum <- sigsInputBTTrans %>% dplyr::group_by(Group) %>% dplyr::summarise(sum = sum(Trunk))
   
   sigsInputBTTrans <- cbind(sigsInputBTTrans, 
                             BranchFrac=rep(0, nrow(sigsInputBTTrans)), 
                             TrunkFrac=rep(0, nrow(sigsInputBTTrans)))
   for (mutationGroup in ls.mutationGroup) {
      groupBSum <- sigsInputBSum$sum[which(sigsInputBSum$Group == mutationGroup)]
      if (groupBSum == 0) {
         sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Branch <- 0
         # sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$BranchFrac <- 0
      }
      groupTSum <- sigsInputTSum$sum[which(sigsInputTSum$Group == mutationGroup)]
      if (groupBSum == 0) {
         sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Trunk <- 0
         # sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$TrunkFrac <- 0
      }
      
      sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$BranchFrac <-
         100*sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Branch/groupBSum
      sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$TrunkFrac <-
         100*sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Trunk/groupTSum
   }
   
   sigsInputBoxplot <- data.frame(matrix(nrow=0, ncol=5))
   colnames(sigsInputBoxplot) <- c("GroupBT", "Group", "BT", "mut.frac", "mut.num")
   for (mutationGroup in ls.mutationGroup) {
      dat.group <- sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]
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
      sigsInputBoxplot <- rbind(sigsInputBoxplot, df.groupT, df.groupB)
   }
   
   df.pValue <- data.frame(matrix(ncol = 2, nrow = 0))
   colnames(df.pValue) <- c("Group", "p.value")
   for (mutationGroup in ls.mutationGroup) {
      branch.mut.num <- sum(sigsInputBoxplot[which(sigsInputBoxplot$Group == mutationGroup &  sigsInputBoxplot$BT == "Branch"), ]$mut.num) 
      branch.mut.num2 <- sum(sigsInputBoxplot[which(sigsInputBoxplot$Group != mutationGroup &  sigsInputBoxplot$BT == "Branch"), ]$mut.num)
      trunk.mut.num <-  sum(sigsInputBoxplot[which(sigsInputBoxplot$Group == mutationGroup & sigsInputBoxplot$BT == "Trunk"), ]$mut.num) 
      trunk.mut.num2 <-  sum(sigsInputBoxplot[which(sigsInputBoxplot$Group != mutationGroup & sigsInputBoxplot$BT == "Trunk"), ]$mut.num) 
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
         # warning(paste0("Patient ", tree.mutSig$patientID, ": There is no enough eligible mutations can be used."))
         pValue <- NA
      }
      row.pValue <- data.frame(mutationGroup, pValue)
      colnames(row.pValue) <- c("Group", "p.value")
      df.pValue <- rbind(df.pValue, row.pValue)
   }
   output <- list(df.pValue=df.pValue, sigsInputBoxplot=sigsInputBoxplot)
   return(output)
}


doPlotTrunkBranch <- function(mutTrunkBranch.list, pvalue = 0.05){
   mutTrunkBranch.list <- lapply(mutTrunkBranch.list,function(x){
      if(sum(x$Trunk) == 0){
          x$trunk.frac <- 0
      }
       else{
           x$trunk.frac <- x$Trunk/sum(x$Trunk)  
       }
       if(sum(x$Branch) == 0){
           x$branch.frac <-  0
       }
       else{
           x$branch.frac <- x$Branch/sum(x$Branch)
       }
      x <- x %>% 
         tidyr::pivot_longer(cols = c("trunk.frac","branch.frac"),names_to = "BT",values_to = "fraction") %>%
         # reshape2::melt(id.vars = c("Group","P_Value","Patient_ID","Branch","Trunk"),variable.name = "BT",value.name = "fraction") %>% 
         dplyr::rowwise() %>% 
         dplyr::mutate(group.name = paste0(Patient_ID,BT)) %>% 
         as.data.frame() 
      return(x)
   })
   # mutTrunkBranch.list$trunk.frac <- mutTrunkBranch.list$Trunk/sum(mutTrunkBranch.list$Trunk)
   # mutTrunkBranch.list$branch.frac <- mutTrunkBranch.list$Branch/sum(mutTrunkBranch.list$Branch)
   # mutTrunkBranch.list <- mutTrunkBranch.list %>%
   #     melt(id.vars = c("Group","p.value","Significance","Patient_ID","Branch","Trunk"),variable.name = "BT",value.name = "fraction") %>% 
   #     dplyr::rowwise() %>% 
   #     dplyr::mutate(group.name = paste0(Patient_ID,BT)) %>% 
   #     as.data.frame() 
   dat <- plyr::rbind.fill(mutTrunkBranch.list)
   dat$group.name <- factor(dat$group.name, levels = unique(dat$group.name))
   
   ## set data table for bar plot
   dat <- as.data.table(dat)
   dat$rect.xmin <- 0
   dat$rect.xmax <- 0
   dat$rect.ymin <- 0
   dat$rect.ymax <- 0
   groups <- unique(dat$group.name)
   patient1 <- unique(dat[group.name == groups[1],]$Patient_ID)
   i <- 0
   for(g in groups){
      ## get position for x axis
      patient2 <- unique(dat[group.name == g,]$Patient_ID)
      if(patient1 != patient2){
         ## The space between the patients
         i <- i + 0.1
      }
      patient1 <- patient2
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
      if(any(fractions == 0)){
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
   all.colors <- c("C>A" = "#E64B35FF",
                   "C>G" = "#4DBBD5FF",
                   "C>T" = "#00A087FF",
                   "T>A" = "#3C5488FF",
                   "T>C" = "#F39B7FFF",
                   "T>G" = "#8491B4FF")
   group.colors <- all.colors[unique(dat$Group)]
   
   ## plot
   pic <- ggplot() + 
      geom_rect(data = dat,
                aes(xmin = rect.xmin, xmax = rect.xmax,
                    ymin = rect.ymin, ymax = rect.ymax,
                    fill = Group),color = "black")+ 
      
      ## lable patientid
      geom_text(data = patient.text.table,
                aes(x = p.x, y = 1.08, label = Patient_ID),size = 5)+
      
      ## label number of mutation
      geom_text(data = num.table,
                aes(x = g.x, y = 1.03, label = num),size = 4)+
      
      
      geom_segment(aes(y = 0 ,yend = 1,x=-Inf,xend=-Inf), size = 1.5)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_text(size = 11,colour = "black"),
            axis.text.x = element_text(angle = 60,size = 11,colour = "black",hjust = 1,margin = margin(t = -10)),
            axis.text.y = element_text(size = 10,colour = "black"),
            axis.ticks.length.y = unit(.25, "cm"),
            axis.ticks.y = element_line(size = 1)) +
      scale_fill_manual(values = group.colors) + 
      xlab("") + 
      ylab("Proportion")+
      scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1))+
      scale_x_continuous(breaks = unique(dat$rect.xmin+(dat$rect.xmax - dat$rect.xmin)/2) ,
                         labels = rep(c("Trunk","Branches"),length(unique(dat$Patient_ID))))
   
   ## label significant pvalue
   if(nrow(pvalue.table) > 0){
      colnames(pvalue.table) <- c("group.name","pv.x","pv.y","pv.label")
      pic <- pic + 
         geom_text(data = pvalue.table,
                   aes(x = pv.x+0.05, y = pv.y, label = pv.label),
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

# getYmax <- function(BT.dat){
#     if(nrow(BT.dat) == 0){
#         return(0)
#     }
#     BT.dat <- as.data.table(BT.dat)
#     
#     ## maximum mut.frac in Branch
#     Branch.ymax <- 0
#     Branch.mutfrac <- BT.dat[BT == "Branch"]$mut.frac
#     if(length(Branch.mutfrac)> 0){
#         Branch.outlier <- boxplot(Branch.mutfrac,range = 100,plot = FALSE)$out
#         Branch.ymax <- max(Branch.mutfrac[!Branch.mutfrac %in% Branch.outlier])  
#     }
#     
#     ## maximum mut.frac in Trunk
#     Trunk.ymax <- 0
#     Trunk.mutfrac <- BT.dat[BT == "Trunk"]$mut.frac
#     if(length(Trunk.mutfrac) > 0){
#         Trunk.outlier <- boxplot(Trunk.mutfrac,range = 100, plot = FALSE)$out
#         Trunk.ymax <- max(Trunk.mutfrac[!Trunk.mutfrac %in% Trunk.outlier])  
#     }
#     
#     return(max(Trunk.ymax,Branch.ymax) + 5)
# }
# 
# 
# doPlotTrunkBranch <- function(tree.mutSig){
#     ## input data from tree.mutSig
#     ls.BT <- .dataProcessBT(tree.mutSig)
#     if(any(is.na(ls.BT)) ){
#         return(NA)
#     }
#     df.pValue <- ls.BT$df.pValue
#     sigsInputBoxplot <- ls.BT$sigsInputBoxplot
#     
#     ## p values of mutational list
#     if(is.na(df.pValue[which(df.pValue$Group == "C>A"), ]$p.value)){
#         CApV <- grid::textGrob(expression(""), 
#                                gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75) 
#     }
#     else{
#         if (df.pValue[which(df.pValue$Group == "C>A"), ]$p.value < 0.01) {
#             CApV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "C>A"), ]$p.value, 
#                       digits = 3)), "**", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)    
#         }else if (df.pValue[which(df.pValue$Group == "C>A"), ]$p.value < 0.05) {
#             CApV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "C>A"), ]$p.value, 
#                       digits = 3)), "*", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)    
#         } else {
#             CApV <- grid::textGrob(expression(""), 
#                                    gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75) 
#         }
#     }
#     if (is.na(df.pValue[which(df.pValue$Group == "C>G"), ]$p.value)){
#         CGpV <- grid::textGrob(expression(""), 
#                                gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75) 
#     }else{
#         if (df.pValue[which(df.pValue$Group == "C>G"), ]$p.value < 0.01) {
#             CGpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "C>G"), ]$p.value, 
#                       digits = 3)), "**", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else if (df.pValue[which(df.pValue$Group == "C>G"), ]$p.value < 0.05) {
#             CGpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "C>G"), ]$p.value, 
#                       digits = 3)), "*", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else {
#             CGpV <- grid::textGrob(expression(""), 
#                                    gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         }
#     }
#     if (is.na(df.pValue[which(df.pValue$Group == "C>T"), ]$p.value)){
#         CTpV <- grid::textGrob(expression(""), 
#                                gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#     }else{
#         if (df.pValue[which(df.pValue$Group == "C>T"), ]$p.value < 0.01) {
#             CTpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "C>T"), ]$p.value, 
#                       digits = 3)), "**", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else if (df.pValue[which(df.pValue$Group == "C>T"), ]$p.value < 0.05) {
#             CTpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "C>T"), ]$p.value, 
#                       digits = 3)), "*", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else {
#             CTpV <- grid::textGrob(expression(""), 
#                                    gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         }
#     }
#     
#     if (is.na(df.pValue[which(df.pValue$Group == "T>A"), ]$p.value)){
#         TApV <- grid::textGrob(expression(""), 
#                                gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#     }else{
#         if (df.pValue[which(df.pValue$Group == "T>A"), ]$p.value < 0.01) {
#             TApV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "T>A"), ]$p.value, 
#                       digits = 3)), "**", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else if (df.pValue[which(df.pValue$Group == "T>A"), ]$p.value < 0.05) {
#             TApV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "T>A"), ]$p.value, 
#                       digits = 3)), "*", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else {
#             TApV <- grid::textGrob(expression(""), 
#                                    gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         }
#     }
#     
#     if (is.na(df.pValue[which(df.pValue$Group == "T>C"), ]$p.value)){
#         TCpV <- grid::textGrob(expression(""), 
#                                gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#     }else{
#         if (df.pValue[which(df.pValue$Group == "T>C"), ]$p.value < 0.01) {
#             TCpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "T>C"), ]$p.value, 
#                       digits = 3)), "**", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else if (df.pValue[which(df.pValue$Group == "T>C"), ]$p.value < 0.05) {
#             TCpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "T>C"), ]$p.value, 
#                       digits = 3)), "*", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else {
#             TCpV <- grid::textGrob(expression(""), 
#                                    gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         }
#     }
#     
#     if (is.na(df.pValue[which(df.pValue$Group == "T>G"), ]$p.value)){
#         TGpV <- grid::textGrob(expression(""), 
#                                gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#     }else{
#         if (df.pValue[which(df.pValue$Group == "T>G"), ]$p.value < 0.01) {
#             TGpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "T>G"), ]$p.value, 
#                       digits = 3)), "**", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else if (df.pValue[which(df.pValue$Group == "T>G"), ]$p.value < 0.05) {
#             TGpV <- grid::textGrob(paste("p = ", as.character(
#                 round(df.pValue[which(df.pValue$Group == "T>G"), ]$p.value, 
#                       digits = 3)), "*", sep=""), 
#                 gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         } else {
#             TGpV <- grid::textGrob(expression(""), 
#                                    gp=grid::gpar(fontsize=12),vjust=0,hjust=0.75)
#         }
#     }
#     
#     
#     ## names of mutational list
#     CA <- grid::textGrob(expression(bold("C > A")),
#                          gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
#     CG <- grid::textGrob(expression(bold("C > G")),
#                          gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
#     CT <- grid::textGrob(expression(bold("C > T")),
#                          gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
#     TA <- grid::textGrob(expression(bold("T > A")),
#                          gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
#     TC <- grid::textGrob(expression(bold("T > C")),
#                          gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
#     TG <- grid::textGrob(expression(bold("T > G")),
#                          gp=grid::gpar(fontsize=12, fontface="bold"),vjust=0,hjust=1)
#     
#     group.colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF",
#                       "#3C5488FF", "#F39B7FFF", "#8491B4FF")
#     
#     sigsInputBoxplot <- dplyr::filter(sigsInputBoxplot,
#                                       !mut.frac %in% c(NA,NaN),
#                                       !mut.num %in% c(NA,NaN))
#     
#     pic <- ggplot(sigsInputBoxplot, aes(x=GroupBT, y=mut.frac, fill=Group)) + 
#         geom_boxplot(coef=100) + 
#         ggtitle(paste0(tree.mutSig$patientID)) + 
#         theme(panel.grid=element_blank(), 
#               panel.border=element_blank(), 
#               panel.background = element_blank(), 
#               legend.position='none', 
#               plot.title = element_text(size = 13, face = "bold", vjust = 0, color = "black"),
#               axis.text.x=element_text(size=10, angle = 90, vjust = 0.5, hjust=1, color = "black"), 
#               axis.ticks.x = element_blank(), 
#               axis.line.y = element_blank(),
#               axis.ticks.length.y = unit(0.3, "cm"),
#               axis.text.y=element_text(size=10, color = "black")) + 
#         annotate("segment", x = 0.3, xend = 0.3, y = 0, yend = 100, size = 0.6) + 
#         ## background colors
#         geom_rect(aes(xmin=0.5, xmax=2.5, ymin=0, ymax=100),
#                   fill="#fce7e4", alpha=0.15) + 
#         geom_rect(aes(xmin=2.5, xmax=4.5, ymin=0, ymax=100),
#                   fill="#ecf8fa", alpha=0.25) + 
#         geom_rect(aes(xmin=4.5, xmax=6.5, ymin=0, ymax=100),
#                   fill="#dbfff9", alpha=0.05) + 
#         geom_rect(aes(xmin=6.5, xmax=8.5, ymin=0, ymax=100),
#                   fill="#e4e8f3", alpha=0.08) + 
#         geom_rect(aes(xmin=8.5, xmax=10.5, ymin=0, ymax=100),
#                   fill="#fdefeb", alpha=0.15) + 
#         geom_rect(aes(xmin=10.5, xmax=12.5, ymin=0, ymax=100),
#                   fill="#e5e8ef", alpha=0.1) + 
#         geom_boxplot(coef=100) + 
#         ## color setting
#         scale_fill_manual(values=group.colors) + 
#         ## axis setting
#         scale_x_discrete(name = "", labels=c( "Trunk","Branch", "Trunk", "Branch",  
#                                               "Trunk", "Branch",  "Trunk", "Branch", 
#                                               "Trunk", "Branch",  "Trunk","Branch")) + 
#         scale_y_continuous(name = "Mutation fraction (%)", limits=c(-5, 100), breaks=seq(0, 100, 25)) + 
#         coord_cartesian(ylim = c(0,100), expand = TRUE) + 
#         ## x axis bar
#         geom_rect(aes(xmin=0.5, xmax=2.5, ymin=-5, ymax=-0.5),
#                   fill=group.colors[1], alpha=1) + 
#         geom_rect(aes(xmin=2.5, xmax=4.5, ymin=-5, ymax=-0.5),
#                   fill=group.colors[2], alpha=0.25) + 
#         geom_rect(aes(xmin=4.5, xmax=6.5, ymin=-5, ymax=-0.5),
#                   fill=group.colors[3], alpha=0.05) + 
#         geom_rect(aes(xmin=6.5, xmax=8.5, ymin=-5, ymax=-0.5),
#                   fill=group.colors[4], alpha=0.08) + 
#         geom_rect(aes(xmin=8.5, xmax=10.5, ymin=-5, ymax=-0.5),
#                   fill=group.colors[5], alpha=0.15) + 
#         geom_rect(aes(xmin=10.5, xmax=12.5, ymin=-5, ymax=-0.5),
#                   fill=group.colors[6], alpha=0.1) + 
#         ## Mutational Type Labels
#         annotation_custom(grob = CA,  xmin = 1, xmax = 3, ymin = -8.5, ymax = -0) + 
#         annotation_custom(grob = CG,  xmin = 3, xmax = 5, ymin = -8.5, ymax = -0) + 
#         annotation_custom(grob = CT,  xmin = 5, xmax = 7, ymin = -8.5, ymax = -0) + 
#         annotation_custom(grob = TA,  xmin = 7, xmax = 9, ymin = -8.5, ymax = -0) + 
#         annotation_custom(grob = TC,  xmin = 9, xmax = 11, ymin = -8.5, ymax = -0) + 
#         annotation_custom(grob = TG,  xmin = 11, xmax = 13, ymin = -8.5, ymax = -0) + 
#         ## Mutational Type p value of wilcox.test
#         
#         annotation_custom(grob = CApV,
#                           xmin = 0.7, xmax = 3,
#                           ymin = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>A"), ]),
#                           ymax = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>A"), ])) +
#         annotation_custom(grob = CGpV,
#                           xmin = 2.7, xmax = 5,
#                           ymin = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>G"), ]),
#                           ymax = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>G"), ])) +
#         annotation_custom(grob = CTpV,  xmin = 4.7, xmax = 7, ymin = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>T"), ]), ymax = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "C>T"), ])) +
#         annotation_custom(grob = TApV,  xmin = 6.7, xmax = 9, ymin = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>A"), ]), ymax = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>A"), ])) +
#         annotation_custom(grob = TCpV,  xmin = 8.7, xmax = 11, ymin = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>C"), ]), ymax = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>C"), ])) +
#         annotation_custom(grob = TGpV,  xmin = 10.7, xmax = 13, ymin = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>G"), ]), ymax = getYmax(sigsInputBoxplot[which(sigsInputBoxplot$Group == "T>G"), ]))
#     #message("Branch-trunk plot generation done!")
#     
#     return(pic)
# }