ccfPair <- function(mafData, pairByType = FALSE){
  
   mafData <- as.data.frame(mafData)
   patientid <- unique(mafData$Patient_ID)
   if(pairByType){
      types <- unique(mafData$Tumor_Type)
      if(length(types) < 2){
         message(paste0("Warnings: Only one tumor type was found of ",patientid,". It you want to compare CCF between regions, pairByType should be set as FALSE"))
         return(NA)
      }
      
      ## get average CCF
      mafData <- mafData %>% 
         # dplyr::mutate(Type_Average_CCF = dplyr::if_else(
         #     Type_Average_CCF > 1,
         #     1,
         #     Type_Average_CCF
         # )) %>% 
          dplyr::mutate(CCF = Type_Average_CCF) %>%
          dplyr::filter(!is.na(CCF))
      pairs <- combn(length(types), 2, simplify = FALSE)  
   }else{
      samples <- unique(mafData$Tumor_Sample_Barcode)
      if(length(samples) < 2){
         message(paste0(unique(mafData$Patient_ID), " have less than 2 samples, return NA"))
         return(NA)
      }
      
      pairs <- combn(length(samples), 2, simplify = FALSE)  
   }
   
   ccf.pair.list <- list()
   pair.name <- c()
   i <- 1
   for (pair in pairs){
      if(pairByType){
         S1 <- types[pair[1]]
         S2 <- types[pair[2]]  
      }else{
         S1 <- samples[pair[1]]
         S2 <- samples[pair[2]]
      }
      if(pairByType){
         ccf.pair <- mafData %>% dplyr::filter(Tumor_Type %in% c(S1, S2)) %>%
            tidyr::unite("Mut_ID", c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), 
                         sep = ":", remove = FALSE) %>% 
            tidyr::unite("Mut_ID2", c("Tumor_Type","Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), 
                         sep = ":", remove = FALSE) %>% 
            dplyr::distinct(Mut_ID2, .keep_all = TRUE) %>%
            dplyr::select(Tumor_Type, Hugo_Symbol, Mut_ID, CCF) %>% 
            tidyr::spread(key = Tumor_Type, value = CCF)
            # dplyr::rename("sample1" = S1, "sample2" = S2) %>% 
            # dplyr::filter(!is.na(S1),!is.na(S2))
         
      }else{
         ccf.pair <- mafData %>% dplyr::filter(Tumor_Sample_Barcode %in% c(S1, S2)) %>%
            tidyr::unite("Mut_ID", c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), 
                         sep = ":", remove = FALSE) %>% 
            dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, Mut_ID, CCF) %>% 
            tidyr::spread(key = Tumor_Sample_Barcode, value = CCF)
            # dplyr::rename("sample1" = S1, "sample2" = S2) %>% 
            # dplyr::filter(!is.na(S1),!is.na(S2)) 
      }
      if(nrow(ccf.pair) == 0){
         message(paste0("Warning: No shared mutaions were detected between ",S1, " and ", S2) )
         p[[i]] <- NA
         next()
      }
      v1 <- as.vector(ccf.pair[,3])
      v2 <- as.vector(ccf.pair[,4])
      ccf.pair <- ccf.pair[(!is.na(v1)&!is.na(v2)),]
      ccf.pair.list[[i]] <- ccf.pair
      pair.name <- c(pair.name, paste(S1,"-",S2, sep = ""))
      i =  i+1
   }
   names(ccf.pair.list) <- pair.name
   return(ccf.pair.list)  
}


# p[[i]] <- ggplot2::ggplot(data= ccf.pair, aes(x=sample1, y=sample2)) +
#    xlab(S1) +
#    ylab(S2) +
#    theme(
#       # axis.text = element_text(size = 12, color = "black"),
#       axis.text = element_text(size = 10, color = "black"),
#       axis.title = element_text(size = 13, color = "black"),
#       legend.text = element_text(size = 10),
#       legend.title = element_text(size = 12),
# 
#       axis.line.x.top = element_blank(),
#       axis.line.y.right = element_blank(),
#       axis.ticks.length.y.right = element_blank(),
#       axis.ticks.length.x.top = element_blank(),
#       panel.background = element_rect(color = "black",size = 0.5),
#       #panel.background = element_rect(colour = "black", size = 0.4),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       plot.title = element_text(hjust = 0.5)
#    )+
#    guides(shape = guide_legend(override.aes = list(size = 0.2))) +
#    # theme_bw() +
#    coord_fixed() +
#    scale_x_continuous(expand = c(0,0), limits = c(0,1))+
#    scale_y_continuous(expand = c(0,0), limits = c(0,1))
# # labs(x = S1, y = S2) +
# if(pairByType){
#    p[[i]] <- p[[i]] + ggtitle(paste("CCF density plot in paired tumor types:\n ", S1, " vs ", S2, sep = ""))
# 
# }
# else{
#    p[[i]] <- p[[i]] + ggtitle(paste("CCF density plot in paired samples:\n ", S1, " vs ", S2, sep = ""))
# }
# # if(show.density){
# p[[i]] <- p[[i]] + stat_density_2d(aes(fill = ..density..), geom = 'tile',contour = FALSE) +
#     # scale_fill_gradient(low = "#8491B499",high = "#E64B35FF")
#     scale_fill_gradientn(colours = c("white",
#                                      "#bcbddc",
#                                      "#9e9ac8",
#                                      "#807dba",
#                                      "#fc9272",
#                                      "#ef3b2c"),
#                          name = "Density") + 
#      geom_point(size = 0.5)

## set table for gene
# genes.table <- ccf.pair %>% 
#    dplyr::rowwise() %>% 
#    dplyr::filter(any(strsplit(Hugo_Symbol,",|;")[[1]] %in% geneList)) %>%
#    as.data.table()
# 
# if(nrow(genes.table) > 0){
#    genes.name <- unlist(lapply(genes.table$Hugo_Symbol, function(x){
#       ns <- strsplit(x,",|;")[[1]]
#       idx <- which(ns %in% geneList)[1]
#       return(ns[idx])
#    }))
#    genes.table$Hugo_Symbol <- genes.name
#    p[[i]] <- p[[i]] +
#       ## add point
#       geom_point(data = genes.table,
#                  aes(x = sample1, y = sample2),
#                  shape = 24,
#                  size = 2,
#                  fill = "red")+
#       ## label gene
#       geom_text_repel(data = genes.table,
#                       aes(x = sample1, y = sample2, label = Hugo_Symbol),
#                       size = 3.5,
#                       force = 10)


# # use smooth scatter
# m <- matrix(c(ccf.pair[,3], ccf.pair[,4]),ncol = 2)
# par(mfrow=c(1,1),mar=c(3,3,3,3),mgp=c(2,1,0),pty="s")
# smoothScatter(m, xlim = c(0,1), ylim = c(0,1),
#               colramp = colorRampPalette(c("white",brewer.pal(9,"BuPu")), space="rgb"),
#               main = patientid,
#               xlab = S1,
#               ylab = S2)
# 
# ## label muation on geneList
# if(!is.null(geneList)){
#     genes.table <- ccf.pair %>%
#         dplyr::rowwise() %>%
#         dplyr::filter(any(strsplit(Hugo_Symbol,",|;")[[1]] %in% geneList)) %>%
#         as.data.table()
#     if(nrow(genes.table) > 0){
#         points(genes.table[,.(sample1,sample2)],cex=0.6,col=2,pch=2)
#         text(genes.table[,.(sample1,sample2)],cex=0.7,pos=1,genes.table$Hugo_Symbol)
#     }else{
#         message(paste0("Warning: None of genes map to pair samples ",S1, " vs ", S2) )
#     }
# }
# x <- recordPlot()