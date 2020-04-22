getHeatmapMatrix <- function(mafData){
   matrix.list <- NA
   binary.mat <- getMutMatrix(mafData, use.ccf = FALSE)
   if("CCF" %in% colnames(mafData)){
      ccf.mat <- getMutMatrix(mafData, use.ccf = TRUE)
   }else{
      ccf.mat <- matrix() 
   }
   return(list(binary.mat, ccf.mat))
}

plotHeatmap <- function(binary.mat,
                        ccf.mat,
                        use.ccf = FALSE,
                        show.class.label = TRUE,
                        geneList = NULL,
                        plot.geneList = FALSE,
                        show.gene = FALSE,
                        show.geneList = TRUE,
                        mut.threshold = 50){
   mut_sort <- binary.mat
   ccf_sort <- ccf.mat
   
   if(all(is.na(ccf_sort)) & use.ccf){
       stop(paste0("Error :Heatmap requires CCF data when use.ccf is True"))
   }
   
   ## delete "NORMAL" 
   if(!1 %in% mut_sort[,"NORMAL"]){
      mut_sort <- mut_sort[,which(colnames(mut_sort)!= "NORMAL")]
      if(use.ccf){
         ccf_sort <- ccf_sort[,which(colnames(ccf_sort)!= "NORMAL")]
      }
   }
   
   mat <- mut_sort
   type  <- "Mutation"
   if(use.ccf){
      type <- "CCF"
      mat <- ccf_sort
   }
   shared.num <- ncol(mut_sort)
   mutation.classes <- apply(mut_sort,1,function(x,shared.num){
      if(sum(x) == shared.num){
         return("Shared")
      }
      else if(sum(x) == 1){
         return("Private")
      }
      else{
         return("P-shared")
      }
   },shared.num = shared.num)
   mut_dat <- heatmap_input(mat,
                            type = type,
                            mutation.classes = mutation.classes,
                            geneList = geneList,
                            plot.geneList = plot.geneList)
   
   ## Do not the row name if the number of mutations is greater than mut.threshold
   mut.num <- nrow(mut_dat)/length(unique(mut_dat$sample))
   
   if(show.gene == TRUE|(!is.null(geneList) & show.geneList == TRUE))
   if(mut.num >= mut.threshold){
       message("Warning: the number of mutations is ", mut.num,
               " which is greater than mut.threthold. Let mut.threshold be larger than ", mut.num," if you want to show gene")
      show.gene = FALSE
      show.geneList = FALSE 
   }
   
   ## the num of each mutation class
   classes.num <- c()
   classes.sum <- length(mut_dat$class)/length(unique(mut_dat$sample))
   
   ## set table for annotation bar
   annotation.bar <- data.frame()
   annotation.bar.width <- (max(mut_dat$xmax)- max(mut_dat$xmin))*0.5
   if(!show.class.label){
      annotation.bar.width <- annotation.bar.width/3
   }
   
   ## position of annotation bar
   xmax <- min(mut_dat$xmin)
   xmin <- xmax - annotation.bar.width
   ymin <- min(mut_dat$ymin)
   classes.level <- unique(mut_dat$class)
   for(class in classes.level){
      ymax <- max(mut_dat[mut_dat$class == class,]$ymax)
      sub <- c(xmin, xmax, ymin, ymax)
      annotation.bar <- rbind(annotation.bar, sub)
      ymin <- ymax
      ## percentage of class
      class.num <- length(which(mut_dat$class == class))/(classes.sum*length(unique(mut_dat$sample)))
      percentage <- paste0("\n",round(class.num,3)*100,"%")
      classes.num[class] <- percentage
   }
   colnames(annotation.bar) <- c("xmin","xmax","ymin","ymax")
   
   ## color of annotation bar
   class.all.colors <- c( "Shared" = "#7fc97f",
                          "P-shared" = "#fdc086",
                          "Private" = "#E64B35FF" )
   class.colors <- class.all.colors[classes.level]
   
   ## label of annotation bar
   classes.label <- paste0(classes.level,classes.num)
   
   p_basic <- ggplot() +
      labs(x = "", y = "") + theme_bw() +
      theme(panel.border = element_blank()) +
      theme(axis.text.y = element_blank())+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      
      ## set label for axis X
      theme(axis.text.x.top = element_text(angle = 90,
                                           hjust = 0,
                                           size = 9,
                                           color = "black",
                                           margin = margin(b = -15)))+
      scale_x_continuous(
         breaks = unique(mut_dat$xmin) + (unique(mut_dat$xmax) - unique(mut_dat$xmin))/2,
         labels = unique(mut_dat$sample),
         position = "top")+
      
      
      theme(axis.ticks = element_blank()) +
      theme(legend.title = element_text(color = "black")) +
      theme(legend.text = element_text( color = "black")) +
      theme(legend.position = "right" )+
      
      ## annotation bar
      geom_rect(data = annotation.bar,
                mapping = aes(xmin = xmin,xmax = xmax,ymin = ymin, ymax = ymax),
                fill = class.colors)
   if(show.class.label){
      p_basic <- p_basic + 
         geom_text(data = annotation.bar,
                   mapping = aes(x = xmin + (xmax-xmin)/2,
                                 y = ymin + (ymax-ymin)/2),
                   label = classes.label,
                   angle = 90)
   }
   if(use.ccf){
      p <- p_basic + 
         geom_rect(data = mut_dat,
                   mapping = aes(xmin = xmin,xmax = xmax,ymin = ymin, ymax = ymax,fill = CCF))+
         scale_fill_gradient(low = "#deebf7", high = "#08306b", na.value="black", limit=c(0, 1))
      
      #ggsave(paste(patientID, "_mut_CCF.pdf", sep = ""), p, width = 4.5, height = 6.5)
   }else if(!use.ccf){
      mut_dat$Mutation <- as.character(mut_dat$Mutation)
      p <- p_basic + 
         geom_rect(data = mut_dat,
                   mapping = aes(xmin = xmin,xmax = xmax,ymin = ymin, ymax = ymax,fill = Mutation))+
         scale_fill_manual(values = c("#deebf7", "#08306b"))
      #ggsave(paste(patientID, "_mut.pdf", sep = ""), p, width = 4.5, height = 6.5)
   }
   if(is.null(geneList) & show.gene){
      breaks.gene <- unique(mut_dat$ymin + (mut_dat$ymax - mut_dat$ymin)/2)
      p <- p + scale_y_continuous(breaks = breaks.gene,
                                  labels = mut_dat[mut_dat$sample==unique(mut_dat$sample)[1],]$Gene,
                                  position = "right")+
         theme(axis.text.y.right = element_text(size = 9,
                                                colour = "black",
                                                face = "italic",
                                                margin = margin(l = -15),
                                                hjust = 0))
   }else if(!is.null(geneList)){
      if(plot.geneList & show.geneList){
         y.breaks <- unique(mut_dat$ymin + (mut_dat$ymax - mut_dat$ymin)/2)
         y.labels <- unique(mut_dat$Gene)
         p <- p + 
            scale_y_continuous(breaks = y.breaks,
                               labels = y.labels,
                               position = "right") +
            theme(axis.text.y.right = element_text(size = 9,
                                                   colour = "black",
                                                   face = "italic",
                                                   margin = margin(l = -15),
                                                   hjust = 0))
      }
      else if(!plot.geneList & show.geneList){
         gene.pos <- unique(which(mut_dat$Gene != "nogene")) 
         y.breaks <- unique((mut_dat$ymin + (mut_dat$ymax - mut_dat$ymin)/2)[gene.pos])
         y.labels <- unique(mut_dat$Gene[gene.pos]) 
         p <- p + 
            scale_y_continuous(breaks = y.breaks,
                               labels = y.labels,
                               position = "right") +
            theme(axis.text.y.right = element_text(size = 9,
                                                   colour = "black",
                                                   face = "italic",
                                                   margin = margin(l = -15),
                                                   hjust = 0))
      }
   }
   
   return(p)
}

heatmap_input <- function(mat,
                          type,
                          mutation.classes,
                          geneList,
                          plot.geneList){
   mat <- as.data.frame(mat)
   sample.num <- ncol(mat)
   mat$class <- factor(mutation.classes, levels = unique(mutation.classes)) 
   ## get gene name
   genes <- lapply(rownames(mat),function(x){
      s <- strsplit(x,":")[[1]][1]
      return(s)
   }) %>% unlist()
   mat$Gene <- as.character(genes) 
   
   ## flit mutation in gene list
   if(!is.null(geneList)){
      if(plot.geneList){
         mat <- mat  %>%
            dplyr::rowwise() %>%
            dplyr::filter(any(strsplit(Gene, ",|;")[[1]] %in% geneList)) %>%
            as.data.frame() 
      }else{mat <- mat  %>%
         dplyr::rowwise() %>%
         dplyr::mutate(Gene = dplyr::if_else(
            any(strsplit(Gene, ",|;")[[1]] %in% geneList),
            Gene,
            "nogene"
         )) %>%
         as.data.frame() 
      }
   }
   mut.num <- nrow(mat)
   mat$mutation <- 1:nrow(mat)
   mat <- dplyr::arrange(mat,class)
   
   ## cumsum postion of axis y
   mat$ymin <- cumsum(c(0, rep(2, mut.num-1)))
   mat$ymax <- mat$ymin + 1.5
   value.name <- "Mutation"
   if(type == "CCF"){
      value.name <- "CCF"
   }
   mut_dat <- reshape2::melt(mat,
                             id.vars = c("ymin","ymax","mutation","class", "Gene"),
                             variable.name = "sample",
                             value.name = value.name)
   mut_dat$xmin <- rep(cumsum(c(0, rep(0.5, sample.num-1))) ,each = nrow(mat))
   mut_dat$xmax <- mut_dat$xmin + 0.49
   mut_dat$sample <- factor(mut_dat$sample, levels =unique(mut_dat$sample))
   return(mut_dat)
}

