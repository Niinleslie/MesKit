plotTree <- function(phyloTree,
                     treeData = NULL,
                     branchCol = "mutType",
                     show.bootstrap = TRUE,
                     signaturesRef = "cosmic_v2",
                     min.mut.count = 15,
                     min.ratio = 1/20,
                     uncommon.col = "red",
                     compare.tree.name = "",
                     use.tumorSampleLabel = FALSE,
                     show.scale.bar = FALSE,
                     scale.bar.x = NULL,
                     scale.bar.y = NULL){
    
    if(min.ratio <= 0|min.ratio > 1){
        stop("min.ratio should be within (0,1]")
    }
    
    ## adjust length of branches by min.ratio
    tree <- getTree(phyloTree)
    min.len <- max(tree$edge.length)*min.ratio
    tree$edge.length[tree$edge.length < min.len] <- min.len
    phyloTree <- new('phyloTree',
                      patientID = getPhyloTreePatient(phyloTree),
                      tree = tree, 
                      binary.matrix = getBinaryMatrix(phyloTree),
                      ccf.matrix = getCCFMatrix(phyloTree), 
                      mut.branches = getMutBranches(phyloTree),
                      branch.type = getBranchType(phyloTree),
                      ref.build = getPhyloTreeRef(phyloTree),
                      bootstrap.value = getBootstrapValue(phyloTree),
                      method = getTreeMethod(phyloTree),
                      tsb.label = getPhyloTreeTsbLabel(phyloTree))
    
    patient <- getPhyloTreePatient(phyloTree)
    if(is.null(treeData)){
        treeData <- getTreeData(phyloTree = phyloTree,
                                branchCol = branchCol,
                                signaturesRef = signaturesRef,
                                min.mut.count = min.mut.count)
        if(identical(treeData, NA)){
            return(NA)
        }
        compare <- FALSE
    }else{
        compare <- TRUE
    }
    
    # set.seed(1234)
    ## get bootstrap value
    boot_value <- getBootstrapValue(phyloTree)
    rootLabel <- "NORMAL"
    ## plot phylotree
    samplePointsSize <- 3
    comparePointsSize <- 3.5
    sampleTextSize <- 3.5
    nodePointsSize <- 1.7
    segmentSize <- 1.5
    nodeStrokeSize <- 0.25
    sampleStrokeSize <- 1.5
    bootLabelSize <- 3
    bootPaddingSize <- 0.35
    bootLabelPaddingSize <- 0.15
    legend.text.size <- 10
    legend.title.size <- legend.text.size + 0.5
    samplesLength <- nrow(treeData[treeData$sample != "internal node",]) 
    if(samplesLength > 7){
        samplePointsSize <- 1.5 
        comparePointsSize <- 2 
        sampleTextSize <- 3
        segmentSize <- 0.8
        nodePointsSize <- 0.8
        nodeStrokeSize <- 0.15
        sampleStrokeSize <- 0.8
        bootLabelSize <- 2.5
        bootPaddingSize <- 0.1
        bootLabelPaddingSize <- 0.1
        legend.text.size <- 9
        legend.title.size <- legend.text.size + 0.5
    }
    rootNode <- treeData[treeData$sample == rootLabel,]$node
    if(length(boot_value) == 1){
        bootsData <- data.frame(x2 = 0, y2 = 0, node = rootNode, end_num = rootNode, boots = boot_value)
    }else{
        sub <- data.table::data.table(x2 = 0, y2 = 0,node = rootNode, end_num = rootNode)
        bootsData <- rbind(treeData[treeData$sample == 'internal node',] %>% 
                           dplyr::select("x2","y2","node","end_num"),sub)
        boots <- c()
        LN <- min(bootsData$node)-1
        
        boots_list <- lapply(seq_len(nrow(bootsData)), function(i){
          if(i == nrow(bootsData)){
            boots <- append(boots, boot_value[rootNode - LN])
          }else{
            boots <- append(boots, boot_value[bootsData$end_num[i] - LN])
          }
          return(boots)
        }) %>% unlist()
        
        bootsData <- cbind(bootsData, boots = boots_list)
    }
    ## initialize
    x2 <- NULL
    y2 <- NULL
    x1 <- NULL
    y1 <- NULL
    ## get the max value of X axis 
    x_max <- max(abs(treeData$x2))
    # print(treeData)
    p <- ggplot(data = treeData) + 
        ## balance the space on the left and right sides
        geom_segment(aes(x = 0, y = 0, xend = x_max, yend = 0),color = "white",size = 0.01)+
        geom_segment(aes(x = 0, y = 0, xend = -x_max, yend = 0),color = "white",size = 0.01)
    
    textAdjust <- mean(as.numeric(treeData$distance))
    
    if(use.tumorSampleLabel){
       tsb.label <- getPhyloTreeTsbLabel(phyloTree)
       if(nrow(tsb.label) == 0){
         stop("Tumor_Sample_Label was not found. Please check clinical data or let use.tumorSampleLabel be FALSE.")
       }
       ctsb <- treeData$sample[!treeData$sample %in% c("NORMAL", "internal node")]
       treeData$sample[!treeData$sample %in% c("NORMAL", "internal node")] <- lapply(as.list(ctsb),
                                                                                     function(x){
                                                                                        tsb.label[which(tsb.label$Tumor_Sample_Barcode == x),]$Tumor_Sample_Label
                                                                                        }) %>% 
                                                                              unlist()
    }
    
    if(!is.null(branchCol)){
        if(branchCol == "mutSig"){
            color_scale <- getSigColors(as.character(unique(treeData$Signature)) )
            sig_level <- levels(treeData$Signature)
            ## initialize
            Signature <- NULL
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Signature), size=segmentSize)
            p <- p + scale_color_manual(breaks = sig_level ,values = color_scale) + 
                theme(legend.title = element_text(size = legend.title.size))
        }else{
            ## sort branch tumor type 
            all.types <- unique(treeData$Mutation_Type) 
            public <- all.types[grep("Public", all.types)] 
            shared <- all.types[grep("Shared", all.types)] 
            private <- all.types[grep("Private", all.types)]
            type.level <- c(public, shared, private)
            
            ## get colors
            type_all_colors <- c("#7fc97f","#fdc086", "#E64B35FF", "#82166E",
                                 "#B77B42","#6349B7","#D5017D","#B77562",
                                 "#88A4FF", "#439F18", "#971D37","#8C9F3C")
            if(length(type.level) > length(type_all_colors)){
                left_colors <- sample(colors(),
                                      length(type.level)-length(type_all_colors),
                                      replace = FALSE)
                type_all_colors <- append(type_all_colors,left_colors)
            }else{
                type_all_colors <- type_all_colors[seq_len(length(type.level)) ] 
            }
            names(type_all_colors) <- type.level
            treeData$Mutation_Type <- factor(treeData$Mutation_Type, levels = type.level)
            ## initialize
            Mutation_Type <- NULL
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = factor(Mutation_Type)), size=segmentSize)
            ## remove legend title
            p <- p +  theme(legend.title = element_blank()) + 
                scale_color_manual(breaks = type.level, values = type_all_colors)
        }
    }else{
        if(compare){
            # p <- p + geom_point(aes(x = x1, y = y1),
            #                       color = uncommon.col,
            #                       data = treeData[treeData$is.match == "NO"],
            #                       size = samplePointsSize)
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                                  color = "black",
                                  data = treeData ,size = segmentSize)
        }else{
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = "black",
                                  data = treeData, 
                                  size=segmentSize, show.legend = FALSE)
            
        }
    }
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   panel.border = element_blank(),
                   legend.text = element_text(size = legend.text.size),
                   legend.position = 'right') + 
        scale_x_discrete(expand = expansion(add = mean(treeData$distance)))+
        coord_fixed(ratio= 1)
    p <- p + geom_text_repel(aes(x = x2, y = y2, label = sample),
                             nudge_y = textAdjust/10,
                             nudge_x = textAdjust/10,
                             segment.color = "grey",
                             segment.size = 0.25,
                             max.overlaps = 200,
                             data = treeData[(!treeData$sample %in% c("internal node",rootLabel)) &
                                                treeData$x2 >= 0, ],
                             size = sampleTextSize ,force = 10)
    p <- p + geom_text_repel(aes(x = x2, y = y2, label = sample),
                             nudge_y = textAdjust/10,
                             nudge_x = -textAdjust/10,
                             segment.color = "grey",
                             segment.size = 0.25,
                             max.overlaps = 200,
                             data = treeData[(!treeData$sample %in% c("internal node",rootLabel)) &
                                                treeData$x2 < 0, ],
                             size = sampleTextSize ,force = 10)
    ## label NORMAL
    p <- p + geom_text(aes(x = x2,y = y2-textAdjust/5),
                       label = rootLabel,
                       data = treeData[treeData$sample == rootLabel,], 
                       size = sampleTextSize)
    
    if(compare & is.null(branchCol)){
       unmatch_dat <- treeData[treeData$is.match == "NO"&treeData$sample == 'internal node',]
       if(nrow(unmatch_dat) > 1){
          unmatch_label <- paste0("Clades absent in ", compare.tree.name)
       }else if(nrow(unmatch_dat) == 1){
          unmatch_label <- paste0("Clade absent in ", compare.tree.name)
       }
       
       if(nrow(unmatch_dat) > 0){
          rootLabel <- "NORMAL"
          numNORMAL <- which(tree$tip.label == rootLabel)
          tree <- ape::root(tree, tree$tip.label[numNORMAL])
          treeEdge <- data.table::data.table(node = tree$edge[,1], endNum = tree$edge[,2])
          rootRow <- which(treeEdge$endNum == numNORMAL)
          rootNode <- treeEdge$node[rootRow]
          
          node_point_list <- c()
          for(i in seq_len(nrow(unmatch_dat))){
             point <- unmatch_dat$end_num[i]
             for(i1 in seq_len(length(tree$tip.label))){
                tip <- tree$tip.label[i1]
                if(tip == rootLabel){
                   next
                }
                start <- i1
                tip_point_list <- c(i1)
                while(TRUE){
                   end <- as.numeric(treeEdge[treeEdge$endNum == start,]$node)
                   tip_point_list <- c(tip_point_list, end)
                   if(end == point){
                      node_point_list <- unique(c(node_point_list, tip_point_list))
                      node_point_list <- node_point_list[node_point_list!=point]
                      break
                   }
                   if(end == rootNode){
                      break
                   }
                   start <- end
                }
             }
          }
          
          p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                                color = uncommon.col,
                                data = treeData[treeData$end_num %in% node_point_list,],
                                size = segmentSize)
          p <- p + geom_point(aes(x = x2, y = y2, color = is.match),
                              # color = uncommon.col,
                              data = unmatch_dat,
                              size = comparePointsSize) + 
             scale_color_manual(label = unmatch_label,values = uncommon.col) + 
             theme(legend.title = element_blank())
       }
    }
    
    p <- p + geom_point(aes(x = x2,y = y2),
                        data = treeData[treeData$sample == 'internal node',],
                        size = nodePointsSize, color = "#8c510a", fill = "white", shape = 21,
                        stroke = nodeStrokeSize)
    p <- p + geom_point(aes(x = x2, y = y2), 
                        data = treeData[treeData$sample != 'internal node',],
                        size = samplePointsSize,color = "#67001F", fill = 'white', shape = 21, stroke = sampleStrokeSize)
    
    
    
    Nd <- treeData[sample == rootLabel,]$distance
    if(length(Nd)!=0){
       if(Nd != 0){
          p <- p + geom_point(aes(x =0 , y = 0), size = nodePointsSize, color = "#8c510a",
                              fill = 'white', shape = 21, stroke = nodeStrokeSize)
       }
    }
    
    if(show.bootstrap){
        p <- p + geom_label_repel(aes(x = x2, y = y2,label = boots),
                                  data = bootsData,
                                  # nudge_y = textAdjust/6,
                                  # fontface = 'bold', 
                                  size = bootLabelSize, 
                                  box.padding = unit(bootPaddingSize, "lines"),
                                  label.padding = bootLabelPaddingSize,
                                  segment.colour = "grey50", 
                                  segment.size = 0.25, 
                                  max.overlaps = 200,
                                  force = 5)
    }
    # if(compare){
    #   is.match <- NULL
    #     p <- p + geom_label_repel(aes(x = x1 + (x2-x1)/2 , y = y1 + (y2 - y1)/2, label = is.match),
    #                               data = treeData[treeData$is.match != "NO",],
    #                               fontface = 'bold', 
    #                               size = bootLabelSize, box.padding = unit(bootPaddingSize, "lines"), point.padding = unit(0.5, "lines"),
    #                               segment.colour = "grey50", segment.size = 0.5, force = 5)
    # }
    # 
    tree.title <- patient

    if(show.scale.bar){
       mean_len <- ceiling(mean(treeData$distance)) 
       if(!is.null(scale.bar.x)){
          x_bar <- as.numeric(scale.bar.x)
       }else{
          x_bar <- min(treeData$x2) + min(treeData$x2)/2
       }
       
       if(!is.null(scale.bar.y)){
          y_bar <- as.numeric(scale.bar.y)
       }else{
          y_bar <- (max(treeData$y2) - min(treeData$y2))/2
       }
       
       
       df_bar <- data.frame(
          x = x_bar,
          xend = x_bar,
          y = y_bar - mean_len/2,
          yend = y_bar + mean_len/2
      )
       p <- p + 
          geom_segment(
             aes(x = x, y = y, xend = xend, yend = yend),
             data = df_bar,
             size = 0.6,
             color = "black"
          ) + 
          geom_text(
             aes(x = x, y = yend + mean_len/10),
             data = df_bar,
             label = as.character(mean_len),
             size = 4,
             color = "black"
          )
    }
    if(!compare){
       p <- p + 
          ggtitle(tree.title)+
          theme(plot.title = element_text(face = "bold",colour = "black", hjust = 0.5,size = 13.5))
    }

    return(p)
}

getTreeData <- function(phyloTree = NULL,
                        branchCol = "mutSig",
                        signaturesRef = "cosmic_v2",
                        min.mut.count = 15,
                        compare = FALSE){
   tree <- getTree(phyloTree)
   rootLabel <- "NORMAL"
   tree <- ape::root(tree, tree$tip.label[which(tree$tip.label == rootLabel)])
   treeEdge <- data.table::data.table(node = tree$edge[,1], endNum = tree$edge[,2], length = tree$edge.length)
   numNORMAL <- which(tree$tip.label == rootLabel)
   #the position of NORMAL in edge 
   rootRow <- which(treeEdge$endNum == numNORMAL)
   # the  Node connected to NORMAL
   rootNode <- treeEdge$node[rootRow]
   rootEdge <- treeEdge$length[rootRow]
   lb <- labelBranch(tree)
   branchLabel <- lb[[1]]
   subnumlist <- lb[[2]]
   mainTrunk <- calMainTrunk(tree, treeEdge, rootNode = rootNode, rootLabel = rootLabel)
   ccn <- calChildNodeNum(tree, treeEdge = treeEdge, mainTrunk, rootNode = rootNode)
   numList <- ccn[[1]]
   pointsList <- ccn[[2]]
   nodeNoOnTree <- ccn[[3]]
   nodeOnTree <- ccn[[4]]
   if(length(mainTrunk) != 0){
      gta <- getNodeAngle(tree = tree, treeEdge = treeEdge,
                            mainTrunk = mainTrunk, numList = numList,
                            pointsList = pointsList, rootNode = rootNode,
                            horizon = pi/2, W = pi, numNORMAL = numNORMAL)
      adjacentWs <- gta[[1]]
      adjacentAngles <- gta[[2]]
      adjacentPoints <- gta[[3]]
   }
   #X1, y1 are the starting point, x2 and y2 are the end point.Horizon stands for the Angle between the line and the positive x axis;W stands for the Angle occupied by branches (explained in the paper)
   #Distance means distance;Node is the internal node with the start;End_num is the number of each node in edge;
   treeData <- data.table::data.table('x1'= 0, 'y1' = 0,
                                      'x2' = 0, 'y2' = 0,
                                      'horizon' = pi/2, 'W' = pi,
                                      'angle' = 0,'distance' = 0,
                                      'node' = rootNode, 'end_num' = rootNode)
   treeData <- setPhyloTree(tree = tree, treeEdge = treeEdge, treeData = treeData, rootNode = rootNode,
                            mainTrunk = mainTrunk, adjacentPoints = adjacentPoints,
                            adjacentWs = adjacentWs, adjacentAngles = adjacentAngles,
                            horizon = pi/2, W = pi)
   # subTrees <- ape::subtrees(tree)
   # subroot <- unlist(lapply(subTrees, function(x){return(x$name)}))
   if(length(nodeNoOnTree)!=0){
      t <- 1
      while(TRUE){
         node <- nodeNoOnTree[t]
         subnodes <- subnumlist[[node]]
         subEdge <- treeEdge[treeEdge$node %in% subnodes, ]
         rootNode <- node
         numNORMAL <- NULL
         mainTrunk <- calMainTrunk(tree, subEdge, rootNode = rootNode, rootLabel = "")
         ccn <- calChildNodeNum(tree, subEdge, mainTrunk, rootNode = rootNode, ft = TRUE)
         numList <- ccn[[1]]
         pointsList <- ccn[[2]]
         nodeNoOnTree <- append(nodeNoOnTree,ccn[[3]])
         nodeOnTree <- append(nodeOnTree,ccn[[4]])
         horizon <- treeData[treeData$end_num == rootNode,]$angle
         W <- treeData[treeData$end_num == rootNode,]$W
         gta <- getNodeAngle(tree = tree, treeEdge = subEdge,
                               mainTrunk = mainTrunk, numList = numList,
                               pointsList = pointsList, rootNode = rootNode,
                               horizon = horizon, W = W, numNORMAL = NULL)
         adjacentWs <- gta[[1]]
         adjacentAngles <- gta[[2]]
         adjacentPoints <- gta[[3]]
         treeData <- setPhyloTree(tree = tree, treeEdge = subEdge, treeData = treeData, rootNode = rootNode,
                                  mainTrunk = mainTrunk, adjacentPoints = adjacentPoints,
                                  adjacentWs = adjacentWs, adjacentAngles = adjacentAngles,
                                  horizon = horizon, W = W)
         if(length(nodeOnTree) == tree$Nnode){
            break
         }
         t <- t + 1
      }
   }
   ## bind NORMAL
   rootNode <- treeData[1,]$node
   treeData[1,]$end_num <- which(tree$tip.label == "NORMAL")
   treeData[1,]$y2 <- -rootEdge
   treeData[1,]$distance <- rootEdge
   treeData$sample <- ""
   
   
   
   sample_list <- vapply(treeData$end_num, function(e){
     if(e > length(tree$tip.label)){
       return("internal node")
     }else{
       return(tree$tip.label[e])
     }
   }, FUN.VALUE = character(1))
   
   treeData$sample <- sample_list
   
   # for(i in seq_len(nrow(treeData))){
   #    if(treeData$end_num[i] > length(tree$tip.label)){
   #       treeData$sample[i] <- "internal node"
   #    }
   #    else{
   #       pos <- treeData$end_num[i]
   #       treeData$sample[i] <- tree$tip.label[pos]
   #    }
   # }
   if(nrow(treeData) == 3){
      angleList <- c(pi/3, 2*pi/3)
      rows <- which(treeData$sample != "NORMAL")
      treeData[rows, ]$angle <- angleList
      treeData[rows, ]$x2 <- treeData[rows, ]$distance*cos(treeData[rows, ]$angle)
      treeData[rows, ]$y2 <- treeData[rows, ]$distance*sin(treeData[rows, ]$angle)
   }
   if(nrow(treeData[treeData$x2 > 0,]) == 0){
      maxy <- which.max(treeData$y2)
      angle <- treeData$angle[maxy] - pi/18
      treeData$angle[maxy] <- angle
      treeData$y2[maxy] <- treeData$distance[maxy]*sin(angle) + treeData$y1[maxy]
      treeData$x2[maxy] <- treeData$distance[maxy]*cos(angle) + treeData$x1[maxy]
   }else if(nrow(treeData[treeData$x2 < 0,]) == 0){
      maxy <- which.max(treeData$y2)
      angle <- treeData$angle[maxy] + pi/18
      treeData$angle[maxy] <- angle
      treeData$y2[maxy] <- treeData$distance[maxy]*sin(angle) + treeData$y1[maxy]
      treeData$x2[maxy] <- treeData$distance[maxy]*cos(angle) + treeData$x1[maxy]
   }
   
   ## label represents the common evolution path of samples
   treeData$label <- ""
   label_list <- vapply(seq_len(nrow(treeData)), function(i){
     if(treeData$sample[i] == "NORMAL"){
       return(branchLabel[[rootNode]])
     }else{
       if(treeData$end_num[i] > length(tree$tip.label)){
         return(branchLabel[[treeData$end_num[i]]])
       }
       else{
         return(as.character(treeData$sample[i])) 
       }
     }
   }, FUN.VALUE = character(1))
   treeData$label <- label_list
   
   if(!is.null(branchCol) & !compare){
       ## add signature
      if(branchCol == "mutSig"){
          tri_matrix <- triMatrix(phyloTree,level = "4")
          fit_out <- fitSignatures(tri_matrix,
                                   signaturesRef = signaturesRef,
                                   min.mut.count = min.mut.count)
          if(is.na(fit_out)){
              message("The number of mutations in all branches are less than min.mut.count in patient ",
                      phyloTree@patientID)
              return(NA)
          }
          cos_sim_matrix <- fit_out[[1]]$cosine.similarity
          sig_level <- colnames(cos_sim_matrix)
          signatures <- apply(cos_sim_matrix,1,function(x)names(which.max(x)))
          if(any(grepl("Signature", signatures))){
              signatures <- gsub('Signature ', '', signatures)
              sig_level <- gsub('Signature ', '', sig_level)
          }else if(any(grepl("SBS", signatures))){
              signatures <- gsub('SBS', '', signatures)
              sig_level <- gsub('SBS', '', sig_level)
          }
          treeData$Signature <- signatures[treeData$label]
          treeData[treeData$Signature == ''| is.na(treeData$Signature)]$Signature <- "Unknown"
          if("Unknown" %in% treeData$Signature){
              sig_level <- c(sig_level,"Unknown")
          }
          treeData$Signature <- factor(treeData$Signature, levels = sig_level)
          treeData <- treeData[order(treeData$Signature), ]
      }else{
          branch_type <- getBranchType(phyloTree)
          types <- as.character(branch_type$Mutation_Type)
          names(types) <- as.character(branch_type$Branch_ID)
          treeData$Mutation_Type <- types[treeData$label]
              # dplyr::mutate(Mutation_Type = types[label]) %>% 
              # as.data.table()
      }
   }
   treeData <- treeData[treeData$distance!= 0|treeData$sample == rootLabel,]
   return(treeData)
}


setPhyloTree <- function(tree, treeEdge, treeData, rootNode, mainTrunk, adjacentPoints, adjacentWs, adjacentAngles,
                         horizon = pi/2 , W = pi){
   trunkPath <- rev(c(mainTrunk,rootNode))
   if(length(trunkPath) > 0){
     # subdat_list <- lapply(2:length(trunkPath), function(i){
     #   x1 <- treeData[treeData$end_num == trunkPath[i-1],]$x2
     #   y1 <- treeData[treeData$end_num == trunkPath[i-1],]$y2
     #   W <- W
     #   angle <- horizon
     #   distance <- treeEdge[treeEdge$endNum == trunkPath[i],]$length
     #   if(horizon == pi/2){
     #     x2 <- x1
     #     y2 <- y1 + distance
     #   }else{
     #     x2 <- x1 + distance*cos(angle)
     #     y2 <- y1 + distance*sin(angle)
     #   }
     #   subdat <- data.table::data.table('x1'= x1, 'y1' = y1,
     #                                    'x2' = x2, 'y2' = y2,
     #                                    'horizon' = horizon,'W' = W,
     #                                    'distance' = distance, 'angle' = horizon,
     #                                    'node' = trunkPath[i-1],'end_num' = trunkPath[i])
     #   return(subdat)
     # })
     # 
     # subdat <- dplyr::bind_rows(subdat_list)
     # treeData <- rbind(treeData, subdat)
     
      for(i in 2:length(trunkPath)){
         x1 <- treeData[treeData$end_num == trunkPath[i-1],]$x2
         y1 <- treeData[treeData$end_num == trunkPath[i-1],]$y2
         W <- W
         angle <- horizon
         distance <- treeEdge[treeEdge$endNum == trunkPath[i],]$length
         if(horizon == pi/2){
            x2 <- x1
            y2 <- y1 + distance
         }else{
            x2 <- x1 + distance*cos(angle)
            y2 <- y1 + distance*sin(angle)
         }
         subdat <- data.table::data.table('x1'= x1, 'y1' = y1,
                                          'x2' = x2, 'y2' = y2,
                                          'horizon' = horizon,'W' = W,
                                          'distance' = distance, 'angle' = horizon,
                                          'node' = trunkPath[i-1],'end_num' = trunkPath[i])
         treeData <- rbind(treeData, subdat)
      }
   }
   if(length(adjacentPoints) > 0){
      # subdat_list <- lapply(seq_len(length(adjacentPoints)), function(i){
      #   point <- adjacentPoints[i]
      #   startnode <- treeEdge[treeEdge$endNum == point, ]$node
      #   x1 <- treeData[treeData$end_num == startnode,]$x2
      #   y1 <- treeData[treeData$end_num == startnode,]$y2
      #   W <- adjacentWs[point]
      #   angle <- adjacentAngles[point]
      #   distance <- treeEdge[treeEdge$endNum == point, ]$length
      #   x2 <- x1 + distance*cos(angle)
      #   y2 <- y1 + distance*sin(angle)
      #   subdat <- data.table::data.table('x1'= x1, 'y1' = y1,
      #                                    'x2' = x2, 'y2' = y2,
      #                                    'horizon' = horizon,'W' = W,
      #                                    'distance' = distance, 'angle' = angle,
      #                                    'node' = startnode,'end_num' = point)
      #   return(subdat)
      # })
      # subdat <- dplyr::bind_rows(subdat_list)
      # treeData <- rbind(treeData, subdat)
      for(i in seq_len(length(adjacentPoints))){
         point <- adjacentPoints[i]
         startnode <- treeEdge[treeEdge$endNum == point, ]$node
         x1 <- treeData[treeData$end_num == startnode,]$x2
         y1 <- treeData[treeData$end_num == startnode,]$y2
         W <- adjacentWs[point]
         angle <- adjacentAngles[point]
         distance <- treeEdge[treeEdge$endNum == point, ]$length
         x2 <- x1 + distance*cos(angle)
         y2 <- y1 + distance*sin(angle)
         subdat <- data.table::data.table('x1'= x1, 'y1' = y1,
                                          'x2' = x2, 'y2' = y2,
                                          'horizon' = horizon,'W' = W,
                                          'distance' = distance, 'angle' = angle,
                                          'node' = startnode,'end_num' = point)
         treeData <- rbind(treeData, subdat)
      }
   }
   return(treeData)
}

calMainTrunk <- function(tree, treeEdge, rootNode, rootLabel = "NORMAL"){
   Ntips <- Ntip(tree)
   subtips <- treeEdge[!treeEdge$endNum > Ntips,]$endNum
   distanceTable <- c(0)
   # for(i in seq_len(length(subtips))){
   #    distanceTable <- cbind(distanceTable, -1)
   # }

   s <- vapply(seq_len(length(subtips)), function(i)-1, FUN.VALUE = numeric(1))
   distanceTable <- c(distanceTable, s)
   names(distanceTable) <- c("x",s)
   distanceTable <- distanceTable[(-1)]
   names(distanceTable) <- tree$tip.label[subtips]
   t <- 1
   path <- list()
   while(t <= length(distanceTable)){
      subPath <- c()
      end <- which(tree$tip.label == names(distanceTable)[t])
      name <- names(distanceTable)[t]
      dis <- 0
      while (TRUE) {
         p <- which(treeEdge$endNum == end)
         dis <- dis + treeEdge$length[p]
         start <- treeEdge$node[p]
         subPath <- append(subPath,end)
         end <- start
         if(start == rootNode){
            break
         }
      }
      path[[name]] <- subPath
      distanceTable[name] <- dis
      t <- t + 1
   }
   if(rootLabel == "NORMAL"){
      distanceTable <- distanceTable[names(distanceTable) != rootLabel]
   }
   result <- path[[names(which.max(distanceTable))]]
}
getNodeAngle <- function(tree, treeEdge, mainTrunk,
                           numList, pointsList, rootNode,
                           horizon = pi/2, W = pi ,numNORMAL = NULL){
    # print(pointsList)
   if(!is.null(numNORMAL)){
      treeEdge <- treeEdge[treeEdge$endNum != numNORMAL,]    
   }
   nodes <- factor(rev(c(mainTrunk,rootNode)), levels = rev(c(mainTrunk,rootNode)))
   adjacentPoints <- c()
   adjacentPoints <- lapply(nodes, function(n){
     p <- treeEdge[treeEdge$node == n,]$endNum %>% setdiff(nodes)
     return(p)
   }) %>% unlist()
   # for(n in nodes){
   #    p <- treeEdge[treeEdge$node == n,]$endNum %>% setdiff(nodes)
   #    adjacentPoints <- append(adjacentPoints,p)
   # }
   left <- c()
   right <- c()
   leftList <- c()
   rightList <- c()
   i <- 1
   
   # adresult <- lapply(seq_len(length(adjacentPoints)),function(i){
   #   point <- adjacentPoints[i]
   #   pos <- which(pointsList == point)
   #   if(i == 1){
   #     if(horizon < pi/2){
   #       return(list(l = numList[pos], ll = point, r = NULL, rl = NULL))
   #     }else{
   #       return(list(l = NULL, ll = NULL, r = numList[pos], rl = point))
   #     }
   #   }else{
   #     if(sum(right) <= sum(left)){
   #       return(list(l = NULL, ll = NULL, r = numList[pos], rl = point))
   #     }else{
   #       return(list(l = numList[pos], ll = point, r = NULL, rl = NULL))
   #     }
   #   }
   # })
   # left <- lapply(adresult, function(x)x$l) %>% unlist()
   # leftList <- lapply(adresult, function(x)x$ll) %>% unlist()
   # right <- lapply(adresult, function(x)x$r) %>% unlist()
   # rightList <- lapply(adresult, function(x)x$rl) %>% unlist()
   for(point in adjacentPoints){
      if(i == 1){
         if(horizon < pi/2){
            leftList <- append(leftList, point)
            pos <- which(pointsList == point)
            left <- append(left,numList[pos])
         }else{
            rightList <- append(rightList, point)
            pos <- which(pointsList == point)
            right <- append(right,numList[pos])
         }
      }else{
         if(sum(right) <= sum(left)){
            rightList <- append(rightList, point)
            pos <- which(pointsList == point)
            right <- append(right,numList[pos])
         }else{
            leftList <- append(leftList, point)
            pos <- which(pointsList == point)
            left <- append(left,numList[pos])
         }
      }
      i <- i+1
   }
   adjacentWs <- rep(0,length(tree$edge.length)+1)
   adjacentAngles <- rep(0,length(tree$edge.length)+1)
   
   if(length(rightList) > 0){
      startr <- horizon - W/2
      wrt <- W/2
      totalR <- sum(right)
      
      for(i in seq_len(length(rightList))){
         if(rightList[i] <=  length(tree$tip.label)){
            n <- 1
         }else{
            n <- right[i]
         }
         wr <- wrt*(n/totalR)
         adjacentWs[rightList[i]] <- wr
         if(i == 1){
            angler <- wr/2 + startr
            c1 <- horizon != pi/2
            c2 <- any(startr != 0,n == 1)
            c3 <- any(startr != pi/2,n == 1)
            if(all(c1,c2,c3)){
                if(length(c(leftList,rightList))!=1){
                    angler <- startr
                }
            }
            adjacentAngles[rightList[i]] <- angler
            startr <- angler
         }else{
            
            if(n == 1 & horizon != pi/2){
               startr <- startr - wr/2
            }
            
            angler <- wr/2 + startr + adjacentWs[rightList[i-1]]/2
            adjacentAngles[rightList[i]] <- angler
            startr <- angler
         }
      }
   }
   if(length(leftList) > 0){
      wlt <- W/2
      startl <- horizon + wlt
      totalL <- sum(left)
      for(i in seq_len(length(leftList))){
         if(leftList[i] <=  length(tree$tip.label)){
            n <- 1
         }else{
            n <- left[i]
         }
         wl <- wlt*(n/totalL)
         adjacentWs[leftList[i]] <- wl
         if(i == 1){
            anglel <- startl - wl/2
            c1 <- horizon != pi/2
            c2 <- any(startl != pi, n == 1)
            c3 <- any(startl != pi/2, n == 1)
            # angle extension
            if(all(c1,c2,c3)){
                if(length(c(leftList,rightList))!=1){
                  anglel <- startl  
                }
               
            }
            adjacentAngles[leftList[i]] <- anglel
            startl <- anglel
         }else{
            
            if(n == 1& horizon != pi/2){
               startl <- startl  + wl/2
            }
            
            anglel <- startl- wl/2 - adjacentWs[leftList[i-1]]/2
            adjacentAngles[leftList[i]] <- anglel
            startl <- anglel
         }
      }
   }
   return(list(adjacentWs,adjacentAngles, adjacentPoints))
}


calChildNodeNum <- function(tree, treeEdge, mainTrunk, rootNode, ft = FALSE){
   pointsList <- sort(unique(c(treeEdge$node,treeEdge$endNum))) 
   numList <- rep(1, length(pointsList))
   lrnumList <- rep(0, length(pointsList))
   # for(i in seq_len(length(pointsList))){
   #    numList[i] <- 1
   #    lrnumList[i] <- 0
   # }
   for(i in seq_len(length(pointsList))){
      if(pointsList[i] == rootNode){
         next
      }
      start <- pointsList[i]
      end <- as.numeric(treeEdge[treeEdge$endNum == start,]$node)
      while(TRUE){
         pos <- which(pointsList == end)
         numList[pos] <- numList[pos] + 1
         if(end == rootNode){
            break
         }
         start <- treeEdge[treeEdge$endNum == end, ]$endNum
         end <- treeEdge[treeEdge$endNum == end, ]$node
      }
      # if(count > 1){
      #     numList[i] <- 0
      # }
   }
   ## samples are not in mainTrunk
   nodeRange <- c(Ntip(tree)+1):(length(tree$edge.length)+1)
   nodeOnTree <- c(mainTrunk,rootNode) %>% intersect(nodeRange)
   nodeNoOnTree <- treeEdge[treeEdge$node %in% c(mainTrunk,rootNode), ]$endNum %>% setdiff(nodeOnTree) %>% intersect(nodeRange)
   return(list(numList,pointsList, nodeNoOnTree, nodeOnTree))
}

getSigColors <- function(signatures){
    
   signature_colors <- c("1" = "#E41A1C", "2" = "#377EB8","3" = "#7F0000",
                        "4" = "#35978f","5" = "#A6D854", "6" = "#FFD92F",
                        "7" = "#E5C494", "8" = "#8DD3C7","9" = "#6E016B",
                        "10" = "#BEBADA", "11" = "#e08214", "12" =  "#80B1D3",
                        "13" = "#d6604d", "14" = "#ffff99","15" = "#FCCDE5",
                        "16" = "#FF6A5A","17" = "#BC80BD","18" = "#CCEBC5",
                        "19" = "#fb9a99", "20"  = "#B6646A", "21" = "#9F994E", 
                        "22" = "#7570B3" , "23" = "#c51b7d" , "24" = "#66A61E" ,
                         "25"  = "#E6AB02" ,"26" = "#003c30",  "27" =  "#666666",
                         "28" = "#524762","29" =  "#8926AC","30" = "#C42669",
                        "1A" = "#47B39F","1B" = "#B84C60", "R1" = "#B93183",
                        "R2" = "#B93183",  "R3" = "#845AD1", "U1" = "#7BA52E",
                        "U2" = "#7A4287","7a" = "#663399", "7b" = "#AD8B5E",
                        "7c" = "#194682", "7d" = "#FF8931",  "10a" = "#AD3D07",
                        "10b" = "#0D5DFF", "17a" = "#F2A200", "17b" = "#8C3B45",
                         "31" = "#009F9A", "32" = "#FF2B32", "33" = "#153E71",
                        "34" = "#FF8C2F", "35" ="#2F8CFF", "36" = "#62BA1F",
                        "37" = "#9745D6",  "38" = "#3826CA", "39" = "#96C4BC",
                        "40"  = "#982B9F", "41" = "#9C3C33", "42" = "#3EB15A",
                        "43" = "#EF1297",  "44" = "#B1098C",  "45" = "#883E93",
                        "46" = "#6ED66B", "47" ="#9F679F" , "48" = "#5B309F",
                        "49" = "#4E7A9F", "50" = "#9F6878",  "51" =  "#154E0F",
                        "52" =  "#9D887D", "53" = "#143045", "54" =  "#8F1B1D",
                        "55" = "#544C9F", "56"  =  "#7F4E10", "57" =  "#429F74",
                        "58" = "#6E6A5F", "59" = "#235F9F", "60" =  "#8A1D1B",
                        "84" = "#4B9F34", "85" = "#598F4D", "Unknown" = "black")
   color_scale <- signature_colors[signatures]
   return(color_scale)
}

labelBranch <- function(tree){
   Root <- which(tree$tip.label == "NORMAL")
   internalNodes <- sort(unique(tree$edge[,1]))
   result <- list()
   subnumList <- list()
   subnumList <- lapply(seq_len(length(tree$edge.length)+1), function(i){
     if(i > Ntip(tree)){
       return(i)
     }
     else{
       return(NA)
     }
   })
   result <- rep(NA, (length(tree$edge.length)+1)) %>% as.list()
   # for(i in seq_len(length(tree$edge.length)+1)){
   #    result[[i]] <- NA
   #    if(i > Ntip(tree)){
   #       subnumList[[i]] <- i
   #    }
   #    else{
   #       subnumList[[i]] <- NA
   #    }
   # }
   # print(result)
   # print(subnumList)
   end <- tree$edge[which(tree$edge[,2] == Root),1]
   for(i in seq_len(length(tree$tip.label))){
      row <- tree$edge[which(tree$edge[,2] == i), ]
      if(i == Root){
         next
      }
      while (TRUE) {
         node <- row[1]
         if(node == end){
            result[[node]] <- append(result[[node]], i)
            break
         }
         else{
            result[[node]] <- append(result[[node]], i)
            row <- tree$edge[which(tree$edge[,2] == node),]
         }
      }
   }
   rootNode <- tree$edge[which(tree$edge[,2] == Root),1]
   for(i in (length(tree$tip.label)+1):(length(tree$edge.length)+1)){
      if(i == rootNode){
         next
      }
      row <- tree$edge[which(tree$edge[,2] == i), ]
      while (TRUE) {
         node <- row[1]
         if(node == rootNode){
            subnumList[[node]] <- append(subnumList[[node]], i)
            break
         }
         else{
            subnumList[[node]] <- append(subnumList[[node]], i)
            row <- tree$edge[which(tree$edge[,2] == node),]
         }
      }
   }
   g <- lapply(result, function(x){return(paste(sort(tree$tip.label[x[-1]], decreasing = TRUE),collapse = "&"))})
   return(list(g,subnumList))
}
