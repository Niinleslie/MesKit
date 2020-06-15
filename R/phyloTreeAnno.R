plotTree <- function(phyloTree,
                     treeData = NULL,
                     branchCol = "mutType",
                     show.bootstrap = TRUE,
                     signaturesRef = "cosmic_v2",
                     min.mut.count = 15,
                     min.ratio = 1/20,
                     common.col = "red"){
    
    if(min.ratio <= 0|min.ratio > 1){
        stop("Error: min.ratio should be within (0,1]")
    }
    
    ## adjust length of branches by min.ratio
    tree <- getTree(phyloTree)
    min.len <- max(tree$edge.length)*min.ratio
    tree$edge.length[tree$edge.length < min.len] <- min.len
    phyloTree@tree <- tree
    
    patient <- getPhyloTreePatient(phyloTree)
    if(is.null(treeData)){
        treeData <- getTreeData(phyloTree = phyloTree,
                                branchCol = branchCol,
                                signaturesRef = signaturesRef,
                                min.mut.count = min.mut.count)
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
    samplesLength <- nrow(treeData[sample != "internal node",]) 
    if(samplesLength > 7){
        samplePointsSize <- 1.5 
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
    rootNode <- treeData[sample == rootLabel,]$node
    if(length(boot_value) == 1){
        bootsData <- data.frame(x2 = 0, y2 = 0, node = rootNode, end_num = rootNode, boots = boot_value)
    }else{
        sub <- data.table::data.table(x2 = 0, y2 = 0,node = rootNode, end_num = rootNode)
        bootsData <- rbind(treeData[sample == 'internal node',][,.(x2,y2,node,end_num)],sub)
        boots <- c()
        LN <- min(bootsData$node)-1
        for(i in seq_len(nrow(bootsData))){      
            if(i == nrow(bootsData)){
                boots <- append(boots, boot_value[rootNode - LN])
                next
            }
            boots <- append(boots, boot_value[bootsData$end_num[i] - LN])
            
        }
        bootsData <- cbind(bootsData, boots = boots)
    }
    ## get the max value of X axis 
    x_max <- max(abs(treeData$x2))
    p <- ggplot(data = treeData) + 
        ## balance the space on the left and right sides
        geom_segment(aes(x = 0, y = 0, xend = x_max, yend = 0),color = "white",size = 0.01)+
        geom_segment(aes(x = 0, y = 0, xend = -x_max, yend = 0),color = "white",size = 0.01)
    
    textAdjust <- mean(as.numeric(treeData$distance))
    
    if(!is.null(branchCol)){
        if(branchCol == "mutSig"){
            color_scale <- getSigColors(as.character(unique(treeData$Signature)) )
            sig_level <- levels(treeData$Signature)
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Signature), size=segmentSize)
            p <- p + scale_color_manual(breaks = sig_level ,values = color_scale) + 
                theme(legend.title = element_text(size = legend.title.size))
        }
        else{
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
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = factor(Mutation_Type)), size=segmentSize)
            ## remove legend title
            p <- p +  theme(legend.title = element_blank()) + 
                scale_color_manual(breaks = type.level, values = type_all_colors)
        }
    }
    else{
        if(compare){
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                                  color = common.col,
                                  data = treeData[is.match != "NO"],
                                  size = segmentSize)
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2),
                                  color = "black",
                                  data = treeData[is.match == "NO"],size = segmentSize)
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
                   panel.grid.minor = element_blank(), panel.background = element_blank(),
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
                             data = treeData[(!sample %in% c("internal node",rootLabel)) &
                                                 x2 >= 0, ],
                             size = sampleTextSize ,force = 10)
    p <- p + geom_text_repel(aes(x = x2, y = y2, label = sample),
                             nudge_y = textAdjust/10,
                             nudge_x = -textAdjust/10,
                             segment.color = "grey",
                             segment.size = 0.25,
                             data = treeData[(!sample %in% c("internal node",rootLabel)) &
                                                 x2 < 0, ],
                             size = sampleTextSize ,force = 10)
    ## label NORMAL
    p <- p + geom_text(aes(x = x2,y = y2-textAdjust/5),
                       label = rootLabel,
                       data = treeData[sample == rootLabel,], 
                       size = sampleTextSize)
    
    p <- p + geom_point(aes(x = x2,y = y2),
                        data = treeData[sample == 'internal node',],
                        size = nodePointsSize, color = "#8c510a", fill = "white", shape = 21,
                        stroke = nodeStrokeSize)
    p <- p + geom_point(aes(x = x2, y = y2), 
                        data = treeData[sample != 'internal node',],
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
                                  force = 5)
    }
    if(compare){
        p <- p + geom_label_repel(aes(x = x1 + (x2-x1)/2 , y = y1 + (y2 - y1)/2,label = is.match),
                                  data = treeData[is.match != "NO",],
                                  fontface = 'bold', 
                                  size = bootLabelSize, box.padding = unit(bootPaddingSize, "lines"), point.padding = unit(0.5, "lines"),
                                  segment.colour = "grey50", segment.size = 0.5, force = 5)
    }
    
    if(compare){
        tree.title <- patient
    }else{
        tree.title <- paste(patient," (n=" ,nrow(getBinaryMatrix(phyloTree)) ,")",sep = "")
    }
    p <- p + 
        ggtitle(tree.title)+
        theme(plot.title = element_text(face = "bold",colour = "black", hjust = 0.5,size = 13.5))
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
         subEdge <- treeEdge[node %in% subnodes, ]
         rootNode <- node
         numNORMAL <- NULL
         mainTrunk <- calMainTrunk(tree, subEdge, rootNode = rootNode, rootLabel = "")
         ccn <- calChildNodeNum(tree, subEdge, mainTrunk, rootNode = rootNode, ft = TRUE)
         numList <- ccn[[1]]
         pointsList <- ccn[[2]]
         nodeNoOnTree <- append(nodeNoOnTree,ccn[[3]])
         nodeOnTree <- append(nodeOnTree,ccn[[4]])
         horizon <- treeData[end_num == rootNode,]$angle
         W <- treeData[end_num == rootNode,]$W
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
   treeData[,sample := ""]
   for(i in seq_len(nrow(treeData))){
      if(treeData$end_num[i] > length(tree$tip.label)){
         treeData$sample[i] <- "internal node"
      }
      else{
         pos <- treeData$end_num[i]
         treeData$sample[i] <- tree$tip.label[pos]
      }
   }
   if(nrow(treeData) == 3){
      angleList <- c(pi/3, 2*pi/3)
      rows <- which(treeData$sample != "NORMAL")
      i = 1
      for(row in rows){
         angle <- angleList[i]
         treeData$angle[row] <- angle
         treeData$x2[row] <- treeData$distance[row]*cos(angle)
         treeData$y2[row] <- treeData$distance[row]*sin(angle)
         i = i+1
      }
   }
   if(nrow(treeData[x2 > 0,]) == 0){
      maxy <- which.max(treeData$y2)
      angle <- treeData$angle[maxy] - pi/18
      treeData$angle[maxy] <- angle
      treeData$y2[maxy] <- treeData$distance[maxy]*sin(angle) + treeData$y1[maxy]
      treeData$x2[maxy] <- treeData$distance[maxy]*cos(angle) + treeData$x1[maxy]
   }else if(nrow(treeData[x2 < 0,]) == 0){
      maxy <- which.max(treeData$y2)
      angle <- treeData$angle[maxy] + pi/18
      treeData$angle[maxy] <- angle
      treeData$y2[maxy] <- treeData$distance[maxy]*sin(angle) + treeData$y1[maxy]
      treeData$x2[maxy] <- treeData$distance[maxy]*cos(angle) + treeData$x1[maxy]
   }
   
   ## label represents the common evolution path of samples
   treeData$label <- ""
   for(i in seq_len(nrow(treeData))){
       if(treeData$sample[i] == "NORMAL"){
           treeData$label[i] <- branchLabel[[rootNode]]
       }else{
           if(treeData$end_num[i] > length(tree$tip.label)){
               treeData$label[i] <- branchLabel[[treeData$end_num[i]]]
           }
           else{
               treeData$label[i] <- as.character(treeData$sample[i]) 
           }
       }
   }
   if(!is.null(branchCol) & !compare){
       ## add signature
      if(branchCol == "mutSig"){
          tri_matrix <- triMatrix(phyloTree,withinTumor = FALSE)
          fit_out <- fitSignatures(tri_matrix,
                                          signaturesRef = signaturesRef,
                                          min.mut.count = min.mut.count)
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
          treeData <- treeData[, Signature:= signatures[label]]
          treeData[Signature == ''|is.na(Signature)]$Signature <- "Unknown"
          if("Unknown" %in% treeData$Signature){
              sig_level <- c(sig_level,"Unknown")
          }
          treeData$Signature <- factor(treeData$Signature, levels = sig_level)
          treeData <- treeData[order(Signature), ]
      }else{
          branch_type <- getBranchType(phyloTree)
          types <- as.character(branch_type$Mutation_Type)
          names(types) <- as.character(branch_type$Branch_ID)
          treeData <- treeData[, Mutation_Type:= types[label]] 
              # dplyr::mutate(Mutation_Type = types[label]) %>% 
              # as.data.table()
      }
   }
   treeData <- treeData[distance!= 0|sample == rootLabel]
   return(treeData)
}


setPhyloTree <- function(tree, treeEdge, treeData, rootNode, mainTrunk, adjacentPoints, adjacentWs, adjacentAngles,
                         horizon = pi/2 , W = pi){
   trunkPath <- rev(c(mainTrunk,rootNode))
   if(length(trunkPath) > 0){
      for(i in 2:length(trunkPath)){
         x1 <- treeData[end_num == trunkPath[i-1],]$x2
         y1 <- treeData[end_num == trunkPath[i-1],]$y2
         W <- W
         angle <- horizon
         distance <- treeEdge[endNum == trunkPath[i],]$length
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
      for(i in seq_len(length(adjacentPoints))){
         point <- adjacentPoints[i]
         startnode <- treeEdge[endNum == point, ]$node
         x1 <- treeData[end_num == startnode,]$x2
         y1 <- treeData[end_num == startnode,]$y2
         W <- adjacentWs[point]
         angle <- adjacentAngles[point]
         distance <- treeEdge[endNum == point, ]$length
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
   subtips <- treeEdge[!endNum > Ntips,]$endNum
   distanceTable <- data.frame(x = 0)
   for(i in seq_len(length(subtips))){
      distanceTable <- cbind(distanceTable, -1)
   }
   distanceTable <- distanceTable[,(-1)]
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
      treeEdge <- treeEdge[endNum != numNORMAL,]    
   }
   nodes <- factor(rev(c(mainTrunk,rootNode)), levels = rev(c(mainTrunk,rootNode)))
   adjacentPoints <- c()
   for(n in nodes){
      p <- treeEdge[node == n,]$endNum %>% setdiff(nodes)
      adjacentPoints <- append(adjacentPoints,p)
   }
   left <- c()
   right <- c()
   leftList <- c()
   rightList <- c()
   i <- 1
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
         wr <- wrt*n/totalR
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
         }
         else{
            n <- left[i]
         }
         wl <- wlt*n/totalL
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
         }
         else{
            
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
   numList <- c()
   lrnumList <- c()
   for(i in seq_len(length(pointsList))){
      numList[i] <- 1
      lrnumList[i] <- 0
   }
   for(i in seq_len(length(pointsList))){
      if(pointsList[i] == rootNode){
         next
      }
      start <- pointsList[i]
      end <- as.numeric(treeEdge[endNum == start,]$node)
      while(TRUE){
         pos <- which(pointsList == end)
         numList[pos] <- numList[pos] + 1
         if(end == rootNode){
            break
         }
         start <- treeEdge[endNum == end, ]$endNum
         end <- treeEdge[endNum == end, ]$node
      }
      # if(count > 1){
      #     numList[i] <- 0
      # }
   }
   ## samples are not in mainTrunk
   nodeRange <- c(Ntip(tree)+1):(length(tree$edge.length)+1)
   nodeOnTree <- c(mainTrunk,rootNode) %>% intersect(nodeRange)
   nodeNoOnTree <- treeEdge[node %in% c(mainTrunk,rootNode), ]$endNum %>% setdiff(nodeOnTree) %>% intersect(nodeRange)
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
   for(i in seq_len(length(tree$edge.length)+1)){
      result[[i]] <- NA
      if(i > Ntip(tree)){
         subnumList[[i]] <- i 
      }
      else{
         subnumList[[i]] <- NA
      }
   }
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
   g <- lapply(result, function(x){return(paste(sort(tree$tip.label[x[-1]], decreasing = TRUE),collapse = "∩"))})
   return(list(g,subnumList))
}
