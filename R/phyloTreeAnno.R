getTreeData <- function(phyloTree = NULL,
                        branchCol = "mutSig",
                        compare = FALSE,
                        ...){
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
   for(i in 1:nrow(treeData)){
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
   for(i in 1:nrow(treeData)){
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
          cos_sim_matrix <- fitSignatures(tri_matrix,...)[[1]]$cosine.similarity
          signatures <- apply(cos_sim_matrix,1,function(x)names(which.max(x)))
          if(any(grepl("Signature", signatures))){
              signatures <- gsub('Signature ', '', signatures)
          }else if(any(grepl("SBS", signatures))){
              signatures <- gsub('SBS', '', signatures)
          }
          treeData <- treeData[, Signature:= signatures[label]]
          # print(treeData$label)
          # print(signatures)
          treeData[Signature == ''|is.na(Signature)]$Signature <- "Unknown"
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
      for(i in 1:length(adjacentPoints)){
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
   for(i in 1:length(subtips)){
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
      
      for(i in 1:length(rightList)){
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
      for(i in 1:length(leftList)){
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
   for(i in 1:length(pointsList)){
      numList[i] <- 1
      lrnumList[i] <- 0
   }
   for(i in 1:length(pointsList)){
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
   for(i in 1:(length(tree$edge.length)+1)){
      result[[i]] <- NA
      if(i > Ntip(tree)){
         subnumList[[i]] <- i 
      }
      else{
         subnumList[[i]] <- NA
      }
   }
   end <- tree$edge[which(tree$edge[,2] == Root),1]
   for(i in 1:length(tree$tip.label)){
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
   g <- lapply(result, function(x){return(paste(sort(tree$tip.label[x[-1]], decreasing = T),collapse = "∩"))})
   return(list(g,subnumList))
}
