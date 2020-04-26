getTreeData <- function(phyloTree = NULL,
                        branchCol = "mutSig",
                        min.mut.count = 15){
   tree <- phyloTree@tree
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
      gta <- getTrunkAngles(tree = tree, treeEdge = treeEdge,
                            mainTrunk = mainTrunk, numList = numList,
                            pointsList = pointsList, rootNode = rootNode,
                            horizon = pi/2, W = pi, numNORMAL = numNORMAL)
      collateralWs <- gta[[1]]
      collateralAngles <- gta[[2]]
      collateralPoints <- gta[[3]]
   }
   #X1, y1 are the starting point, x2 and y2 are the end point.Horizon stands for the Angle between the line and the positive x axis;W stands for the Angle occupied by branches (explained in the paper)
   #Distance means distance;Node is the internal node with the start;End_num is the number of each node in edge;
   treeData <- data.table::data.table('x1'= 0, 'y1' = 0,
                                      'x2' = 0, 'y2' = 0,
                                      'horizon' = pi/2, 'W' = pi,
                                      'angle' = 0,'distance' = 0,
                                      'node' = rootNode, 'end_num' = rootNode)
   treeData <- setPhyloTree(tree = tree, treeEdge = treeEdge, treeData = treeData, rootNode = rootNode,
                            mainTrunk = mainTrunk, collateralPoints = collateralPoints,
                            collateralWs = collateralWs, collateralAngles = collateralAngles,
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
         gta <- getTrunkAngles(tree = tree, treeEdge = subEdge,
                               mainTrunk = mainTrunk, numList = numList,
                               pointsList = pointsList, rootNode = rootNode,
                               horizon = horizon, W = W, numNORMAL = NULL)
         collateralWs <- gta[[1]]
         collateralAngles <- gta[[2]]
         collateralPoints <- gta[[3]]
         treeData <- setPhyloTree(tree = tree, treeEdge = subEdge, treeData = treeData, rootNode = rootNode,
                                  mainTrunk = mainTrunk, collateralPoints = collateralPoints,
                                  collateralWs = collateralWs, collateralAngles = collateralAngles,
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
   if(!is.null(branchCol)){
      signature <- doTreeMutSig(phyloTree, min.mut.count = min.mut.count)$mutSigsOutput %>%
         dplyr::group_by(branch) %>% 
         dplyr::filter(sig.prob == max(sig.prob)) %>% 
         dplyr::ungroup() %>% 
         as.data.frame()
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
      treeData <- addSignature(tree, treeData, signature)
   }
   treeData <- treeData[distance!= 0|sample == rootLabel,]
   return(treeData)
}


setPhyloTree <- function(tree, treeEdge, treeData, rootNode, mainTrunk, collateralPoints, collateralWs, collateralAngles,
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
   if(length(collateralPoints) > 0){
      for(i in 1:length(collateralPoints)){
         point <- collateralPoints[i]
         startnode <- treeEdge[endNum == point, ]$node
         x1 <- treeData[end_num == startnode,]$x2
         y1 <- treeData[end_num == startnode,]$y2
         W <- collateralWs[point]
         angle <- collateralAngles[point]
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
getTrunkAngles <- function(tree, treeEdge, mainTrunk,
                           numList, pointsList, rootNode,
                           horizon = pi/2, W = pi ,numNORMAL = NULL){
   if(!is.null(numNORMAL)){
      treeEdge <- treeEdge[endNum != numNORMAL,]    
   }
   nodes <- factor(rev(c(mainTrunk,rootNode)), levels = rev(c(mainTrunk,rootNode)))
   collateralPoints <- c()
   for(n in nodes){
      p <- treeEdge[node == n,]$endNum %>% setdiff(nodes)
      collateralPoints <- append(collateralPoints,p)
   }
   left <- c()
   right <- c()
   leftList <- c()
   rightList <- c()
   i <- 1
   for(point in collateralPoints){
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
   collateralWs <- rep(0,length(tree$edge.length)+1)
   collateralAngles <- rep(0,length(tree$edge.length)+1)
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
         collateralWs[rightList[i]] <- wr
         if(i == 1){
            angler <- wr/2 + startr
            if(all(
               horizon != pi/2,
               (startr != 0|n == 1),
               (startr != pi/2 | n == 1)
            )){
               angler <- startr
            }
            collateralAngles[rightList[i]] <- angler
            startr <- angler
         }else{
            
            if(n == 1 & horizon != pi/2){
               startr <- startr - wr/2
            }
            
            angler <- wr/2 + startr + collateralWs[rightList[i-1]]/2
            collateralAngles[rightList[i]] <- angler
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
         collateralWs[leftList[i]] <- wl
         if(i == 1){
            anglel <- startl - wl/2
            # angle extension
            if(all(
               horizon != pi/2,
               (startl != pi|n == 1),
               (startl != pi/2 | n == 1)
            )){
               anglel <- startl
            }
            
            collateralAngles[leftList[i]] <- anglel
            startl <- anglel
         }
         else{
            
            if(n == 1& horizon != pi/2){
               startl <- startl  + wl/2
            }
            
            anglel <- startl- wl/2 - collateralWs[leftList[i-1]]/2
            collateralAngles[leftList[i]] <- anglel
            startl <- anglel
         }
      }
   }
   return(list(collateralWs,collateralAngles, collateralPoints))
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



addSignature <- function(tree, treeData, signature){
   #add signature to treeData
   treeData$Signature <- ''
   treeData$Branch_Tumor_Type <- ''
   sigs <- strsplit(as.character(signature$branch),"∩")
   sigs <- lapply(sigs, function(x){return(paste(sort(x,decreasing = T),collapse = "∩"))})
   t <- 1
   while(t<=length(sigs)){
      pos <- which(treeData$label == sigs[[t]])
      treeData$Signature[pos] <- as.character(signature$sig[t]) 
      treeData$Branch_Tumor_Type[pos] <- as.character(signature$Branch_Tumor_Type[t])
      t <- t + 1
   }
   if(treeData$Signature[which(treeData$sample == 'NORMAL')] == ''){
      treeData$Signature[which(treeData$sample == 'NORMAL')] = as.character(signature$sig[1])
      treeData$Branch_Tumor_Type[which(treeData$sample == 'NORMAL')] = as.character(signature$Branch_Tumor_Type[1])
      
   }
   treeData[treeData$Signature == '',]$Signature <- "noMapSig"
   treeData <- treeData[order(treeData$Signature), ]
   treeData$Signature <- gsub('Signature ', '', treeData$Signature)
   return(treeData)
}

getColors <- function(signatures){
   allColorScale <- c("#E41A1C","#377EB8","#7F0000",
                      "#35978f","#FC8D62","#2166ac",
                      "#E78AC3","#A6D854","#FFD92F",
                      "#E5C494","#8DD3C7", "#6E016B" ,
                      "#BEBADA", "#e08214", "#80B1D3",
                      "#d6604d","#ffff99","#FCCDE5",
                      "#FF6A5A","#BC80BD","#CCEBC5" ,
                      "#fb9a99","#B6646A", "#9F994E", 
                      "#7570B3" ,"#c51b7d" ,"#66A61E" ,
                      "#E6AB02" ,"#003c30", "#666666",'black')
   allSignature <- append(gsub('Signature.', '',row.names(deconstructSigs::signatures.cosmic)), 'noMapSig')
   colorScale <- c()
   
   for(i in 1:length(signatures)){
      color <- allColorScale[which(allSignature == signatures[i])]
      colorScale <- append(colorScale, color)
   }
   
   if('black' %in% colorScale){
      num <- which(colorScale == 'black')
      colorScale <- colorScale[-num]
      colorScale <- append(colorScale, 'black')
   }
   
   return(colorScale)
}

getPrivateMutation <- function(phyloTree){
   mut.branches <- phyloTree@mut.branches
   private.idx <- unlist(lapply(unique(mut.branches$Branch_ID),
                                 function(x){return(length(strsplit(x,"∩")[[1]]) == 1)}))
   samples <- unique(mut.branches$Branch_ID)[private.idx]
   private.muts <- mut.branches[mut.branches$Branch_ID %in% samples,]
   totalMutSum <- nrow(mut.branches)
   privateMutSum <- nrow(private.muts)
   privateMutProportion <- paste(round((privateMutSum/totalMutSum)*100,1),"%",sep = "")
   return(list(totalMutSum, privateMutProportion))
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


drawPhyloTree <- function(phyloTree = NULL,
                          treeData = NULL,
                          branchCol = TRUE,
                          show.bootstrap = TRUE,
                          use.box = TRUE,
                          show.heatmap = TRUE,
                          use.ccf = FALSE,
                          compare = FALSE,
                          common.lty = "solid",
                          min.ratio = 1/30,
                          show.class.label = FALSE,
                          geneList = NULL,
                          plot.geneList = FALSE,
                          show.gene = FALSE,
                          show.geneList = TRUE,
                          min.mut.count = 15,
                          mut.threshold = 150,
                          use.shiny = FALSE){
    
    patientID <- phyloTree@patientID
    
    ## shiny progression
    if(use.shiny){
        incProgress(amount=1)
        setProgress(message = paste('Drawing ', "phylogenetic tree - ", patientID, sep=""))
    }
    
   if(!is.null(min.ratio)){
      if(min.ratio > 0){
         min.len <- max(phyloTree@tree$edge.length)*min.ratio
         phyloTree@tree$edge.length[phyloTree@tree$edge.length < min.len] <- min.len
      }
   }
   if(is.null(treeData)){
      treeData <- getTreeData(phyloTree = phyloTree,
                              branchCol = branchCol,
                              min.mut.count = min.mut.count)
   }
   set.seed(1234)
   myBoots <- phyloTree@bootstrap.value
   rootLabel <- "NORMAL"
   ## plot phylotree
   samplePointsSize <- 3
   sampleTextSize <- 3
   nodePointsSize <- 1.7
   segmentSize <- 1.5
   # nodeStrokeSize <- 0.5
   # sampleStrokeSize <- 1
   nodeStrokeSize <- 0.25
   sampleStrokeSize <- 1.5
   bootLabelSize <- 2.2
   bootTextSize <- 2.2
   bootPaddingSize <- 0.35
   samplesLength <- nrow(treeData[sample != "internal node",]) 
   if(samplesLength > 7){
      samplePointsSize <- 1.5 
      sampleTextSize <- 2
      segmentSize <- 0.8
      nodePointsSize <- 0.8
      # nodeStrokeSize <- 0.25
      # sampleStrokeSize <- 0.5
      nodeStrokeSize <- 0.15
      sampleStrokeSize <- 0.8
      bootLabelSize <- 1.5
      bootTextSize <- 1.5
      bootPaddingSize <- 0.1
   }
   rootNode <- treeData[sample == rootLabel,]$node
   if(length(myBoots) == 1){
      bootsData <- data.frame(x2 = 0, y2 = 0, node = rootNode, end_num = rootNode, boots = myBoots)
   }else{
      sub <- data.table::data.table(x2 = 0, y2 = 0,node = rootNode, end_num = rootNode)
      bootsData <- rbind(treeData[sample == 'internal node',][,.(x2,y2,node,end_num)],sub)
      boots <- c()
      LN <- min(bootsData$node)-1
      for(i in 1:nrow(bootsData)){      
         if(i == nrow(bootsData)){
            boots <- append(boots,  myBoots[rootNode - LN])
            next
         }
         boots <- append(boots, myBoots[bootsData$end_num[i] - LN])
         
      }
      bootsData <- cbind(bootsData, boots = boots)
   }
   p <- ggplot(data = treeData)
   textAdjust <- mean(as.numeric(treeData$distance))
   # maxy <- max(abs(treeData$y2))
   # maxx <- max(abs(treeData$x2)))
   # ratecoord <- maxy/maxx
   if(!is.null(branchCol)){
      if(branchCol == "mutSig"){
         colorScale <- getColors(unique(treeData$Signature))
         if(compare){
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Signature),
                                  data = treeData[is.match != "NO"],
                                  linetype = common.lty,
                                  size = segmentSize)
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Signature),
                                  data = treeData[is.match == "NO"],size = segmentSize)
         }
         else{
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Signature), size=segmentSize)
         }
         p <- p + scale_color_manual(values = colorScale)
      }
      else{
         
         ## sort branch tumor type 
         all.types <- unique(treeData$Branch_Tumor_Type) 
         public <- all.types[grep("Public", all.types)] 
         shared <- all.types[grep("Shared", all.types)] 
         private <- all.types[grep("Private", all.types)]
         type.level <- c(public, shared, private)
         treeData$Branch_Tumor_Type <- factor(treeData$Branch_Tumor_Type, levels = type.level)
         
         if(compare){
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Branch_Tumor_Type),
                                  data = treeData[is.match != "NO"],
                                  linetype = common.lty,
                                  size = segmentSize)
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Branch_Tumor_Type),
                                  data = treeData[is.match == "NO"],size = segmentSize)
         }
         else{
            p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = Branch_Tumor_Type), size=segmentSize)
         }
         p <- p + scale_color_discrete(breaks = type.level)
      }
   }
   else{
      if(compare){
         p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = sample),
                               data = treeData[!sample %in% c("internal node",rootLabel),], 
                               size = segmentSize, show.legend = T)
         p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = 'black',
                               data = treeData[sample == rootLabel,], 
                               size = segmentSize, show.legend = F )
         p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = '#67001F',
                               data = treeData[sample == 'internal node'&is.match != "NO",],
                               linetype = common.lty,
                               size = segmentSize, show.legend = F )
         p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = '#67001F',
                               data = treeData[sample == 'internal node'&is.match == "NO",],
                               size = segmentSize, show.legend = F)
      }else{
         p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = sample),
                               data = treeData[!sample %in% c("internal node",rootLabel),], 
                               size=segmentSize, show.legend = T)
         p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = 'black',
                               data = treeData[sample == rootLabel,], 
                               size = segmentSize, show.legend = F )
         p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = '#67001F',
                               data = treeData[sample == 'internal node',], 
                               size = segmentSize, show.legend = F )
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
                  legend.title = element_text(),
                  legend.position = 'top',
                  legend.direction = "horizontal") + 
      guides(color = guide_legend(nrow=1))+
      coord_fixed(ratio= 1) +
      scale_x_discrete(expand = expansion(add = mean(treeData$distance)))+
      scale_y_discrete(expand = expansion(add = mean(treeData$distance)/5))
   
   ## label Sample ID
   # if(samplesLength > 10000){
   #    p <- p + geom_text(aes(x = x2 + cos(angle)*textAdjust/10,
   #                           y = y2 + sin(angle)*textAdjust/10,
   #                           label = sample,
   #                           angle = angle*180/pi),
   #                       hjust = 0,
   #                       fontface = "bold",
   #                       data = treeData[(!sample %in% c("internal node",rootLabel)) &x2 > 0 & y2 !=max(treeData$y2), ],
   #                       size = sampleTextSize)
   #    p <- p + geom_text(aes(x = x2 + cos(angle)*textAdjust/10,
   #                           y = y2 + sin(angle)*textAdjust/10,
   #                           label = sample,
   #                           angle = angle*180/pi + 180),
   #                       hjust = 1,
   #                       fontface = "bold",
   #                       data = treeData[(!sample %in% c("internal node",rootLabel)) &x2 < 0& y2 !=max(treeData$y2), ],
   #                       size = sampleTextSize)
   #    ## label NORMAL
   #    p <- p + geom_text(aes(x = x2,y = y2-textAdjust/5),
   #                       label = rootLabel,
   #                       data = treeData[sample == rootLabel,], 
   #                       size = sampleTextSize)
   #    ## lable sample on the top
   #    if(nrow(treeData[which.max(treeData$y2), ]) > 0){
   #       p <- p + geom_text(aes(x = x2,
   #                              y = y2 + textAdjust/5,
   #                              label = sample),
   #                          data = treeData[which.max(treeData$y2), ],
   #                          size = sampleTextSize)
   #    } 
   # }
   # else{
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
   # }
   
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
         # p <- p + geom_point(aes(x =0 , y = 0), size = 1.7, color = "grey50",
         #                     fill = 'grey50', shape = 21, stroke = 0.5)
         p <- p + geom_point(aes(x =0 , y = 0), size = nodePointsSize, color = "#8c510a",
                             fill = 'white', shape = 21, stroke = nodeStrokeSize)
      }
   }
   if(show.bootstrap){
      if(use.box){
         p <- p + geom_label_repel(aes(x = x2, y = y2,label = boots),
                                   data = bootsData,
                                   # nudge_y = textAdjust/6,
                                   fontface = 'bold', size = bootLabelSize, box.padding = unit(bootPaddingSize, "lines"), point.padding = unit(0.5, "lines"),
                                   segment.colour = "grey50", segment.size = 0.25, force = 5)
      }else{
         p <- p + geom_text_repel(aes(x = x2, y = y2,label = boots),
                                  data = bootsData,
                                  fontface = 'bold', size = bootTextSize, box.padding = unit(bootPaddingSize, "lines"), point.padding = unit(0.5, "lines"),
                                  segment.colour = "grey50", segment.size = 0.25, force = 5)
      }
   }
   if(compare){
      p <- p + geom_label_repel(aes(x = x1 + (x2-x1)/2 , y = y1 + (y2 - y1)/2,label = is.match),
                                data = treeData[is.match != "NO",],
                                fontface = 'bold', size = bootLabelSize, box.padding = unit(bootPaddingSize, "lines"), point.padding = unit(0.5, "lines"),
                                segment.colour = "grey50", segment.size = 0.5, force = 5)
   }
   if(show.heatmap){
      ## combind heatmap and tree
      h <- plotHeatmap(binary.mat =  phyloTree@binary.matrix,
                       ccf.mat = phyloTree@ccf.matrix,
                       use.ccf = use.ccf,
                       show.class.label = show.class.label,
                       geneList = geneList,
                       plot.geneList = plot.geneList,
                       show.gene = show.gene,
                       show.geneList = show.geneList,
                       mut.threshold = mut.threshold)
      PH <- cowplot::plot_grid(p,
                               h,
                               rel_widths = c(1,0.5))
      # PH <- ggdraw(xlim = c(0.1,0.7)) +
      #     draw_plot(p, x = -0.05,y = 0, width = 0.7) +
      #     draw_plot(h, x = 0.48,y = -0.01, width = 0.2,height = 1)
      # PH <- ggdraw(xlim = c(0,1),ylim = c(0,1)) +
      #       draw_plot(p, x = 0,y = 0, width = 0.7) +
      #       draw_plot(h, x = 0.6,y = 0, width = 0.3)
      
      ## get title 
      pm <- getPrivateMutation(phyloTree = phyloTree)
      totalMutSum <- pm[[1]]
      privateMutProportion <- pm[[2]]
      title <- ggdraw() + 
         draw_label(paste(patientID,"\n(n = " ,totalMutSum ,"; ",privateMutProportion,")",sep = ""),fontface = "bold")
      
      ## combind heatmap,tree and title
      PH <- plot_grid(title,PH,ncol = 1,rel_heights=c(0.09, 1)) + 
         theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
      
      # ggsave(filename = paste0(patientID,".pdf"),plot = PH,width = 10,height = 6.5)
      return(PH)
   }
   else{
      return(p)
   }
}