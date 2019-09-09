#' A  phylogenetic tree painter
#' 
#' @import reshape2 ape ggplot2 deconstructSigs RColorBrewer ggrepel
#' 
#' @param njtree NJtree object
#' @param phylotree.type Phylotree format,you can choose "njtree","newick","beast","PAML" with root 
#' @param show.mutSig if show Mutational Signature in Images
#' @param sig.min.mut.number minimum mutation number in each branch
#' @param show.heatmap if plot heatmap that show mutation distribution in each branch
#' @param heatmap.type type of heatmap
#' @param phylotree.dat If the format of the phylotree is not "njtree", upload the path of the file to be analyzed with this parameter
#' 
#' @return Images of Phylotree
#' 
#' @export plotPhyloTree
#' 
#' @examples
#' plotPhyloTree(njtree)
#' plotPhyloTree(njtree)
#' # if use ccf 
#' plotPhyloTree(njtree, heatmap.type = '')
#' # use other tree format
#' newick.file <- system.file("extdata/newick", "1.nwk", package="MesKit")
#' plotPhyloTree(phylotree.dat = newick.file1, phylotree.type = 'newick')
#' beast.file <- system.file("extdata/BEAST", "sample.beast", package="MesKit")
#' plotPhyloTree(phylotree.dat = beast.file , phylotree.type = 'beast')
#' PAML.file <- system.file("extdata/PAML", "sample.paml", package="MesKit")
#' plotPhyloTree(phylotree.dat = PAML.file , phylotree.type = 'PAML')

## main  function
plotPhyloTree <- function(njtree = NULL, phylotree.type = 'njtree', 
                          show.mutSig = TRUE, sig.min.mut.number = 50, 
                          show.heatmap = TRUE, heatmap.type = 'binary', phylotree.dat = NULL){
  if(heatmap.type == 'binary'){
    use.ccf = FALSE
  }else{
    use.ccf = TRUE
    if (is.null(njtree@ccf_sort)){
      stop("Missing ccf file. Cannot generate heatmap of 'CCF' type")
    } 
  }
  
  if(phylotree.type != 'njtree'){
    show.heatmap = FALSE
    show.mutSig = FALSE
    if(phylotree.type == 'newick'){
      Phylo <- read.tree(phylotree.dat)
    }else if(phylotree.type == 'beast'){
      beast <- read.beast(phylotree.dat)
      Phylo <- as.phylo(beast)
    }
    else if(phylotree.type == 'PAML'){
      PAML <- read.paml_rst(phylotree.dat)
      Phylo <- as.phylo(PAML)
    }
    else{
      stop("the form of the tree file is not supported")
    }
  }
  else{
    if(is.null(njtree)){
      stop("njtree is null, you need to generate njtree using NJtree.R")
    }
    # PhyloTree input data
    Phylo <- njtree@nj
    refBuild <- njtree@refBuild
    signature <- suppressMessages(treeMutationalSig(njtree))
    njtree@patientID <- paste(njtree@patientID, ".NJtree", sep = "")
  }
  # generate phylotree data
  phylotree.input.data <- phylotreeInput(Phylo, signature, show.mutSig ,phylotree.type)
  phylotree.input.data <- phylotree.input.data[(phylotree.input.data$distance!=0|
                                                  phylotree.input.data$sample == 'NORMAL'),]
  
  
  if(show.mutSig){
    #Set the color
    color.scale <- colorSet(unique(phylotree.input.data$signature))
    njtree@patientID <- paste(njtree@patientID, ".mutsig", sep = "")
  }
  #plot phylotree
  phylotree <- generatePlotObject(phylotree.input.data, color.scale, show.mutSig, 
                                  phylotree.type)
  
  if(show.heatmap){
    heatmap <- mut.heatmap(njtree, use.ccf)
    plot.njtree <- ggdraw() + draw_plot(phylotree, x = 0,y = 0, width = 0.8) + draw_plot(heatmap, x = 0.8,y = -0.035, width = 0.2)
    return(plot.njtree)
  }
  else{
    return(phylotree)
  }
}
##generate plot data 
phylotreeInput <- function(phylo, signature = '', show.mutSig, phylotree.type){
  #Generate the adjacency matrix
  edge <-  phylo$edge
  #Add the third column of the matrix to the node distance
  distance <- phylo$edge.length
  edge <- matrix(c(edge, distance), nrow = length(edge[, 1]))
  # which root
  if(phylotree.type == 'njtree'){
    NO.Root <- which(phylo$tip.label == 'NORMAL')
    Root.tip <- NO.Root
    #the position of NORMAL in edge 
    Root.row <- which(edge[,2] == NO.Root)
    # the  Node connected to NORMAL
    Root.node <- edge[Root.row, 1]
    Root.label <- 'NORMAL'
  }
  else{
    NO.Root <- length(phylo$tip.label) + 1
    #the position of NORMAL in edge 
    Root.node <- NO.Root
    Root.label <- 'Root'
    Root.tip <- 0
    if(is.null(phylo$root.edge)){
      message('root egde is 0')
      phylo <- root(phylo, phylo$tip.label[floor(length(phylo$tip.label)/2)])
      if("Root" %in% phylo$node.label){
        stop("plot.PhyloTree draw tree with root only")
      }
      edge <-  phylo$edge
      distance <- phylo$edge.length
      edge <- matrix(c(edge, distance), nrow = length(edge[, 1]))
      NO.Root <- floor(length(phylo$tip.label)/2)
      #the position of NORMAL in edge 
      Root.row <- which(edge[,2] == NO.Root)
      # the  Node connected to NORMAL
      Root.node <- edge[Root.row, 1]
      Root.label <- 'Root'
      Root.tip <- 0
    }
    Root.edge <- phylo$root.edge
    edge <- rbind(edge, c(Root.node, Root.tip, Root.edge))
  }
  #store the sample name or node in a list
  name.list <- c('')
  #Store the list of nodes by plot order
  node.list <- c(Root.node)
  #X1, y1 are the starting point, x2 and y2 are the end point.Horizon stands for the Angle between the line and the positive x axis;W stands for the Angle occupied by branches (explained in the paper)
  #Distance means distance;Node is the internal node with the start;End_num is the number of each node in edge;
  plot.data <- data.frame('x1'=0, 'y1' = 0, 'x2' = 0, 'y2' = 0, 'horizon' = pi/2,
                          'w' = 1/6*pi, 'distance' = 0, 'node' = 0, 'end_num' = 0,
                          'horizon_angle' = 0)
  
  SetPhylotree <- function(sub.edge, x0, y0, w, horizon, name.list,
                           plot.data, point.list, target.node, Root.node, edge.left){
    #total vertex of tree
    vertex.total = 0
    #Stores the list of angles
    angle.list = c()
    #Stores the list of vertexes
    vertex.list = c()
    vertex.sub = 0
    angle = 0
    for(i in 1 : length(point.list)){
      #if it is sample, the number of vertices plus one
      if(point.list[i] <= length(phylo$tip.label)){
        vertex.total <- vertex.total + 1
        vertex.list <- append(vertex.list, 1)
      }
      #If it is an internal node, it counts several vertices connected to the internal node
      else{
        if(length(point.list[point.list > length(phylo$tip.label )]) == 1){
          vertex.total <- vertex.total+length(edge.left[, 1])
          vertex.sub <- length(edge.left[, 1])
        }else{
          for(i1 in 1 : length(edge.left[, 1])){
            if(point.list[i] %in% edge.left[i1, ]){
              vertex.total <- vertex.total + 1
              vertex.sub <- vertex.sub + 1
            } 
          }
        }
        vertex.list <- append(vertex.list, vertex.sub)
        vertex.sub <- 0
      }
    }
    #if(w <= pi/6 ){
    #w <- w + pi/6
    #}
    #The Angle assignment is based on the proportion of the number of vertices
    for(i in 1 : length(vertex.list)){
      angle.sub <- (vertex.list[i]/vertex.total)*w
      #Assign angles and put them into a list
      angle.list  <- append(angle.list, angle.sub)
    }
    #Calculate the Angle with the X axis
    angle.list2 <- c()
    for(i in 1:length(angle.list )){
      if(i == 1){
        angle=(pi-w)/2 + angle.list [i]/2 + horizon - pi/2
      }else{
        angle <- angle+angle.list [(i-1)]/2 + angle.list [i]/2 + horizon - pi/2
      }
      angle.list2 <- append(angle.list2, angle)
    }
    #Angle adjustment (select a branch Angle to align with the previous internal node)
    if(target.node != Root.node){
      a <- which(vertex.list == max(vertex.list))
      angle.list2[a[1]] <- horizon
    }
    i = 1
    while (i <= length(sub.edge[,1])) {
      if(point.list[i] > length(phylo$tip.label)){
        n <- 'internal node'
      }
      #If the number is greater than the length label of labels, it is an internal node, otherwise it is a sample
      else {n <- phylo$tip.label[point.list[i]]}
      if(phylotree.type!= 'njtree' & point.list[i]== Root.tip){
        n <- 'Root'
      }else if(phylotree.type == 'njtree' & point.list[i]== Root.tip){
        n <- 'NORMAL'
      }
      #Merge the sample name into a list for coloring
      name.list <- append(name.list, n)
      #get the distance 
      distance <- as.numeric(sub.edge[i,3])
      #Calculate the coordinate
      angle <- angle.list2[i]
      # if(angle <= 34*pi/36 & angle >= 30*pi/36){
      # angle = angle - pi/12
      # }
      # if(angle <= 8*pi/36 & angle >= 4*pi/36){
      #  angle = angle + pi/12
      # }
      x <- cos(angle)*distance + x0
      y <- sin(angle)*distance + y0
      #Merge data
      row <- list(x0, y0, x, y, angle, angle.list[i],
                  distance, target.node, point.list[i], horizon)
      plot.data <- rbind(plot.data, row)
      i <- i+1
    }
    #Return the data of the plotting and the list of sample names
    result <- list(plot.data, name.list)
    return(result)
  }
  t=1
  while(t <= phylo$Nnode){
    #The first diagram is of the internal node connected to NORMAL
    target.node <- node.list[t]
    #he position of the target t node in the matrix (the first data point is to find the NORMAL)
    row.list <- c()
    for(i1 in 1:length(edge[, 1])){
      if(target.node %in% edge[i1,1:2]){
        row.list <- append(row.list, i1)
      }
    }
    #Filter out the submatrix with target_node
    sub.edge <- matrix(edge[row.list, ], nrow = length(row.list))
    #Edge removes the sub_edge matrix
    edge <- edge[-row.list, ]
    #Extract the point that is concatenated with targe_node and put it in the list
    point.list <- c()
    for( i2 in 1 : length(sub.edge[, 1])){
      point <- sub.edge[i2, 1 : 2][sub.edge[i2, 1 : 2] != target.node]
      point.list <- append(point.list, point)
    }
    #first internal node(NORMAL)
    if(t == 1){
      p.n <- SetPhylotree(sub.edge, 0, 0, pi*2/3, pi/2, name.list, 
                          plot.data, point.list,target.node, Root.node, edge)
      plot.data <- p.n[[1]]
      name.list <- p.n[[2]]
    }
    else{
      #Find the current starting coordinate in the plot data
      row.num <- which(plot.data$end_num == target.node)
      p.n <- SetPhylotree(sub.edge, plot.data$x2[row.num], plot.data$y2[row.num],
                          plot.data$w[row.num], plot.data$horizon[row.num], name.list, plot.data, point.list,
                          target.node, Root.node, edge)
      plot.data <- p.n[[1]]
      name.list <- p.n[[2]]
    }
    node.list<-append(node.list, point.list[point.list>length(phylo$tip.label)])
    t=t+1
  }
  #Add the sample_name to the plot data
  plot.data <- cbind(plot.data, sample=name.list)
  #Transform NORMAL branches into vertical roots
  if (phylotree.type == 'njtree'){
    where.root <- which(plot.data$sample == Root.label)
    plot.data$x2[where.root] <- 0
    plot.data$y2[where.root] <- -plot.data$distance[where.root]
  }
  else{
    where.root <- which(plot.data$sample == Root.label)
    plot.data$x2 [where.root] <- 0
    plot.data$y2[where.root] <- -plot.data$distance[where.root]
  }
  
  plot.data <- plot.data[-1, ]
  #Adjust Angle of non-NORMAL sample and internal nodes
  t=1
  while(t <= length(plot.data[, 1])){
    if (all(plot.data$end_num[t] <= length(phylo$tip.label),
            plot.data$sample[t] != 'internal node',
            plot.data$sample[t] != Root.label,
            plot.data$horizon_angle[t] != plot.data$horizon[t])){
      #the angles less than ninety degrees are adjust to thirty degrees
      if(plot.data$horizon[t] < pi/2){
        plot.data$horizon[t] <- pi/6
        angle <- plot.data$horizon[t]
        plot.data$x2[t] <- plot.data$distance[t]*cos(angle) + plot.data$x1[t]
        plot.data$y2[t] <- plot.data$distance[t]*sin(angle) + plot.data$y1[t]}
      # above 90 are adjusted to 150
      else if(plot.data$horizon[t] > pi/2){
        plot.data$horizon[t] <- pi*5/6
        angle <- plot.data$horizon[t]
        plot.data$x2[t] <- plot.data$distance[t]*cos(angle) + plot.data$x1[t]
        plot.data$y2[t] <- plot.data$distance[t]*sin(angle) + plot.data$y1[t]}
    }
    t = t + 1
  }
  
  #right part of tree
  t=1
  sample.right <- plot.data[which(plot.data$sample != 'internal node' & 
                                    plot.data$sample != Root.label &
                                    plot.data$horizon <= pi/2 & 
                                    plot.data$horizon > 0), ]
  while(t <= length(sample.right[, 1])){
    u=1
    while(u <= (length(sample.right[, 1]) - 1)){
      if (sample.right$y1[u] == sample.right$y1[(u + 1)]){
        if (sample.right$horizon[u] < sample.right$horizon[(u + 1)]){
          u <- u+1
          next
        }else{
          data.1 <- sample.right[u, ]
          data.2 <- sample.right[(u+1), ]
          sample.right[(u+1), ] <- data.1
          sample.right[(u), ] <- data.2
          u <- u+1; next
        }
      }
      else if(sample.right$y1[u]>sample.right$y1[(u+1)]){
        data.1 <- sample.right[u, ]
        data.2 <- sample.right[(u + 1),]
        sample.right[(u+1), ] <- data.1
        sample.right[(u), ] <- data.2
        u <- u+1
        next
      }else{
        u <- u+1
        next
      }
      u <- u + 1
    }
    t <- t+ 1
  }
  list.right <- as.character(sample.right$sample)
  #left part of tree
  t <- 1
  sample.left <- plot.data[which(plot.data$sample!='internal node' & 
                                   plot.data$sample!=Root.label & 
                                   plot.data$horizon >= pi/2 & 
                                   plot.data$horizon < pi), ]
  while(t<=length(sample.left[, 1])){
    u=1
    while(u<=(length(sample.left[, 1]) - 1)){
      if(sample.left$y1[u] == sample.left$y1[(u+1)]){
        if(sample.left$horizon[u]<sample.left$horizon[(u+1)]){
          u=u+1
          next
        }else{
          data.1 <- sample.left[u, ]
          data.2 <- sample.left[(u+1), ]
          sample.left[(u+1), ] <- data.1
          sample.left[(u), ] <- data.2
          u <- u+1
          next
        }
      }
      else if (sample.left$y1[u] < sample.left$y1[(u + 1)]){
        data.1 <- sample.left[u, ]
        data.2 <- sample.left[(u+1), ]
        sample.left[(u+1), ] <- data.1
        sample.left[(u), ] <- data.2
        u <- u+1
        next
      }else{
        u <- u+1
        next
      }
      u <- u+1
    }
    t <- t+1
  }
  list.left <- as.character(sample.left$sample)
  
  #Adjust the Angle again
  t=1
  adjust.node.list <- c()
  while(t <= length(plot.data[,1])){
    row <- ''
    row1 <- ''
    if(all(plot.data$end_num[t] <= length(phylo$tip.label),
           plot.data$sample[t] != 'internal node',
           plot.data$sample[t] != Root.label,
           plot.data$horizon_angle[t] != plot.data$horizon[t])){
      internal.node <- plot.data$node[t]
      row1 <- which(plot.data$end_num == internal.node)
      if (length(row1) != 0){
        adjust.node <- plot.data$node[row1]
        if(adjust.node == Root.node){
          t <- t+1
          next
        }
        if(adjust.node %in% adjust.node.list){
          t <- t+1
          next
        }else{
          if(length(NO.Root)== 0){
            break
          }
          if(plot.data$distance[row1] < mean(phylo$edge.length)/10){
            row <- which(plot.data$node == adjust.node&
                           plot.data$end_num <= length(phylo$tip.label))
            if(length(row)!=0){
              if(plot.data$horizon[row] < pi/2 & plot.data$sample[row] == list.right[1]){
                plot.data$horizon[row] <- plot.data$horizon[row] - pi/18
                angle <- plot.data$horizon[row]
                plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
                plot.data$y2[row] <- plot.data$distance[row]*sin(angle) + plot.data$y1[row]
              }
              else if(plot.data$horizon[row]>pi/2 & 
                      plot.data$sample[row] == list.left[length(list.left)]){
                if(plot.data$horizon[row] == max(plot.data$horizon)){
                  plot.data$horizon[row] <- plot.data$horizon[row]+pi/18
                  angle <- plot.data$horizon[row]
                  plot.data$x2[row] <- plot.data$distance[row]*cos(angle)+plot.data$x1[row]
                  plot.data$y2[row] <- plot.data$distance[row]*sin(angle)+plot.data$y1[row]
                }
              }
            }
          }
          adjust.node.list <- append(adjust.node.list, adjust.node)
        }
      }
    }
    t = t + 1
  }
  ##Bisect the branches angle
  bisectBranchAngle <- function(list.part, plot.data){
    t = F
    while(t == F){
      if(length(list.part) == 1|length(list.part) == 2){
        break
      }
      if(length(list.part)==3){
        first.angle <-plot.data$horizon[which(plot.data$sample==list.part[1])] 
        final.angle <- plot.data$horizon[which(plot.data$sample==list.part[length(list.part)])]
        part.total.nodes <- length(list.part)
        average.angle <- (final.angle - first.angle)/part.total.nodes
        angle <- first.angle+average.angle
        for(i in 2 : (length(list.part) - 1)){
          row <- which(plot.data$sample == list.part[i])
          plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
          plot.data$y2[row] <- plot.data$distance[row]*sin(angle) + plot.data$y1[row]
          angle <- angle + average.angle
        }
        break
      }else{
        #Judge whether the adjustment area is on the left or right
        judge.angle <- plot.data$horizon[which(plot.data$sample == list.part[2])]
        if(judge.angle > pi/2){
          first.angle <-plot.data$horizon[which(plot.data$sample == list.part[1])] 
          final.angle <- plot.data$horizon[which(plot.data$sample == list.part[length(list.part) - 2])]
          part.total.nodes <- length(list.part)
          average.angle <- (final.angle - first.angle)/(part.total.nodes - 2)
          angle <- first.angle+average.angle
          for(i in 2:(length(list.part)-2)){
            row <- which(plot.data$sample == list.part[i])
            plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
            plot.data$y2[row] <- plot.data$distance[row]*sin(angle) + plot.data$y1[row]
            angle = angle+average.angle
          }
          break
        }
        else{
          first.angle <- plot.data$horizon[which(plot.data$sample == list.part[2])] 
          final.angle <- plot.data$horizon[which(plot.data$sample == list.part[length(list.part)])]
          part.total.nodes <- length(list.part)
          average.angle <- (final.angle - first.angle)/(part.total.nodes - 2)
          angle <- first.angle+average.angle
          for(i in 3 : (length(list.part) - 1)){
            row <- which(plot.data$sample == list.part[i])
            plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
            plot.data$y2[row] <- plot.data$distance[row]*sin(angle) + plot.data$y1[row]
            angle <- angle+average.angle
          }
          break
        }
      }
    }
    return(plot.data)
  }
  plot.data <- bisectBranchAngle(list.left, plot.data)
  plot.data <- bisectBranchAngle(list.right, plot.data)
  # Adjust the Angle of data with only two samples
  if(nrow(plot.data) == 3){
    angle.list <- c(pi/6, 5*pi/6)
    for(i in 2:nrow(plot.data)){
      angle <- angle.list[i - 1]
      plot.data$horizon[i] <- angle
      plot.data$x2[i] <- plot.data$distance[i]*cos(angle)
      plot.data$y2[i] <- plot.data$distance[i]*sin(angle)
    }
  }
  if(show.mutSig){
    plot.data <- addSignature(phylo, plot.data, signature)
  }
  
  return(plot.data)
}
##add signature
addSignature <- function(phylo, plot.data, signature){
  #add signature to plot.data
  plot.data$signature <- ''
  t <- 1
  while(t<=length(signature$branch)){
    if(length(signature$branch[[t]]) == 1){
      row <-which(plot.data$sample == signature$branch[[t]])
      plot.data$signature[row] <- as.character(signature$sig[[t]])     
      t <- t + 1
      next
    }
    #find Find the row corresponding to the sample or internal in the plot.data
    row.list <- c()
    for(i in 1:length(signature$branch[[t]])){
      sample <- signature$branch[[t]][i]
      row <- which(plot.data$sample == sample)
      row.list <- append(row.list, row)
    }
    # the sample in sample_list share one common node(previous node  to the lowest row )
    # signature of NORMAL
    if(length(signature$branch[[t]]) == length(phylo$tip.label ) - 1){
      row <- which(plot.data$sample == 'NORMAL')
      plot.data$signature[row] <-  as.character(signature$sig[(t)]) 
      t <- t+1
      next
    }
    #Internal nodes connected to NORMAL
    else if(length(signature$branch[[t]]) == length(phylo$tip.label ) - 2){
      for(i in 1 : length(plot.data$sample)){
        if(plot.data$end_num[i] > length(phylo$tip.label)){
          command.row <- i
          break
        }
      }
      plot.data$signature[command.row] <-  as.character(signature$sig[(t)]) 
      t=t+1
      next
    }
    # if signature does not belong to non- NORMAL sample
    if(length(signature$branch[[t]]) != 1){
      lowest.row <- min(row.list)
      command.row <- which(plot.data$end_num == plot.data$node[lowest.row]) 
      if(length(command.row) != 0){
        plot.data$signature[command.row] <- as.character(signature$sig[(t)]) 
      }
    }else{
      row <- which(plot.data$sample == signature$branch[[t]])
      plot.data$signature[row] <-  as.character(signature$sig[(t)]) 
    }
    t = t + 1
  }
  if(plot.data$signature[which(plot.data$sample == 'NORMAL')] == ''){
    plot.data$signature[which(plot.data$sample == 'NORMAL')] = as.character(signature$sig[1]) 
  }
  plot.data <- plot.data[order(plot.data$signature), ]
  plot.data$signature <- gsub('No.Signature', 'No signature', plot.data$signature)
  plot.data$signature <- gsub('Signature.', '', plot.data$signature)
  return(plot.data)
}
##color scale set
colorSet <- function(signatures){
  all.color.scale <- c("#E41A1C","#377EB8","#4DAF4A","#66C2A5","#FC8D62",
                       "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
                       "#8DD3C7", "#6E016B" ,"#BEBADA", "#FB8072", "#80B1D3",
                       "#FDB462","#B3DE69","#FCCDE5","#91003F","#BC80BD",
                       "#CCEBC5" ,"#FFED6F","#1B9E77", "#D95F02", "#7570B3" ,
                       "#E7298A" ,"#66A61E" ,"#E6AB02" ,"#A6761D", "#666666",'black')
  all.signature <- append(gsub('Signature.', '',row.names(signatures.cosmic)), 'No signature')
  color.scale <- c()
  for(i in 1:length(signatures)){
    color <- all.color.scale[which(all.signature == signatures[i])]
    color.scale <- append(color.scale, color)
  }
  if('black' %in% color.scale){
    num <- which(color.scale == 'black')
    color.scale <- color.scale[-num]
    color.scale <- append(color.scale, 'black')
  }
  return(color.scale)
}
## plot PhyloTree 
generatePlotObject <- function(plot.data, color.scale = '', show.mutSig, phylotree.type){
  p <- ggplot(data = plot.data)
  text.adjust <- mean(as.numeric(plot.data$distance))
  if(phylotree.type == 'njtree'){
    Root.label <- 'NORMAL'
  }
  else{
    Root.label <- 'Root'
  }
  if(show.mutSig){
    p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = signature), size=1.5)
    #the color of signature.1~signature.30 is
    #   ["#E41A1C" "#377EB8" "#4DAF4A""#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
    #   "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" "#D9D9D9" "#BC80BD"
    #    "#CCEBC5" "#FFED6F""#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"]
    p <- p + scale_color_manual(values = color.scale)
  }else{
    p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = sample),
                          data = plot.data[plot.data$sample != Root.label&plot.data$sample != 'internal node',], 
                          size=1.5, show.legend = T)
    p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = 'black',
                          data = plot.data[plot.data$sample == Root.label,], 
                          size = 1.5, show.legend = F )
    p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = '#67001F',
                          data = plot.data[plot.data$sample == 'internal node',], 
                          size = 1.5, show.legend = F )
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
                 legend.position = 'right',
                 legend.title = element_text(face="bold")
  )
  p <- p + geom_text_repel(aes(x = x2*0.95, y = y2, label = sample),vjust = 1,
                           nudge_y = text.adjust/7, segment.alpha = 0,
                           data = plot.data[(plot.data$sample != 'internal node'&plot.data$sample != Root.label), ],
                           fontface = 'bold', size = 4)
  p <- p + geom_text_repel(aes(x = x2,y = y2), label = Root.label, vjust = 0, 
                           nudge_y = -text.adjust/7, segment.alpha = 0,
                           data = plot.data[(plot.data$sample == Root.label),], 
                           fontface = 'bold', size = 4)
  #Leaf nodes and internal nodes are distinguished by point size
  p <- p + geom_point(aes(x = x2,y = y2), 
                      data = plot.data[plot.data$sample == 'internal node',],
                      size = 1.7, color = "#67001F", fill = '#67001F', shape = 21, 
                      stroke = 0.5)
  p <- p + geom_point(aes(x = x2, y = y2), 
                      data = plot.data[plot.data$sample != 'internal node',],
                      size = 3, color = "#67001F", fill = 'white', shape = 21, stroke = 1)
  Nd <- plot.data$distance[which(plot.data$sample == Root.label)]
  if(length(Nd)!=0){
    if(Nd != 0){
      p <- p + geom_point(aes(x =0 , y = 0), size = 1.7, color = "#67001F",
                          fill = '#67001F', shape = 21, stroke = 0.5)
    }
  }
  return(p)
}


