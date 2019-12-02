#' A  phylogenetic tree painter
#' 
#' @import reshape2 ape ggplot2 deconstructSigs RColorBrewer ggrepel
#' 
#' @param njtree NJtree object
#' @param show.mutSig If show Mutational Signature on tree.Default is "True"
#' @param sig.name Choose "alias"  when you want showing alias in branch.Default is "default"
#' @param show.heatmap If plot heatmap showing mutation distribution in each branch.Default is True
#' @param heatmap.type Type of heatmap,choose 'CCF' if use ccf.Default is True
#' 
#' @return Images of Phylotree
#' 
#' @export plotPhyloTree
#' 
#' @examples
#' plotPhyloTree(njtree)
#' plotPhyloTree(njtree)
#' ## Use ccf 
#' plotPhyloTree(njtree, heatmap.type = 'CCF')

plotPhyloTree <- function(njtree = NULL, show.mutSig = TRUE,sig.name = "default",show.heatmap = TRUE, heatmap.type = 'binary'){
  if(heatmap.type == 'binary'){
    use.ccf = FALSE
  }else{
    use.ccf = TRUE
  }
  # if(is.null(njtree)){
  #     stop("You need to generate njtree using NJtree.R")
  # }
  phylo <- njtree@nj
  refBuild <- njtree@refBuild
  signature <- treeMutationalSig(njtree)
  patientID <- njtree@patientID
  fileID <- paste(njtree@patientID, ".NJtree", sep = "")
  Root.label <- 'NORMAL'
  # generate phylotree data
  phylotree.input.data <- phylotreeInput(phylo, signature, show.mutSig ,Root.label)
  phylotree.input.data <- phylotree.input.data[(phylotree.input.data$distance!=0|phylotree.input.data$sample == Root.label),]
  if(show.mutSig){
    #Set the color
    color.scale <- colorSet(unique(phylotree.input.data$signature))
    fileID <- paste(fileID, ".mutsig", sep = "")
  }
  #plot phylotree
  phylotree <- generatePlotObject(phylotree.input.data, color.scale, show.mutSig, sig.name = sig.name, Root.label = Root.label)
  if(show.heatmap){
    heatmap <- mut.heatmap(njtree, use.ccf)
    pm <- getPrivateMutation(njtree)
    totalMut.sum <- pm[[1]]
    privateMut.proportion <- pm[[2]]
    PH <- ggdraw(xlim = c(0.1,0.7)) + draw_plot(phylotree, x = -0.05,y = 0, width = 0.7) + draw_plot(heatmap, x = 0.46,y = -0.12, width = 0.15)
    title <- ggdraw() + draw_label(paste(patientID,"\n(n = " ,totalMut.sum ,"; ",privateMut.proportion,")",sep = ""),fontface = "bold")
    PH <- plot_grid(title,PH,ncol = 1,rel_heights=c(0.09, 1))+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
    # ggsave(filename = "./tree.pdf",plot = PH, width = 10, height = 6.5)
    return(PH)
  }
  else{
    return(phylotree)
  }
}
##generate plot data 
phylotreeInput <- function(phylo, signature = '', show.mutSig, Root.label){
  phylo <- ape::root(phylo, phylo$tip.label[which(phylo$tip.label == Root.label)])
  # phylo$edge[which(phylo$edge == phylo$edge[1,2])] = phylo$edge[1,1]
  # phylo$edge <- phylo$edge[-1,]
  # phylo$edge.length <- phylo$edge.length[-1]
  # phylo$Nnode <- phylo$Nnode - 1
  phylo$edge.length <- phylo$edge.length/100
  edge <-  phylo$edge
  distance <- phylo$edge.length
  edge <- matrix(c(edge, distance), nrow = length(edge[, 1]))
  NO.Root <- which(phylo$tip.label == Root.label)
  Root.tip <- NO.Root
  #the position of NORMAL in edge 
  Root.row <- which(edge[,2] == NO.Root)
  # the  Node connected to NORMAL
  Root.node <- edge[Root.row, 1]
  Root.edge <- edge[Root.row, 3]
  verticalPath <- calVerticalPath(phylo,Root.label)
  branch.label <- labelBranch(phylo)
  #store the sample name or node in a list
  name.list <- c('')
  #Store the list of nodes by plot order
  node.list <- c(Root.node)
  #X1, y1 are the starting point, x2 and y2 are the end point.Horizon stands for the Angle between the line and the positive x axis;W stands for the Angle occupied by branches (explained in the paper)
  #Distance means distance;Node is the internal node with the start;End_num is the number of each node in edge;
  plot.data <- data.frame('x1'=0, 'y1' = 0, 'x2' = 0, 'y2' = 0, 'horizon' = pi/2,
                          'w' = 1/6*pi, 'distance' = 0, 'node' = 0, 'end_num' = 0,
                          'angle' = 0)
  
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
        }
        else{
          for(i1 in 1 : length(edge.left[, 1])){
            if(point.list[i] %in% edge.left[i1, ]){
              vertex.total <- vertex.total + 1
              vertex.sub <- vertex.sub + 1
            } 
          }
        }
        vertex.list <- append(vertex.list, (vertex.sub))
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
    for(i in 1:length(angle.list)){
      if(i == 1){
        angle=(pi-w)/2 + angle.list [i]/2 + horizon - pi/2
      }else{
        angle <- angle+angle.list [(i-1)]/2 + angle.list [i]/2 + horizon - pi/2
      }
      # if(angle == horizon){
      #   if(horizon >= pi/2){
      #     angle <- angle + pi/144
      #   }
      #   else{
      #     angle <- angle - pi/144
      #   }
      # }
      angle.list2 <- append(angle.list2, angle)
    }
    # x <- which(!(sub.edge[,2] %in% verticalPath))
    # if(horizon == pi/2){
    #   if(angle.list2[x] < pi/2 &
    #      as.numeric(sub.edge[x,2]) > length(phylo$tip.label)){
    #     angle.list2[x] <- angle.list2[x] - pi/18
    #   }
    #   else if(angle.list2[x] > pi/2 &
    #           as.numeric(sub.edge[x,2]) > length(phylo$tip.label)){
    #     angle.list2[x] <- angle.list2[x] + pi/18
    #   }
    # }
    ## Angle adjustment (select a branch Angle to align with the previous internal node)
    a <- which(sub.edge[,2] %in% verticalPath)
    if(length(a) != 0){
      if(target.node == Root.node){
        # if(sub.edge[a,2] > length(phylo$tip.label)){
        #   angle.list2[a[1]] <- horizon
        # }
        # else{
        #   a <- which(vertex.list == max(vertex.list))
        #   if(length(a) > 1){
        #     a <- which(sub.edge[,3] == max(sub.edge[,3]))
        #   }
        #   angle.list2[a[1]] <- horizon
        # }
        angle.list2[a[1]] <- horizon
      }
      else{
        angle.list2[a[1]] <- horizon
        ##xxx
        left <- nrow(plot.data[(plot.data$angle > pi/2&plot.data$end_num!= NO.Root),])
        right <- nrow(plot.data[(plot.data$angle < pi/2&plot.data$end_num!= NO.Root),])-1
        if(left > right){
          if(angle.list2[which(angle.list2 != horizon)]> pi/2){
            angle.list2[which(angle.list2 != horizon)] <- abs(pi - angle.list2[which(angle.list2 != horizon)])
          }
        }
        else if(left < right){
          if(angle.list2[which(angle.list2 != horizon)] < pi/2){
            angle.list2[which(angle.list2 != horizon)] <- abs(pi - angle.list2[which(angle.list2 != horizon)])
          }
        }
      }
    }
    else{
      a <- which(vertex.list == max(vertex.list))
      if(length(a) > 1){
        a <- which(sub.edge[,3] == max(sub.edge[,3]))
      }
      angle.list2[a[1]] <- horizon
    }
    x <- which(angle.list2!= pi/2)
    if(length(x) == 1){
      if(horizon == pi/2){
        if(angle.list2[x] < pi/2 &
           as.numeric(sub.edge[x,2]) > length(phylo$tip.label)){
          angle.list2[x] <- angle.list2[x] - pi/18
        }
        else if(angle.list2[x] > pi/2 &
                as.numeric(sub.edge[x,2]) > length(phylo$tip.label)){
          angle.list2[x] <- angle.list2[x] + pi/18
        }
      }
    }
    i = 1
    while (i <= length(sub.edge[,1])) {
      if(point.list[i] > length(phylo$tip.label)){
        n <- 'internal node'
      }
      #If the number is greater than the length label of labels, it is an internal node, otherwise it is a sample
      else {
        n <- phylo$tip.label[point.list[i]]
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
      if(angle == pi/2){
        x <- x0
      }
      else{
        x <- cos(angle)*distance + x0
      }
      y <- sin(angle)*distance + y0
      #Merge data
      row <- list(x0, y0, x, y, horizon, angle.list[i],
                  distance, target.node, point.list[i], angle)
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
        if(edge[i1,2] != NO.Root){
          row.list <- append(row.list, i1)
        }
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
                          plot.data, point.list, target.node, Root.node, edge)
      plot.data <- p.n[[1]]
      name.list <- p.n[[2]]
      Root.info <- list(0, 0, 0, -(Root.edge), 0, pi/2, Root.edge, Root.node, NO.Root, 0)
      plot.data <- rbind(plot.data, Root.info)
      name.list <- append(name.list, "NORMAL")
    }
    else{
      #Find the current starting coordinate in the plot data
      row.num <- which(plot.data$end_num == target.node)
      p.n <- SetPhylotree(sub.edge, plot.data$x2[row.num], plot.data$y2[row.num],
                          plot.data$w[row.num], plot.data$angle[row.num], name.list, plot.data, point.list,
                          target.node, Root.node, edge)
      plot.data <- p.n[[1]]
      name.list <- p.n[[2]]
    }
    node.list<-append(node.list, point.list[point.list>length(phylo$tip.label)])
    t=t+1 
  }
  #Add the sample_name to the plot data
  plot.data <- cbind(plot.data, sample = name.list)
  plot.data <- plot.data[-1, ]
  #Adjust Angle of non-NORMAL sample and internal nodes
  t=1
  while(t <= length(plot.data[, 1])){
    if (all(plot.data$end_num[t] <= length(phylo$tip.label),
            plot.data$sample[t] != 'internal node',
            plot.data$sample[t] != Root.label,
            plot.data$angle[t] != plot.data$horizon[t])){
      #the angles less than ninety degrees are adjust to thirty degrees
      if(plot.data$x1[t] > 0){
        plot.data$angle[t] <- pi*(1/6-1/36) 
        angle <- plot.data$angle[t]
        plot.data$x2[t] <- plot.data$distance[t]*cos(angle) + plot.data$x1[t]
        plot.data$y2[t] <- plot.data$distance[t]*sin(angle) + plot.data$y1[t]}
      else if(plot.data$x1[t] < 0){
        plot.data$angle[t] <- pi*(5/6+1/36) 
        angle <- plot.data$angle[t]
        plot.data$x2[t] <- plot.data$distance[t]*cos(angle) + plot.data$x1[t]
        plot.data$y2[t] <- plot.data$distance[t]*sin(angle) + plot.data$y1[t]}
      else{
        if(plot.data$angle[t] < pi/2){
          plot.data$angle[t] <- pi*(1/6-1/36) 
          angle <- plot.data$angle[t]
          plot.data$x2[t] <- plot.data$distance[t]*cos(angle) + plot.data$x1[t]
          plot.data$y2[t] <- plot.data$distance[t]*sin(angle) + plot.data$y1[t]}
        else if(plot.data$angle[t] > pi/2){
          plot.data$angle[t] <- pi*(5/6+1/36)
          angle <- plot.data$angle[t]
          plot.data$x2[t] <- plot.data$distance[t]*cos(angle) + plot.data$x1[t]
          plot.data$y2[t] <- plot.data$distance[t]*sin(angle) + plot.data$y1[t]}
      }
    }
    t = t + 1
  }
  #right part of tree
  t=1
  sample.right <- plot.data[which(plot.data$sample != 'internal node' & 
                                    plot.data$sample != Root.label &
                                    plot.data$angle <= pi/2 & 
                                    plot.data$angle > 0&
                                    plot.data$angle != plot.data$horizon), ]
  while(t <= length(sample.right[, 1])){
    u=1
    while(u <= (length(sample.right[, 1]) - 1)){
      if (sample.right$y1[u] == sample.right$y1[(u + 1)]){
        if (sample.right$angle[u] < sample.right$angle[(u + 1)]){
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
                                   plot.data$angle >= pi/2 & 
                                   plot.data$angle < pi&
                                   plot.data$angle != plot.data$horizon), ]
  while(t<=length(sample.left[, 1])){
    u=1
    while(u<=(length(sample.left[, 1]) - 1)){
      if(sample.left$y1[u] == sample.left$y1[(u+1)]){
        if(sample.left$angle[u]<sample.left$angle[(u+1)]){
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
           plot.data$angle[t] != plot.data$horizon[t])){
      # if(nrow(plot.data[plot.data$x2 > 0,]) == 2|nrow(plot.data[plot.data$x2 < 0,]) == 2){
      #   break
      # }
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
          if(plot.data$distance[row1] < mean(phylo$edge.length)/5){
            row <- which(plot.data$node == adjust.node&
                           plot.data$end_num <= length(phylo$tip.label))
            row2 <- which(plot.data$node == plot.data$end_num[row1] &
                            plot.data$angle != plot.data$horizon &  
                            plot.data$end_num <= length(phylo$tip.label))
            if(length(row)!=0 & length(row2)!=0){
              if(plot.data$angle[row] < pi/2 &
                 plot.data$angle[row2] < pi/2 &
                 plot.data$sample[row] == list.right[1]){
                if(nrow(plot.data[plot.data$x2 > 0,]) == 2){
                  break
                }
                plot.data$angle[row] <- plot.data$angle[row] - pi*1/10
                angle <- plot.data$angle[row]
                plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
                plot.data$y2[row] <- plot.data$distance[row]*sin(angle) + plot.data$y1[row]
              }
              else if(plot.data$angle[row] > pi/2 &
                      plot.data$angle[row2] > pi/2 &
                      plot.data$sample[row] == list.left[length(list.left)]){
                if(nrow(plot.data[plot.data$x2 < 0,]) == 2){
                  break
                }
                if(plot.data$angle[row] == max(plot.data$angle)){
                  plot.data$angle[row] <- plot.data$angle[row] + pi*1/10
                  angle <- plot.data$angle[row]
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
      # if(length(list.part) == 1|length(list.part) == 2){
      #   break
      # }
      if(length(list.part) == 1){
        break
      }
      else if(length(list.part) == 2){
        judge.angle <- plot.data$angle[which(plot.data$sample == list.part[1])]
        if(judge.angle > 0 & judge.angle < pi/2 & nrow(plot.data[plot.data$x2 > 0,]) == 2){
          first.angle <- plot.data$angle[which(plot.data$sample==list.part[1])] 
          final.angle <- pi/2
          row <- which(plot.data$sample == list.part[2])
          angle <- (final.angle - first.angle)/2 + first.angle
          plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
          plot.data$y2[row] <- plot.data$distance[row]*sin(angle) + plot.data$y1[row]
          break
        }
        else if(judge.angle > pi/2 & judge.angle < pi & nrow(plot.data[plot.data$x2 < 0,]) == 2){
          first.angle <- pi/2
          final.angle <- plot.data$angle[which(plot.data$sample==list.part[2])] 
          row <- which(plot.data$sample == list.part[1])
          angle <- (final.angle - first.angle)/2 + first.angle
          plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
          plot.data$y2[row] <- plot.data$distance[row]*sin(angle) + plot.data$y1[row]
          break
        }
        else{
          break
        }
      }
      else if(length(list.part)==3){
        first.angle <-plot.data$angle[which(plot.data$sample==list.part[1])] 
        final.angle <- plot.data$angle[which(plot.data$sample==list.part[length(list.part)])]
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
        judge.angle <- plot.data$angle[which(plot.data$sample == list.part[2])]
        if(length(judge.angle != 0)){
          if(judge.angle > pi/2){
            first.angle <-plot.data$angle[which(plot.data$sample == list.part[1])] 
            final.angle <- plot.data$angle[which(plot.data$sample == list.part[length(list.part) - 2])]
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
            first.angle <- plot.data$angle[which(plot.data$sample == list.part[2])] 
            final.angle <- plot.data$angle[which(plot.data$sample == list.part[length(list.part)])]
            part.total.nodes <- length(list.part)
            average.angle <- (final.angle - first.angle)/(part.total.nodes - 2)
            angle <- first.angle+average.angle
            for(i in 3 : (length(list.part) - 1)){
              row <- which(plot.data$sample == list.part[i])
              plot.data$x2[row] <- plot.data$distance[row]*cos(angle) + plot.data$x1[row]
              plot.data$y2[row] <- plot.data$distance[ row]*sin(angle) + plot.data$y1[row]
              angle <- angle+average.angle
            }
            break
          }
        }
        else{break}
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
      plot.data$angle[i] <- angle
      plot.data$x2[i] <- plot.data$distance[i]*cos(angle)
      plot.data$y2[i] <- plot.data$distance[i]*sin(angle)
    }
  }
  if(show.mutSig){
    plot.data$label <- ""
    for(i in 1:nrow(plot.data)){
      if(plot.data$sample[i] == "NORMAL"){
        plot.data$label[i] <- branch.label[[Root.node]]
      }
      else{
        if(plot.data$end_num[i] > length(phylo$tip.label)){
          plot.data$label[i] <- branch.label[[plot.data$end_num[i]]]
        }
        else{
          plot.data$label[i] <- as.character(plot.data$sample[i]) 
        }
      }
    }
    plot.data <- addSignature(phylo, plot.data, signature)
  }
  return(plot.data)
}
##add signature
addSignature <- function(phylo, plot.data, signature){
  #add signature to plot.data
  plot.data$signature <- ''
  plot.data$alias <- ''
  sigs <- strsplit(as.character(signature$Branch),"∩")
  sigs <- lapply(sigs, function(x){return(paste(sort(x,decreasing = T),collapse = "∩"))})
  t <- 1
  while(t<=length(sigs)){
    pos <- which(plot.data$label == sigs[t])
    plot.data$signature[pos] <- as.character(signature$Signature[t]) 
    plot.data$alias[pos] <- as.character(signature$Alias[t])
    t <- t + 1
  }
  if(plot.data$signature[which(plot.data$sample == 'NORMAL')] == ''){
    plot.data$signature[which(plot.data$sample == 'NORMAL')] = as.character(signature$Signature[1])
    plot.data$signature[which(plot.data$alias == 'NORMAL')] = as.character(signature$Alias[1])
  }
  plot.data <- plot.data[order(plot.data$signature), ]
  plot.data$signature <- gsub('No.Signature', 'No signature', plot.data$signature)
  plot.data$signature <- gsub('Signature.', '', plot.data$signature)
  return(plot.data)
}
##color scale set
colorSet <- function(signatures){
  ## FF6A5A:Sig19
  all.color.scale <- c("#E41A1C","#377EB8","#7F0000",
                       "#35978f","#FC8D62","#2166ac",
                       "#E78AC3","#A6D854","#FFD92F",
                       "#E5C494","#8DD3C7", "#6E016B" ,
                       "#BEBADA", "#e08214", "#80B1D3",
                       "#d6604d","#ffff99","#FCCDE5",
                       "#FF6A5A","#BC80BD","#CCEBC5" ,
                       "#fb9a99","#B6646A", "#9F994E", 
                       "#7570B3" ,"#c51b7d" ,"#66A61E" ,
                       "#E6AB02" ,"#003c30", "#666666",'black')
  all.signature <- append(gsub('Signature.', '',row.names(deconstructSigs::signatures.cosmic)), 'No signature')
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
generatePlotObject <- function(plot.data, color.scale = '', show.mutSig, sig.name, Root.label){
  p <- ggplot(data = plot.data)
  text.adjust <- mean(as.numeric(plot.data$distance))
  dy <- max(plot.data$y2)-min(plot.data$y2)
  dx <- max(plot.data$x2)-min(plot.data$x2)
  if(show.mutSig){
    p <- p + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = signature), size=1.5)
    #the color of signature.1~signature.30 is
    #   ["#E41A1C" "#377EB8" "#4DAF4A""#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
    #   "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" "#D9D9D9" "#BC80BD"
    #    "#CCEBC5" "#FFED6F""#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"]
    p <- p + scale_color_manual(values = color.scale)
  }
  else{
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
                 legend.title = element_text(face="bold"),
                 legend.position = 'top',
                 legend.direction = "horizontal") + guides(color = guide_legend(nrow=1))+coord_fixed(ratio= 1)               
  if(sig.name == "default"){
    p <- p + geom_text_repel(aes(x = x2*1.5, y = y2, label = sample),vjust = 1,
                             nudge_y = text.adjust/5, segment.alpha = 0,
                             data = plot.data[(plot.data$sample != 'internal node'&plot.data$sample != Root.label), ],
                             fontface = 'bold', size = 3)
  }
  else{
    p <- p + geom_text_repel(aes(x = x2*1.2, y = y2, label = alias),vjust = 1,
                             nudge_y = text.adjust/8, segment.alpha = 0,
                             data = plot.data[(plot.data$sample != 'internal node'&plot.data$sample != Root.label), ],
                             fontface = 'bold', size = 4) 
  }
  if(sig.name == "default"){
    p <- p + geom_text_repel(aes(x = x2,y = y2), label = Root.label, vjust = 0, 
                             nudge_y = -text.adjust/5, segment.alpha = 0,
                             data = plot.data[(plot.data$sample == Root.label),], 
                             fontface = 'bold', size = 3)
  }
  else{
    p <- p + geom_text_repel(aes(x = x2,y = y2), label = "T", vjust = 0, 
                             nudge_y = -text.adjust/5, segment.alpha = 0,
                             data = plot.data[(plot.data$sample == Root.label),], 
                             fontface = 'bold', size = 4)
  }
  #Leaf nodes and internal nodes are distinguished by point size
  p <- p + geom_point(aes(x = x2,y = y2), 
                      data = plot.data[plot.data$sample == 'internal node',],
                      size = 1.7, color = "#8c510a", fill = "#8c510a", shape = 21, 
                      stroke = 0.5)
  p <- p + geom_point(aes(x = x2, y = y2), 
                      data = plot.data[plot.data$sample != 'internal node',],
                      size = 3,color = "#67001F", fill = 'white', shape = 21, stroke = 1)
  # if(phylotree.type!= "njtree"){
  #   p <- p + geom_point(x = 0, y = 0, 
  #                       size = 3, color = "#8c510a", fill = 'white', shape = 21, stroke = 1)
  # }
  Nd <- plot.data$distance[which(plot.data$sample == Root.label)]
  if(length(Nd)!=0){
    if(Nd != 0){
      p <- p + geom_point(aes(x =0 , y = 0), size = 1.7, color = "#8c510a",
                          fill = '#8c510a', shape = 21, stroke = 0.5)
    }
  }
  return(p)
}

getPrivateMutation <- function(njtree){
  totalMut <- njtree@mut_branches
  private.order <- unlist(lapply(names(totalMut),
                                 function(x){return(length(strsplit(x,"b)")[[1]]) == 1)})) 
  privateMut <- totalMut[private.order]
  countMutation <- function(mut){
    sum <- 0
    for(i in 1:length(mut)){
      count <- nrow(mut[[i]])
      sum <- sum + count
    }
    return(sum)
  }
  totalMut.sum <- countMutation(totalMut)
  privateMut.sum <- countMutation(privateMut)
  privateMut.proportion <- paste(round((privateMut.sum/totalMut.sum)*100,1),"%",sep = "")
  return(list(totalMut.sum, privateMut.proportion))
}

calVerticalPath <- function(phylo,Root.label){
  #Generate the adjacency matrix
  edge <-  phylo$edge
  #Add the third column of the matrix to the node distance
  distance <- phylo$edge.length
  edge <- matrix(c(edge, distance), nrow = length(edge[, 1]))
  if(Root.label == "NORMAL"){
    NO.Root <- which(phylo$tip.label == Root.label)
    NO.Root <- which(phylo$tip.label == Root.label)
    #the position of NORMAL in edge 
    Root.row <- which(edge[,2] == NO.Root)
    # the  Node connected to NORMAL
    Root.node <- edge[Root.row, 1]
    edge <- edge[-Root.row,]
  }
  else{
    NO.Root <- length(phylo$tip.label) + 1
    # the  Node connected to NORMAL
    Root.node <- NO.Root
    Root.tip <- NO.Root
  }
  distance.table <- data.frame(x = 0)
  for(i in 1:(length(phylo$tip.label)-1)){
    distance.table <- cbind(distance.table, -1)
  }
  distance.table <- distance.table[,(-1)]
  if(Root.label == "NORMAL"){
    names(distance.table) <- phylo$tip.label[which(phylo$tip.label != Root.label)]
  }
  else{
    distance.table <- cbind(distance.table, -1)
    names(distance.table) <- phylo$tip.label
  }
  t <- 1
  path <- list()
  while(t <= length(distance.table)){
    subpath <- c()
    end <- t
    dis <- 0
    while (TRUE) {
      p <- which(edge[,2] == end)
      dis <- dis + edge[p,3]
      start <- edge[p,1]
      subpath <- append(subpath,end)
      end <- start
      if(start == Root.node){
        break
      }
    }
    path[[t]] <- subpath
    distance.table[t] <- dis
    t <- t + 1
  }
  result <- path[[which.max(distance.table)]]
}
## label sample in each branch                            
labelBranch <- function(phylo){
  Root <- which(phylo$tip.label == "NORMAL")
  internalNodes <- sort(unique(phylo$edge[,1]))
  result <- list()
  for(i in 1:(length(phylo$edge.length)+1)){
    result[[i]] <- NA
  }
  end <- phylo$edge[which(phylo$edge[,2] == Root),1]
  for(i in 1:length(phylo$tip.label)){
    row <- phylo$edge[which(phylo$edge[,2] == i), ]
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
        row <- phylo$edge[which(phylo$edge[,2] == node),]
      }
    }
  }
  g <- lapply(result, function(x){return(paste(sort(phylo$tip.label[x[-1]], decreasing = T),collapse = "∩"))})
  return(g)
}

