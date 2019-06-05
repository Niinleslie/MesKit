#' A object of NJtree painter
#' 
#' @import reshape2 reshape ape ggplot2 deconstructSigs RColorBrewer
#' 
#' @param maf MAF object
#' @return NJtree plot (phylogenetic tree and heatmap)
#' 
#' @examples
#' maf <- read.Maf("311252",dat.dir = './data/multi_lesion', BSG = "BSgenome.Hsapiens.UCSC.hg19")
#' plot.NJtree(maf)



#pakages
library(cowplot)
library(reshape2)
library(reshape)
library(ape)
library(ggplot2)
library(deconstructSigs)
library(RColorBrewer)
library(methods)
library(reshape2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomeInfoDb)
library(grDevices)
library(graphics)
library(utils)
# data frame needed
library(plyr)
library(ggrepel)

##generate plot data 
PhyloTree.input <-function(njtree){
  #Generate the adjacency matrix
  edge <- njtree@nj$edge
  #Add the third column of the matrix to the node distance
  distance <- njtree@nj$edge.length
  edge <- matrix(c(edge, distance), nrow = length(edge[,1]))
  #num of normal 
  NO_normal <- which(njtree@nj$tip.label == 'normal')
  #the position of normal in edge 
  normal_row <- which(edge[,2] == NO_normal)
  # the  Node connected to normal
  normal_node <- edge[normal_row, 1]
  #store the sample name or node in a list
  name_list <- c('')
  #Store the list of nodes by plot order
  node_list <- c(normal_node)
  #X1, y1 are the starting point, x2 and y2 are the end point.Horizon stands for the Angle between the line and the positive x axis;W stands for the Angle occupied by branches (explained in the paper)
  #Distance means distance;Node is the internal node with the start;End_num is the number of each node in edge;
  plot.data <- data.frame('x1'=0, 'y1' = 0, 'x2' = 0, 'y2' = 0, 'horizon' = pi/2,
                          'w' = 1/6*pi, 'distance' = 0, 'node' = 0,'end_num' = 0,
                          'horizon_angle' = 0)
  Set.PhyloTree <- function(sub_edge, x0, y0, w, horizon, name_list,
                            plot.data, point_list, target_node, normal_node, edge_left){
    #total vertex of tree
    vertex_total=0
    #Stores a list of angles
    angle_list=c()
    #Stores a list of vertexes
    vertex_list=c()
    vertex_sub=0
    angle=0
    for(i in 1:length(point_list)){
      #if it is sample, the number of vertices plus one
      if(point_list[i] <= length(njtree@nj$tip.label )){
        vertex_total <- vertex_total+1
        vertex_list <- append(vertex_list, 1)
      }
      #If it is an internal node, it counts several vertices connected to the internal node
      else{
        if(length(point_list[point_list>length(njtree@nj$tip.label )]) == 1){
          vertex_total <- vertex_total+length(edge_left[, 1])
          vertex_sub <- length(edge_left[, 1])
        }
        else{
          for(i1 in 1:length(edge_left[,1])){
            if(point_list[i] %in% edge_left[i1,]){
              vertex_total <- vertex_total+1
              vertex_sub <- vertex_sub+1
            } 
          }
        }
        vertex_list <- append(vertex_list,vertex_sub)
        vertex_sub <- 0
      }
    }
    #The Angle assignment is based on the proportion of the number of vertices
    for(i in 1:length(vertex_list)){
      angle_sub <- (vertex_list[i]/vertex_total)*w
      #Assign angles and put them into a list
      angle_list <- append(angle_list, angle_sub)
    }
    #Calculate the Angle with the X axis
    angle_list2 <- c()
    for(i in 1:length(angle_list)){
      if(i == 1){angle=(pi-w)/2+angle_list[i]/2+horizon-pi/2}
      else{angle <- angle+angle_list[(i-1)]/2+angle_list[i]/2+horizon-pi/2}
      angle_list2 <- append(angle_list2,angle)
    }
    #Angle adjustment (select a branch Angle to align with the previous internal node)
    if(target_node != normal_node){
      a <- which(vertex_list == max(vertex_list))
      angle_list2[a[1]] <- horizon
    }
    i=1
    while (i <= length(sub_edge[,1])) {
      if(point_list[i] > length(njtree@nj$tip.label )){
        n <- 'internal node'
      }
      #If the number is greater than the length label of labels, it is an internal node, otherwise it is a sample
      else {n <- njtree@nj$tip.label[point_list[i]]}
      #Merge the sample name into a list for coloring
      name_list <- append(name_list,n)
      #get the distance 
      distance <- as.numeric(sub_edge[i,3])
      #Calculate the coordinate
      angle <- angle_list2[i]
      x <- cos(angle)*distance+x0
      y <- sin(angle)*distance+y0
      #Merge data
      row <- list(x0, y0, x, y, angle, angle_list[i], distance, target_node, point_list[i], horizon)
      plot.data <- rbind(plot.data, row)
      i <- i+1
    }
    #Return the data of the plotting and the list of sample names
    result <- list(plot.data, name_list)
    return(result)
  }
  t=1
  while(t <= njtree@nj$Nnode){
    #The first diagram is of the internal node connected to normal
    target_node <- node_list[t]
    #he position of the target t node in the matrix (the first data point is to find the normal)
    row_list <- c()
    for(i1 in 1:length(edge[, 1])){
      if(target_node %in% edge[i1, ]){
        row_list <- append(row_list, i1)
      }
    }
    #Filter out the submatrix with target_node
    sub_edge <- edge[row_list,]
    #Edge removes the sub_edge matrix
    edge <- edge[-row_list,]
    #Extract the point that is concatenated with targe_node and put it in the list
    point_list <- c()
    for( i2 in 1:length(sub_edge[,1])){
      point <- sub_edge[i2,1:2][sub_edge[i2,1:2] != target_node]
      point_list <- append(point_list,point)
    }
    #first internal node(normal)
    if(t == 1){
      p_n<- Set.PhyloTree(sub_edge,0,0,pi*2/3,pi/2,name_list,plot.data,point_list,target_node,normal_node,edge)#º¯Êý·µ»ØµÄÁÐ±í
      plot.data <- p_n[[1]]
      name_list <- p_n[[2]]
    }
    else{
      #Find the current starting coordinate in the plot data
      row_num <- which(plot.data$end_num == target_node)
      p_n <- Set.PhyloTree(sub_edge, plot.data$x2[row_num], plot.data$y2[row_num]
                           ,plot.data$w[row_num], plot.data$horizon[row_num], name_list, plot.data, point_list,
                           target_node, normal_node, edge)
      plot.data <- p_n[[1]]
      name_list <- p_n[[2]]
    }
    node_list<-append(node_list, point_list[point_list>length(njtree@nj$tip.label)])
    t=t+1
  }
  #Add the sample_name to the plot data
  plot.data <- cbind(plot.data, sample=name_list)
  #Transform normal branches into vertical roots
  row_normal<- which(plot.data$sample == 'normal')
  plot.data$x2 [row_normal] <- 0
  plot.data$y2[row_normal] <- -plot.data$distance[row_normal]
  plot.data <- plot.data[-1,]
  #Adjust Angle of non-normal sample and internal nodes
  t=1
  while(t<= length(plot.data[, 1])){
    if(all(plot.data$end_num[t] <= length(njtree@nj$tip.label), plot.data$sample[t] != 'internal node', plot.data$sample[t] != 'normal'
           ,plot.data$horizon_angle[t] != plot.data$horizon[t]))
    {
      #the angles less than ninety degrees are adjust to thirty degrees
      if(plot.data$horizon[t]<pi/2){
        plot.data$horizon[t] <- pi/6
        angle=plot.data$horizon[t]
        plot.data$x2[t] <- plot.data$distance[t]*cos(angle)+plot.data$x1[t]
        plot.data$y2[t] <- plot.data$distance[t]*sin(angle)+plot.data$y1[t]
      }
      # above 90 are adjusted to 150
      else if(plot.data$horizon[t]>pi/2){
        plot.data$horizon[t] <- pi*5/6
        angle=plot.data$horizon[t]
        plot.data$x2[t] <- plot.data$distance[t]*cos(angle)+plot.data$x1[t]
        plot.data$y2[t] <- plot.data$distance[t]*sin(angle)+plot.data$y1[t]
      }
    }
    t=t+1
  }
  #Adjust the Angle again
  t=1
  adjust_node_list <- c()
  while(t <= length(plot.data[,1])){
    row <- ''
    row1 <- ''
    if(all(plot.data$end_num[t] <= length(njtree@nj$tip.label),plot.data$sample[t] != '',plot.data$sample[t] != 'normal'
           ,plot.data$horizon_angle[t] != plot.data$horizon[t])){
      internal_node <- plot.data$node[t]
      row1 <- which(plot.data$end_num == internal_node)
      if(length(row1)!=0){
        adjust_node <- plot.data$node[row1]
        if(adjust_node == normal_node){t <- t+1;next}
        if(adjust_node %in% adjust_node_list){t <- t+1;next}
        else{
          if(plot.data$distance[row1] < 10){
            row <- which(plot.data$node == adjust_node&plot.data$end_num <= length(njtree@nj$tip.label))
            if(length(row)!=0){
              if(plot.data$horizon[row] < pi/2){
                plot.data$horizon[row] <- plot.data$horizon[row]-pi/18
                angle=plot.data$horizon[row]
                plot.data$x2[row] <- plot.data$distance[row]*cos(angle)+plot.data$x1[row]
                plot.data$y2[row] <- plot.data$distance[row]*sin(angle)+plot.data$y1[row]
              }
              else if(plot.data$horizon[row]>pi/2){
                plot.data$horizon[row] <- plot.data$horizon[row]+pi/18
                angle <- plot.data$horizon[row]
                plot.data$x2[row] <- plot.data$distance[row]*cos(angle)+plot.data$x1[row]
                plot.data$y2[row] <- plot.data$distance[row]*sin(angle)+plot.data$y1[row]
              }
            }
          }
          adjust_node_list <- append(adjust_node_list, adjust_node)
        }
      }
    }
    t=t+1
  }
  ##Bisect the branches angle
  #right part of tree
  t=1
  sample_order_data <- plot.data[which(plot.data$sample != 'internal node' & plot.data$sample != 'normal'
                                       & plot.data$horizon <= pi/2 & plot.data$horizon > 0),]
  while(t <= length(sample_order_data[, 1])){
    u=1
    while(u <= (length(sample_order_data[, 1])-1)){
      if(sample_order_data$y1[u] == sample_order_data$y1[(u+1)]){
        if(sample_order_data$horizon[u] < sample_order_data$horizon[(u+1)]){
          u <- u+1;next
        }
        else{
          data_1 <- sample_order_data[u,]
          data_2 <- sample_order_data[(u+1),]
          sample_order_data[(u+1),] <- data_1
          sample_order_data[(u),] <- data_2
          u <- u+1; next
        }
      }
      else if(sample_order_data$y1[u]>sample_order_data$y1[(u+1)]){
        data_1 <- sample_order_data[u,]
        data_2 <- sample_order_data[(u+1),]
        sample_order_data[(u+1),] <- data_1
        sample_order_data[(u),] <- data_2
        u <- u+1;next
      }
      else{u <- u+1;next}
      u <- u+1
    }
    t <- t+1
  }
  list_right <- as.character(sample_order_data$sample)
  #left part of tree
  t <- 1
  sample_order_data <- plot.data[which(plot.data$sample!='internal node'&plot.data$sample!='normal'&plot.data$horizon>=pi/2&plot.data$horizon<pi),]
  while(t<=length(sample_order_data[,1])){
    u=1
    while(u<=(length(sample_order_data[,1])-1)){
      if(sample_order_data$y1[u] == sample_order_data$y1[(u+1)]){
        if(sample_order_data$horizon[u]<sample_order_data$horizon[(u+1)]){u=u+1;next}
        else{
          data_1 <- sample_order_data[u,]
          data_2 <- sample_order_data[(u+1),]
          sample_order_data[(u+1),] <- data_1
          sample_order_data[(u),] <- data_2
          u <- u+1;next
        }
      }
      else if(sample_order_data$y1[u] < sample_order_data$y1[(u+1)]){
        data_1 <- sample_order_data[u,]
        data_2 <- sample_order_data[(u+1),]
        sample_order_data[(u+1),] <- data_1
        sample_order_data[(u),] <- data_2
        u <- u+1;next
      }
      else{u <- u+1;next}
      u <- u+1
    }
    t <- t+1
  }
  list_left <- as.character(sample_order_data$sample)
  Bisect_branche_angle <- function(list_part, plot.data){
    t=F
    while(t == F){
      if(length(list_part) == 1|length(list_part) == 2){
        break
      }
      if(length(list_part)==3){
        first_angle <-plot.data$horizon[which(plot.data$sample==list_part[1])] 
        final_angle <- plot.data$horizon[which(plot.data$sample==list_part[length(list_part)])]
        total_part_nodes <- length(list_part)
        average_angle <- (final_angle-first_angle)/total_part_nodes
        angle <- first_angle+average_angle
        for(i in 2:(length(list_part)-1)){
          row=which(plot.data$sample==list_part[i])
          plot.data$x2[row] <- plot.data$distance[row]*cos(angle)+plot.data$x1[row]
          plot.data$y2[row] <- plot.data$distance[row]*sin(angle)+plot.data$y1[row]
          angle <- angle+average_angle
        }
        break
      }
      else{
        #Judge whether the adjustment area is on the left or right
        judge_angle <- plot.data$horizon[which(plot.data$sample == list_part[2])]
        if(judge_angle > pi/2){
          first_angle <-plot.data$horizon[which(plot.data$sample == list_part[1])] 
          final_angle <- plot.data$horizon[which(plot.data$sample == list_part[length(list_part)-2])]
          total_part_nodes <- length(list_part)
          average_angle <- (final_angle-first_angle)/(total_part_nodes-2)
          angle <- first_angle+average_angle
          for(i in 2:(length(list_part)-2)){
            row=which(plot.data$sample == list_part[i])
            plot.data$x2[row] <- plot.data$distance[row]*cos(angle)+plot.data$x1[row]
            plot.data$y2[row] <- plot.data$distance[row]*sin(angle)+plot.data$y1[row]
            angle=angle+average_angle
          }
          break
        }
        else{
          first_angle <-plot.data$horizon[which(plot.data$sample==list_part[2])] 
          final_angle <- plot.data$horizon[which(plot.data$sample==list_part[length(list_part)])]
          total_part_nodes <- length(list_part)
          average_angle <- (final_angle-first_angle)/(total_part_nodes-2)
          angle <- first_angle+average_angle
          for(i in 3:(length(list_part)-1)){
            row <- which(plot.data$sample == list_part[i])
            plot.data$x2[row] <- plot.data$distance[row]*cos(angle)+plot.data$x1[row]
            plot.data$y2[row] <- plot.data$distance[row]*sin(angle)+plot.data$y1[row]
            angle <- angle+average_angle
          }
          break
        }
      }
    }
    return(plot.data)
  }
  plot.data <- Bisect_branche_angle(list_left, plot.data)
  plot.data <- Bisect_branche_angle(list_right, plot.data)
  plot.data <- add_signature(plot.data,njtree)
  return(plot.data)
}
##add signature
add_signature <- function(plot.data, njtree){
  #add signature to plot.data
  plot.data$signature <- ''
  t <- 1
  while(t<=length(njtree@signature$branch)){
    if(length(njtree@signature$branch[[t]]) == 1){
      row <-which(plot.data$sample == njtree@signature$branch[[t]])
      plot.data$signature[row] <-as.character(njtree@signature$sig[[t]])     
      t <- t+1
      next
    }
    #find Find the row corresponding to the sample or internal in the plot.data
    row_list <- c()
    for(i in 1:length(njtree@signature$branch[[t]])){
      sample <- njtree@signature$branch[[t]][i]
      row <- which(plot.data$sample == sample)
      row_list <- append(row_list,row)
    }
    # the sample in sample_list share one common node(previous node  to the lowest row )
    # signature of normal
    if(length(njtree@signature$branch[[t]]) == length(njtree@nj$tip.label )-1){
      row <- which(plot.data$sample == 'normal')
      plot.data$signature[row] <-  as.character(njtree@signature$sig[(t)]) 
      t <- t+1
      next
    }
    #Internal nodes connected to normal
    else if(length(njtree@signature$branch[[t]]) == length(njtree@nj$tip.label )-2){
      for(i in 1:length(plot.data$sample)){
        if(plot.data$end_num[i] > length(njtree@nj$tip.label )){
          commnode_row <- i
          break
        }
      }
      plot.data$signature[commnode_row] <-  as.character(njtree@signature$sig[(t)]) 
      t=t+1
      next
    }
    # if signature does not belong to non- normal sample
    if(length(njtree@signature$branch[[t]]) != 1){
      lowest_row <- min(row_list)
      commnode_row <- which(plot.data$end_num == plot.data$node[lowest_row]) 
      if(length(commnode_row) != 0){
        plot.data$signature[commnode_row] <- as.character(njtree@signature$sig[(t)]) 
      }
    }
    else{
      row <- which(plot.data$sample == njtree@signature$branch[[t]])
      plot.data$signature[row] <-  as.character(njtree@signature$sig[(t)]) 
    }
    t=t+1
  }
  if(plot.data$signature[which(plot.data$sample == 'normal')] == ''){
    plot.data$signature[which(plot.data$sample == 'normal')] = as.character(njtree@signature$sig[1]) 
  }
  plot.data <- plot.data[order(plot.data$signature),]
  plot.data$signature <- gsub('No.Signature', 'No signature', plot.data$signature)
  plot.data$signature <- gsub('Signature.', '', plot.data$signature)
  return(plot.data)
}
##color scale set
color_set <- function(signatures){
  all_color_scale <- c("#E41A1C","#377EB8","#4DAF4A","#66C2A5","#FC8D62",
                       "#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494",
                       "#8DD3C7", "#6E016B" ,"#BEBADA", "#FB8072", "#80B1D3",
                       "#FDB462","#B3DE69","#FCCDE5","#91003F","#BC80BD",
                       "#CCEBC5" ,"#FFED6F","#1B9E77", "#D95F02", "#7570B3" ,
                       "#E7298A" ,"#66A61E" ,"#E6AB02" ,"#A6761D", "#666666",'black')
  all_signature <- append(gsub('Signature.','',row.names(signatures.cosmic)),'No signature')
  color_scale <- c()
  for(i in 1:length(signatures)){
    color <- all_color_scale[which(all_signature == signatures[i])]
    color_scale <- append(color_scale, color)
  }
  if('black' %in% color_scale){
    num <- which(color_scale == 'black')
    color_scale <- color_scale[-num]
    color_scale <- append(color_scale, 'black')
  }
  return(color_scale)
}
## plot PhyloTree 
PhyloTree <- function(plot.data, color_scale, show.mutSig){
  p <- ggplot(data = plot.data)
  text_just <- mean(as.numeric(plot.data$distance))
  if(show.mutSig){
    p <- p+geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = signature), size=1.5)
    #the color of signature.1~signature.30 is
    #   ["#E41A1C" "#377EB8" "#4DAF4A""#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494"
    #   "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" "#D9D9D9" "#BC80BD"
    #    "#CCEBC5" "#FFED6F""#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D" "#666666"]
    p <- p+scale_color_manual(values = color_scale)
  }
  else{
    p <- p+geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, color = sample),
                        data = plot.data[plot.data$sample != 'normal'&plot.data$sample != 'internal node',], 
                        size=1.5, show.legend = T)
    p <- p+geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = 'black',
                        data = plot.data[plot.data$sample == 'normal',], 
                        size = 1.5, show.legend = F )
    p <- p+geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), color = '#67001F',
                        data = plot.data[plot.data$sample == 'internal node',], 
                        size = 1.5, show.legend = F )
  }
  p <- p+theme(axis.title.x = element_blank(),
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
  p <- p+geom_text_repel(aes(x = x2*0.95, y = y2, label = sample),vjust = 1 , nudge_y = text_just/7, segment.alpha = 0,
                         data = plot.data[(plot.data$sample != 'internal node'&plot.data$sample != 'normal'),],
                         fontface = 'bold', size = 4)
  p <- p+geom_text_repel(aes(x = x2,y = y2), label = 'NORMAL', vjust = 0, nudge_y = -text_just/7, segment.alpha = 0,
                         data = plot.data[(plot.data$sample == 'normal'),], fontface = 'bold', size = 4)
  #Leaf nodes and internal nodes are distinguished by point size
  p <- p+geom_point(aes(x = x2,y = y2), data = plot.data[plot.data$sample == 'internal node',],
                    size = 1.7, color = "#67001F", fill = '#67001F', shape = 21, stroke = 0.5)
  p <- p+geom_point(aes(x = x2, y = y2), data = plot.data[plot.data$sample != 'internal node',],
                    size = 3, color = "#67001F", fill = 'white', shape = 21, stroke = 1)
  if(plot.data$distance[which(plot.data$sample == 'normal')] != 0){
      p <- p+geom_point(aes(x =0 , y = 0), size = 1.7, color = "#67001F",
                    fill = '#67001F', shape = 21, stroke = 0.5)
    }
return(p)
}
## plot.NJtree
plot.NJtree <- function(maf, use.indel = FALSE, show.mutSig = TRUE, sig.min.mut.number = 50,
                        show.heatmap = TRUE, heatmap.type = 'binary',
                        ccf.mutation.id = c("Hugo_Symbol","Chromosome","Start_Position"), 
                        ccf.mutation.sep = ":"){
  if(heatmap.type == 'binary'){
    use.ccf = FALSE
  }
  else{
    use.ccf = TRUE
  }
  # set NJtree object(njtree)
  njtree <- NJtree(maf, use.indel, use.ccf, ccf.mutation.id = ccf.mutation.id, ccf.mutation.sep = ccf.mutation.sep)
  # PhyloTree input data
  PhyloTree.input_data<- PhyloTree.input(njtree)
  PhyloTree.input_data <- PhyloTree.input_data[(PhyloTree.input_data$distance!=0|
                                                  PhyloTree.input_data$sample == 'normal'),]
  #Set the color
  color_scale <- color_set(unique(PhyloTree.input_data$signature))
  tree <- PhyloTree(PhyloTree.input_data, color_scale ,show.mutSig)
  heatmap <- mut.heatmap(maf, use.indel, use.ccf,ccf.mutation.id = ccf.mutation.id, ccf.mutation.sep = ccf.mutation.sep)
  maf@patientID <- paste(maf@patientID, ".NJtree", sep = "")
  if(show.mutSig){
    maf@patientID <- paste(maf@patientID, ".mutsig", sep = "")
  }
  if(show.heatmap){
    plot.njtree <- ggdraw() + draw_plot(tree, x = 0,y = 0, width = 0.8) + draw_plot(heatmap, x = 0.8,y = -0.035, width = 0.2) 
    if(!use.ccf){
      if(use.indel){
        ggsave(filename = paste("./Figures/", maf@patientID, ".useindel.pdf", sep = ""),
               plot = plot.njtree, width = 14, height = 7)}
      else{
        ggsave(filename = paste("./Figures/", maf@patientID, ".pdf", sep = ""),
               plot = plot.njtree, width = 14, height = 7)}
      }
      else{
        if(use.indel){
          ggsave(filename = paste("./Figures/", maf@patientID, ".useindel.ccf.pdf", 
                  sep = ""),plot = plot.njtree, width = 14, height = 7)}
        else{ggsave(filename = paste("./Figures/", maf@patientID, ".ccf.pdf", sep = ""),
                 plot = plot.njtree, width = 14, height = 7) }

     }
    return(plot.njtree)
  }
  else{
    return(tree)
  }
  
}



