#' Get sampleID of private branchs and shared branches/trunk in NJtree
#' @param matTree a nj tree object generated from ape
#' @return a dataframe with sampleID of each branch and trunk
#'
#' @examples
#  branches <- readPhyloTree(matTree)

readPhyloTree <- function(matTree){
  # edges <- matTree$edge
  # tips <- matTree$tip.label
  # leaf_nodes <- seq(1, length(tips))
  # mid_nodes <- seq(length(tips)+1, max(edges))
  # root <- which(tips == "NORMAL")
  # 
  # #return the adjacent nodes of the specific node
  # get_adjacent_node <- function(node){
  #   
  #   N_adjacent <- as.vector(edges[which(edges[,1] == node | edges[,2] == node),])
  #   N_adjacent_mid <- N_adjacent[N_adjacent != node & !(N_adjacent %in% leaf_nodes)]
  #   
  #   if(!is.null(N_adjacent_mid)){
  #     return(N_adjacent_mid)
  #     
  #   }
  # }
  # 
  # root_R <- get_adjacent_node(root)
  # dist <- list()
  # dist[1] <- root_R
  # Node_dist <- rep(0, times = max(edges))
  # Node_dist[root_R] <- 1
  # Node_after <- list(NULL)
  # length(Node_after) <- max(edges)
  # for (i in 2:matTree$Nnode){
  #   i_B <- i-1
  #   if(length(unlist(dist)) < matTree$Nnode){
  #     for(N in dist[[i_B]]){
  #       R_Nodes <- get_adjacent_node(N)
  #       R_Nodes <- R_Nodes[!(R_Nodes %in% unlist(dist))]
  #       dist[[i]] <- R_Nodes
  #       Node_dist[R_Nodes] <- i
  #       if(!is.null(R_Nodes)){
  #         Node_after[[N]] <- R_Nodes
  #       }
  #     }
  #   }
  #   else{break}
  # }
  # 
  # #return leaf nodes which on the same path of the specific node
  # get_leaf_node <- function(nodes){
  #   
  #   leafs=c()
  #   for(node in nodes){
  #     N_adjacent <- as.vector(edges[which(edges[,1] == node | edges[,2] == node),])
  #     N_adjacent <- N_adjacent[N_adjacent != node]
  #     N_adjacent_after <- c()
  #     if(length(intersect(N_adjacent, length(tips)+1:max(edges)))){
  #       N_adjacent_after <- N_adjacent[!(N_adjacent %in% leaf_nodes) & N_adjacent %in% Node_after[[node]]]
  #     }
  #     N_adjacent_leaf <- intersect(N_adjacent, leaf_nodes)
  #     
  #     if(length(N_adjacent_leaf)){
  #       leafs <- unique(c(leafs, N_adjacent_leaf))
  #     }
  #     
  #     if(length(N_adjacent_after)){
  #       for(mid_N in N_adjacent_after){
  #         leafs <- c(leafs, get_leaf_node(mid_N))
  #       }
  #     }
  #     
  #     return(leafs)
  #   }
  # }
  # 
  # 
  # edges_label <- c()
  # for (i in 1:dim(edges)[1]){
  #   edge <- edges[i,]
  #   if(root %in% edge){
  #     edges_label <- c(edges_label, paste(tips[-root], collapse = "∩"))
  #   }
  #   else if(min(edge) > length(tips)){
  #     N_end <- edge[2]
  #     if(Node_dist[edge[1]] > Node_dist[N_end]){
  #       N_end <- edge[1]
  #     }
  #     label <- paste(tips[get_leaf_node(N_end)], collapse = "∩")
  #     edges_label <- c(edges_label, label)
  #   }
  #   else{
  #     edges_label <- c(edges_label, tips[min(edge)])
  #   }
  # }
  Root <- which(matTree$tip.label == "NORMAL")
  matTree <- root(matTree, Root)
  result <- list()
  for(i in 1:(length(matTree$edge.length)+1)){
      result[[i]] <- NA
  }
  end <- matTree$edge[which(matTree$edge[,2] == Root),1]
  for(i in 1:length(matTree$tip.label)){
      row <- matTree$edge[which(matTree$edge[,2] == i), ]
      if(i == Root){
          next
      }
      result[[i]][2] <- i
      while (TRUE) {
          node <- row[1]
          if(node == end){
              result[[node]] <- append(result[[node]], i)
              break
          }
          else{
              result[[node]] <- append(result[[node]], i)
              row <- matTree$edge[which(matTree$edge[,2] == node),]
          }
      }
  }
  g <- lapply(result, function(x){return(paste(sort(matTree$tip.label[x[-1]], decreasing = T),collapse = "∩"))})
  edges_label <- unlist(g)[-Root]

  nodes.num <- seq_along(matTree$tip.label)-1
  NJtree_branch <- edges_label[order(unlist(lapply(edges_label, function(x) length(unlist(strsplit(x, split = "∩"))))))]
  edges.num <- length(strsplit(tail(NJtree_branch,1), split = "∩")[[1]])
  NJtree_branch_Alias <- c(paste(rep("B", length(NJtree_branch)-1), (length(NJtree_branch)-1):1, sep = ""), "T")
  NJtree_sampleID <- data.frame(Branch = NJtree_branch, Alias = NJtree_branch_Alias)
  
  return(NJtree_sampleID)
}