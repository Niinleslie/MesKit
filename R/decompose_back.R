(S:633.375,TP2-2:718.625,((TC1:1028.416667,normal:332.5833333):45.5,((TP1:735.6,TC2:713.4):33.5,V:611.5):108.75):36);


mid_node <- function(edges, tip_num, edge){
	edges <- matrix(edges)
	N1 <- edge[1]
	N2 <- edge[1]

	N1_related <- as.vector(apply(edges,1,function(x) any(x == N1)),))
	N1_related <- N1_related[node_related != N1 & node_related != N2]

	N2_related <- as.vector(apply(edges,1,function(x) any(x == node_id)),))
	N2_related <- N2_related[node_related != N1 & node_related != N2]
		
	if(N1_related > tip_num & N2_related > tip_num){
		return("BOTH_MID")
	}
	else if(N1_related > tip_num){
		return(list(N2,N2_related))
	}
	else{return list(N1,N1_related)}
}



read.NJtree <- function(newick){
	edges <- mat.nj$edge
	tips <- mat.nj$tip.label
	mid_nodes <- seq(length(tips)+1,max(edge))


	nodes <- vector(length = dim(edges)[1])
	nodes[1:length(tips)] <- tips
	length(nodes) <- dim(edges)[1]

	for(i in dim(edges)[1]){
		edge <- edge[i,]
		min_node <- min(edge)
		max_node <- max(edge)

		if(min_node %in% 1:length(tips)){
			edges_label[i] <- tips[min_node]
		}
		else if(mid_node(edges, tip_num, edge) == "BOTH_MID"){			
		}
		else{
			related_leaf_nodes <- mid_leaf_node(edges, tip_num, edge)[[2]]
			edges_label[i] <- paste(related_leaf_nodes, collapse = "∩")

		}
		}


			}
		}

	}


	for (i in 1:dim(edges)[1]){
		mid_node <- max(edges[i,])
		min_node <- min(edges[i,])
		if (is.null(nodes[mid_node]) == FALSE){
			nodes[mid_node] <- tips[min_node]
		}
		else if(min_node %in% 1:length(tips)){
			nodes[mid_node] <- paste(nodes[mid_node],tips[min_node],sep="∩")
		}
	}
}
