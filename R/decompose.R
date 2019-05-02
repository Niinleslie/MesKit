#' Get sampleID of private branchs and shared branches/trunk in NJtree
#' @param mat.nj a nj tree object generated from ape
#' @param patientID represent the ID of the specific patient
#' @return a dataframe with sampleID of each branch and trunk
#' @export
#'
#' @examples
#  NJtree_sampleID <- read.NJtree(mat.nj)

read.NJtree <- function(mat.nj, patientID){
	edges <- mat.nj$edge
	tips <- mat.nj$tip.label
	leaf_nodes <- seq(1, length(tips))
	mid_nodes <- seq(length(tips)+1, max(edges))
	root <- which(tips == "normal")

	#return the adjacent nodes of the specific node
	get_adjacent_node <- function(node){

		N_adjacent <- as.vector(edges[which(edges[,1] == node | edges[,2] == node),])
		N_adjacent_mid <- N_adjacent[N_adjacent != node & !(N_adjacent %in% leaf_nodes)]

		if(!is.null(N_adjacent_mid)){
			return(N_adjacent_mid)

		}
	}


	root_R <- get_adjacent_node(root)
	dist <- list()
	dist[1] <- root_R
	Node_dist <- rep(0, times = max(edges))
	Node_dist[root_R] <- 1
	Node_after <- list(NULL)
	length(Node_after) <- max(edges)
	for (i in 2:mat.nj$Nnode){
		i_B <- i-1
		if(length(unlist(dist)) < mat.nj$Nnode){
			for(N in dist[[i_B]]){
				R_Nodes <- get_adjacent_node(N)
				R_Nodes <- R_Nodes[!(R_Nodes %in% unlist(dist))]
				dist[[i]] <- R_Nodes
				Node_dist[R_Nodes] <- i
				if(!is.null(R_Nodes)){
					Node_after[[N]] <- R_Nodes
				}
			}
		}
		else{break}
	}

	#return leaf nodes which on the same path of the specific node
	get_leaf_node <- function(nodes){

		leafs=c()
		for(node in nodes){
			N_adjacent <- as.vector(edges[which(edges[,1] == node | edges[,2] == node),])
			N_adjacent <- N_adjacent[N_adjacent != node]
			N_adjacent_after <- c()
			if(length(intersect(N_adjacent, length(tips)+1:max(edges)))){
			N_adjacent_after <- N_adjacent[!(N_adjacent %in% leaf_nodes) & N_adjacent %in% Node_after[[node]]]
			}
			N_adjacent_leaf <- intersect(N_adjacent, leaf_nodes)

			if(length(N_adjacent_leaf)){
				leafs <- unique(c(leafs, N_adjacent_leaf))
			}

			if(length(N_adjacent_after)){
				for(mid_N in N_adjacent_after){
					leafs <- c(leafs, get_leaf_node(mid_N))
				}
			}

			return(leafs)
		}
	}


	edges_label <- c()
	for (i in 1:dim(edges)[1]){
		edge <- edges[i,]
		if(root %in% edge){
			edges_label <- c(edges_label, paste(tips[-root], collapse = "∩"))
		}
		else if(min(edge) > length(tips)){
			N_end <- edge[2]
			if(Node_dist[edge[1]] > Node_dist[N_end]){
				N_end <- edge[1]
			}
			label <- paste(tips[get_leaf_node(N_end)], collapse = "∩")
			edges_label <- c(edges_label, label)
		}
		else{
			edges_label <- c(edges_label, tips[min(edge)])
		}
	}

	NJtree_sampleID <- edges_label[order(unlist(lapply(edges_label, function(x) length(unlist(strsplit(x, split = "∩"))))))]
	write.table(NJtree_sampleID, file = paste(patientID, ".NJtree.edges", sep = ""), quote = F, row.names = F, col.names = F, fileEncoding = "UTF-8")
}

