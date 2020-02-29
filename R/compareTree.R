#' compareTree
#' @description compares two phylogenetic trees and returns a detailed report of several distance methods
#' 
#' 
#' @param tree1 a phylogenetic tree
#' @param tree2 a phylogenetic tree
#' @param plot a logical value. If TRUE, the two trees are plotted on the same device and their similarities are shown.
#' @param min.ratio a numeric value.If not NULL, min.length = max(tree$edge.length)*min.ratio.
#'                 Length of branches which are less than min.length will be adjusted to min.length.
#' @param compare.linetype line type of common branch.
#' @return returns a vector containing the following tree distance methods by R package phangorn
#' Symmetric.difference  Robinson-Foulds distance
#' KF-branch distance  the branch score distance ( (Kuhner & Felsenstein 1994))
#' Path.difference  difference in the path length, counted as the number of branches 
#' Weighted.path.difference	 difference in the path length, counted using branches lengths
#' @examples
#' tree1 <- getPhyloTree(maf, method = "NJ")
#' tree2 <- getPhyloTree(maf, method = "MP")
#' compareTree(tree1, tree2)
#' compareTree(tree1, tree2, plot = TRUE)
#' @export compareTree


compareTree <- function(tree1, tree2, plot = FALSE, min.ratio = 1/30, compare.linetype = "solid"){
	phylo.tree1 <- tree1@tree
	phylo.tree2 <- tree2@tree

	dist <- phangorn::treedist(phylo.tree1, phylo.tree2)
	names(dist) <- c("Symmetric.difference", "KF-branch distance", "Path difference", "Weighted path difference")
	if(plot){
	    if(!is.null(min.ratio)){
	        min1 <- max(tree1@tree$edge.length)*min.ratio
	        min2 <- max(tree2@tree$edge.length)*min.ratio
	        tree1@tree$edge.length[tree1@tree$edge.length < min1] <- min1
	        tree2@tree$edge.length[tree2@tree$edge.length < min2] <- min2
	    }
	    treedat1 <- getTreeDat(tree1, show.mutSig = T)
	    treedat2 <- getTreeDat(tree2, show.mutSig = T)
	    m12 <- match(treedat1[sample == "internal node",]$label, treedat2[sample == "internal node",]$label)
	    if(length(m12[!is.na(m12)]) > 0){
	        cat(paste0("Both tree have ",length(m12[!is.na(m12)]), " same branch"))
	        compare <- TRUE
	        treedat1$is.match <- 'NO'
	        treedat2$is.match <- 'NO'
	        x <- 1
	        for(i in 1:length(m12)){
	            if(is.na(m12[i])){
	                next
	            }
	            else{
	                pos1 <- which(treedat1$end_num == treedat1[sample == "internal node",]$end_num[i])
	                pos2 <- which(treedat2$end_num == treedat2[sample == "internal node",]$end_num[m12[i]])
	                treedat1$is.match[pos1] <- paste0("Common ",x)
	                treedat2$is.match[pos2] <- paste0("Common ",x)
	                x <- x + 1
	            }
	        }
	    }else{
	        cat("Both tree have not the same branch")
	        compare <- FALSE
	    }
	    p1 <- drawPhyloTree(treeData = treedat1, myBoots = tree1@bootstrap.value,
	                        patientID = tree1@patientID, show.mutSig = TRUE,
	                        show.bootstrap = FALSE, use.box = FALSE,
	                        show.heatmap = FALSE, use.ccf = FALSE, compare = compare, compare.linetype = compare.linetype)
	    p2 <- drawPhyloTree(treeData = treedat2, myBoots = tree2@bootstrap.value,
	                        patientID = tree2@patientID, show.mutSig = TRUE,
	                        show.bootstrap = FALSE, use.box = FALSE,
	                        show.heatmap = FALSE, use.ccf = FALSE, compare = compare, compare.linetype = compare.linetype)
	    p <- ggpubr::ggarrange(p1, p2, nrow =1, common.legend = TRUE, legend="top",labels = c(tree1@method,tree2@method))
	    return(p)
	}
    
	return(dist)
}
