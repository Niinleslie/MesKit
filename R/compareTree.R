#' compareTree
#' @description compares two phylogenetic trees and returns a detailed report of several distance methods
#' 
#' 
#' @param tree1 a single class PhyloTree
#' @param tree2 a single class PhyloTree
#' @param plot a logical value, default FALSE. If TRUE, 
#' the two trees will be plotted on the same device and their similarities will be shown.
#' @param min.ratio double (Default:1/30). If min.ratio is not NULL,
#'  all edge length of a phylogenetic tree should be greater than
#'  min.ratio*the longest edge length.
#'  If not, the edge length will be reset as min.ratio*longest edge length.
#' 
#' @param common.lty line type of common branch.
#' 
#' @return returns a vector containing the following tree distance methods by R package phangorn
#' Symmetric.difference  Robinson-Foulds distance
#' KF-branch distance  the branch score distance (Kuhner & Felsenstein 1994)
#' Path.difference  difference in the path length, counted as the number of branches 
#' Weighted.path.difference	 difference in the path length, counted using branches lengths
#' 
#' @examples
#' tree1 <- getPhyloTree(maf, method = "NJ")[['patientID']]
#' tree2 <- getPhyloTree(maf, method = "MP")[['patientID']]
#' compareTree(tree1, tree2)
#' compareTree(tree1, tree2, plot = TRUE)
#' @export compareTree

compareTree <- function(tree1,
                        tree2,
                        plot = FALSE,
                        branchCol = "mutSig",
                        show.heatmap = FALSE,
                        min.ratio = 1/30,
                        show.heatmap = FALSE, 
                        use.ccf = FALSE,
                        show.bootstrap = FALSE,
                        use.box = FALSE,
                        common.lty = "solid"){
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
	    treedat1 <- getTreeData(tree1, show.mutSig = T)
	    treedat2 <- getTreeData(tree2, show.mutSig = T)
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
	    p1 <- drawPhyloTree(phyloTree = tree1,
	                        treeData = treedat1,
	                        branchCol = branchCol,
	                        show.heatmap = show.heatmap,
	                        show.bootstrap = show.bootstrap,
	                        use.box = use.box,
	                        use.ccf = use.ccf,
	                        compare = compare,
	                        common.lty = common.lty,
	                        min.ratio = min.ratio)
	    p2 <- drawPhyloTree(phyloTree = tree2,
	                        treeData = treedat2,
	                        branchCol = branchCol,
	                        show.bootstrap = show.bootstrap,
	                        use.box = use.box,
	                        show.heatmap = show.heatmap,
	                        use.ccf = use.ccf
	                        compare = compare,
	                        common.lty = common.lty,
	                        min.ratio = min.ratio)
	    ptree <- cowplot::plot_grid(p1,
	                                p2,
	                                labels = c(tree1@method,tree2@method))
	    # p <- ggpubr::ggarrange(p1, p2, nrow =1, common.legend = TRUE, legend="top",labels = c(tree1@method,tree2@method))
	    return(ptree)
	}
    
	return(dist)
}
