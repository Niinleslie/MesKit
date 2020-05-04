#' compareTree
#' @description compares two phylogenetic trees and returns a detailed report of several distance methods
#' 
#' 
#' @param phyloTree1 a phyloTree object
#' @param phyloTree2 a phyloTree object
#' @param plot FALSE(Default). If TRUE, two trees will be plotted on the same device and their similarities will be shown.
#' @param min.ratio double (Default: 1/30). If min.ratio is not NULL,
#' all edge length which are smaller than min.ratio*the longest edge length will be reset as min.ratio*longest edge length. 
#' @param common.lty line type of common branch.
#' 
#' @return a vector containing the following tree distance methods by R package phangorn
#' Symmetric.difference  Robinson-Foulds distance
#' KF-branch distance  the branch score distance (Kuhner & Felsenstein 1994)
#' Path.difference  difference in the path length, counted as the number of branches 
#' Weighted.path.difference	 difference in the path length, counted using branches lengths
#' 
#' @examples
#' phyloTree1 <- getPhyloTree(maf, method = "NJ")[['patientID']]
#' phyloTree2 <- getPhyloTree(maf, method = "MP")[['patientID']]
#' compareTree(phyloTree1, phyloTree2)
#' compareTree(phyloTree1, phyloTree2, plot = TRUE)
#' @export compareTree

compareTree <- function(phyloTree1,
                        phyloTree2,
                        plot = FALSE,
                        min.ratio = 1/30,
                        show.bootstrap = FALSE,
                        use.box = FALSE,
                        common.lty = "solid"){
    
    if(min.ratio <= 0){
        stop("Error: min.ratio must greater than 0")
    }
    
	tree1 <- phyloTree1@tree
	tree2 <- phyloTree2@tree
	dist <- phangorn::treedist(tree1, tree2)
	names(dist) <- c("Symmetric.difference", "KF-branch distance", "Path difference", "Weighted path difference")
	if(plot){
	    if(!is.null(min.ratio)){
	        min1 <- max(phyloTree1@tree$edge.length)*min.ratio
	        min2 <- max(phyloTree2@tree$edge.length)*min.ratio
	        phyloTree1@tree$edge.length[phyloTree1@tree$edge.length < min1] <- min1
	        phyloTree2@tree$edge.length[phyloTree2@tree$edge.length < min2] <- min2
	    }
	    treedat1 <- getTreeData(phyloTree1, compare = TRUE)
	    treedat2 <- getTreeData(phyloTree2, compare = TRUE)
	    m12 <- match(treedat1[sample == "internal node",]$label, treedat2[sample == "internal node",]$label)
	    if(length(m12[!is.na(m12)]) > 0){
	        cat(paste0("Both tree have ",length(m12[!is.na(m12)]), " same branches"))
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
	        return(dist)
	        compare <- FALSE
	    }
	    p1 <- drawPhyloTree(phyloTree = phyloTree1,
	                        treeData = treedat1,
	                        show.bootstrap = show.bootstrap,
	                        use.box = use.box,
	                        compare = compare,
	                        common.lty = common.lty,
	                        min.ratio = min.ratio)
	    p2 <- drawPhyloTree(phyloTree = phyloTree2,
	                        treeData = treedat2,
	                        show.bootstrap = show.bootstrap,
	                        use.box = use.box,
	                        compare = compare,
	                        common.lty = common.lty,
	                        min.ratio = min.ratio)
	    ptree <- cowplot::plot_grid(p1,
	                                p2,
	                                labels = c(phyloTree1@method,phyloTree2@method))
	    # p <- ggpubr::ggarrange(p1, p2, nrow =1, common.legend = TRUE, legend="top",labels = c(phyloTree1@method,phyloTree2@method))
	    return(list(compare.dist = dist, compare.plot = ptree))
	}
    
	return(dist)
}
