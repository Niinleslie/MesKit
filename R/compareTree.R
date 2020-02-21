#' compareTree
#' @description compares two phylogenetic trees and returns a detailed report of several distance methods
#' 
#' 
#' @param tree1 a phylogenetic tree
#' @param tree2 a phylogenetic tree
#' @return returns a vector containing the following tree distance methods by R package phangorn
#' Symmetric.difference  Robinson-Foulds distance
#' KF-branch distance  the branch score distance ( (Kuhner & Felsenstein 1994))
#' Path.difference  difference in the path length, counted as the number of branches 
#' Weighted.path.difference	 difference in the path length, counted using branches lengths
#' @examples
#' tree1 <- getPhyloTree(maf, method = "NJ")
#' tree2 <- getPhyloTree(maf, method = "MP")
#' compareTree(tree1, tree2)
#' @export compareTree


compareTree <- function(tree1, tree2){
	phylo.tree1 <- tree1@tree
	phylo.tree2 <- tree2@tree

	dist <- phangorn::treedist(phylo.tree1, phylo.tree2)
	names(dist) <- c("Symmetric.difference", "KF-branch distance", "Path difference", "Weighted path difference")

	return(dist)
}
