#' mutTrunkBranch
#' @description Summarize and conduct paired Wilcoxon test of mutations of trunk/branches in a phylogenetic tree.
#' 
#' @param tree.mutSig the output of treeMutSig function.
#' @param conf.level confidence level of the interval for wilcox.test. Default: 0.95. Option: on the scale of 0 to 1.
#' 
#' @examples
#' mutTrunkBranch(tree.mutSig, conf.level = 0.95)
#' @return Box plots based on mutational categories
#' @export mutTrunkBranch

mutTrunkBranch <- function(tree.mutSig, conf.level = 0.95) {
    mutTB.list <- suppressWarnings(lapply(tree.mutSig, doMutTrunkBranch,
                                   conf.level = conf.level)) 
    return(mutTB.list)
}

doMutTrunkBranch <- function(tree.mutSig, conf.level = 0.95){
    ## input data from tree.mutSig
    ls.BT <- .dataProcessBT(tree.mutSig)
    df.pValue <- ls.BT$df.pValue
    sigsInputBoxplot <- ls.BT$sigsInputBoxplot
    ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    
    ## generate output data.frame with quantity of mutations in different categories
    output <- data.frame(matrix(nrow=6, ncol=2))
    colnames(output) <- c("Trunk", "Branch")
    output <- cbind(Group=ls.mutationGroup, output)
    output <- merge(output, df.pValue, by=c("Group"))
    output <- cbind(output, Significance=rep("-", nrow(output)))
    for (mutationGroup in ls.mutationGroup) {
        output$Branch[which(output$Group == mutationGroup)] <- sum(sigsInputBoxplot[which(
            sigsInputBoxplot$Group == mutationGroup & sigsInputBoxplot$BT == "Branch"),]$mut.num)
        output$Trunk[which(output$Group == mutationGroup)] <- sum(sigsInputBoxplot[which(
            sigsInputBoxplot$Group == mutationGroup & sigsInputBoxplot$BT == "Trunk"),]$mut.num)
        if (!is.null(output[which(output$p.value < (1-conf.level)), ]$Significance)) {
            output[which(output$p.value < (1-conf.level)), c("Significance")] <- "*"
        }
    }
    return(output)
}

.dataProcessBT <- function(tree.mutSig) {
    ## input data from tree.mutSig
    sigsInput <- tree.mutSig$sigsInput
    mutSigsOutput <- tree.mutSig$mutSigsOutput
    
    ## label the Trunk
    if (any(mutSigsOutput$alias == "T")){
        trunkName <- mutSigsOutput[which(mutSigsOutput$alias == "T"), ]$branch
    } else {
        stop("Trunk ERROR: There is no trunk mutation and the branch-trunk plot could not be plotted.
             Warnings and outputs from function getNJtree should be checked.")
    } 
    
    ## separate trunk and branch data
    sigsInput.trunk <- sigsInput[which(rownames(sigsInput) == trunkName), ]
    sigsInput.branch <- sigsInput[which(rownames(sigsInput) != trunkName), ]
    sigsInput.branch <- colSums(sigsInput.branch)
    sigsInputBT <- rbind(Trunk=sigsInput.trunk, Branch=sigsInput.branch)
    sigsInputBTTrans <- data.frame(Mutational_Type=colnames(sigsInputBT), t(sigsInputBT))
    ls.mutationGroup <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    
    ## generate Mutation Type for every column
    for (mutationGroup in ls.mutationGroup) {
        sigsInputBTTrans$Group[which(grepl(mutationGroup, sigsInputBTTrans$Mutational_Type))] <- mutationGroup
    }
    
    sigsInputBSum <- sigsInputBTTrans %>% dplyr::group_by(Group) %>% dplyr::summarise(sum = sum(Branch))
    sigsInputTSum <- sigsInputBTTrans %>% dplyr::group_by(Group) %>% dplyr::summarise(sum = sum(Trunk))
    
    sigsInputBTTrans <- cbind(sigsInputBTTrans, 
                              BranchFrac=rep(0, nrow(sigsInputBTTrans)), 
                              TrunkFrac=rep(0, nrow(sigsInputBTTrans)))
    for (mutationGroup in ls.mutationGroup) {
        groupBSum <- sigsInputBSum$sum[which(sigsInputBSum$Group == mutationGroup)]
        if (groupBSum == 0) {
            sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Branch <- 0
        }
        groupTSum <- sigsInputTSum$sum[which(sigsInputTSum$Group == mutationGroup)]
        if (groupBSum == 0) {
            sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Trunk <- 0
        }
        
        sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$BranchFrac <- 
            100*sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Branch/groupBSum
        sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$TrunkFrac <- 
            100*sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]$Trunk/groupTSum
    }
    
    sigsInputBoxplot <- data.frame(matrix(nrow=0, ncol=5))
    colnames(sigsInputBoxplot) <- c("GroupBT", "Group", "BT", "mut.frac", "mut.num")
    for (mutationGroup in ls.mutationGroup) {
        dat.group <- sigsInputBTTrans[which(sigsInputBTTrans$Group == mutationGroup), ]
        df.groupB <- data.frame(rep(paste(mutationGroup, "Branch", sep=" "), nrow(dat.group)), 
                                rep(mutationGroup, nrow(dat.group)), 
                                rep("Branch", nrow(dat.group)), 
                                dat.group$BranchFrac, 
                                dat.group$Branch)
        df.groupT <- data.frame(rep(paste(mutationGroup, "Trunk", sep=" "), nrow(dat.group)), 
                                rep(mutationGroup, nrow(dat.group)), 
                                rep("Trunk", nrow(dat.group)),
                                dat.group$TrunkFrac, 
                                dat.group$Trunk)
        colnames(df.groupB) <- c("GroupBT", "Group", "BT", "mut.frac", "mut.num")
        colnames(df.groupT) <- c("GroupBT", "Group", "BT", "mut.frac", "mut.num")
        sigsInputBoxplot <- rbind(sigsInputBoxplot, df.groupB, df.groupT)
    }
    
    df.pValue <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(df.pValue) <- c("Group", "p.value")
    for (mutationGroup in ls.mutationGroup) {
        pValue <- wilcox.test(
            sigsInputBoxplot[
                which(sigsInputBoxplot$Group == mutationGroup & 
                          sigsInputBoxplot$BT == "Branch"), ]$mut.frac, 
            sigsInputBoxplot[
                which(sigsInputBoxplot$Group == mutationGroup & 
                          sigsInputBoxplot$BT == "Trunk"), ]$mut.frac, 
            paired=TRUE, alternative = "two.sided", conf.level=conf.level, 
            exact=FALSE
        )$p.value
        row.pValue <- data.frame(mutationGroup, pValue)
        colnames(row.pValue) <- c("Group", "p.value")
        df.pValue <- rbind(df.pValue, row.pValue)
    }
    output <- list(df.pValue=df.pValue, sigsInputBoxplot=sigsInputBoxplot)
    return(output)
}