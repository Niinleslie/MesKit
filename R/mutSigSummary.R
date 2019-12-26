#' mutSigSummary
#' @description Provide a summary data frame for signatures in different branches.
#' 
#' @param treeMSOutput The output of function treeMutationalSig.
#' 
#' @return Data frame of each set/branch's mutational signature.
#' 
#' @examples
#' mutSig.summary <- mutSigSummary(treeMSOutput)
#' 
#' @export mutSigSummary

mutSigSummary <- function(treeMSOutput){
    mutSigsOutput <- treeMSOutput$mutSigsOutput
    df.aetiology <- treeMSOutput$df.aetiology
    
    ## list of branch names
    ls.branchesName <- mutSigsOutput$branch
    ls.aeti <- c()
    
    for (branch in ls.branchesName) {
        signature <- as.character(mutSigsOutput[which(mutSigsOutput$branch == branch), ]$sig)
        aetiology <- as.character(df.aetiology[which(df.aetiology$sig == signature), ]$aeti)
        ls.aeti <- c(ls.aeti, aetiology)
    }
    
    ## rearrange the order of columns
    if (is.null(mutSigsOutput$putative_driver_genes)) {
        mutSigsOutput <- data.frame(branch=mutSigsOutput$branch,
                                    alias=mutSigsOutput$alias, 
                                    mut.num=mutSigsOutput$mut.num, 
                                    sig=mutSigsOutput$sig, 
                                    sig.prob=mutSigsOutput$sig.prob, 
                                    aeti=ls.aeti)
        colnames(mutSigsOutput) <- c("Branch", "Alias", "Mutation quantity", "Signature", "Signature weight", "Aetiology")
    } else {
        mutSigsOutput <- data.frame(branch=mutSigsOutput$branch,
                                    alias=mutSigsOutput$alias, 
                                    mut.num=mutSigsOutput$mut.num, 
                                    sig=mutSigsOutput$sig, 
                                    sig.prob=mutSigsOutput$sig.prob,
                                    aeti=ls.aeti,
                                    putative_driver_genes=mutSigsOutput$putative_driver_genes)
        colnames(mutSigsOutput) <- c("Branch", "Alias", "Mutation quantity", "Signature", "Signature weight", "Aetiology", "Oncogene list")
    }
    
    return(mutSigsOutput)
}