#' Check Mutational Signature for each branch of phylogenetic tree
#' @description Read maf file as data.frame. Define branches' set relationship
#'  by re-labeling their tumor sample barcode from the smallest set. Calcualte 
#'  each branch's mutational signature weight according to cosmic reference and
#'  pick the maxium. Return a data frame of each set/branch's mutational 
#'  signature.
#' 
#' @import reshape2 BSgenome BSgenome.Hsapiens.UCSC.hg19 GenomeInfoDb grDevices
#' graphics utils deconstructSigs
#' @importFrom plyr match_df
#' 
#' @param maf_file specify a maf document/directory as the input of the 
#' function
#' @param branch_file specify a txt document/directory as the input of 
#' the branches (needed to be refined)
#' @return data frame of each set/branch's mutational signature.
#' 
#' @export treeMutationalSig
#'
#' @examples
#' \dontrun{
#' treeMutationalSig(maf_file, branch_file)
#' treeMutationalSig(maf_file, branch_file, driver_genes_dir)
#' treeMutationalSig(maf_file, branch_file, 
#' driver_genes_dir, mut.threshold=30)
#'}

## main function
treeMutationalSig <- function(mafData, branch, 
                                 patientID, refBuild, 
                                 driverGenesDir=FALSE, mutThreshold=50){
    maf_input <- mafData
    ## get mutationalSigs-related  infomation
    datSample <- data.frame(as.character(maf_input$Tumor_Sample_Barcode), 
                             stringsAsFactors=FALSE)
    datChr <- data.frame(chr=as.character(maf_input$Chromosome), 
                          stringsAsFactors=FALSE)
    datChr$chr <- paste("chr", datChr$chr, sep="")
    datPosStart <- maf_input$Start_Position
    datPosEnd <- maf_input$End_Position
    datRef <- maf_input$Reference_Allele
    datAlt <- maf_input$Tumor_Seq_Allele2
    datNum <- seq_along(datAlt)
    datMutgene <-  maf_input$Hugo_Symbol
    mutId <- select(tidyr::unite(maf_input, "mut.id", 
                                  Hugo_Symbol, Chromosome, 
                                  Start_Position, End_Position, 
                                  Reference_Allele, Tumor_Seq_Allele2, 
                                  sep=":"), mut.id)
    mutSigRef <- data.frame(datNum, datSample, 
                              datChr, datPosStart, 
                              datPosEnd, datRef, 
                              datAlt, datMutgene, 
                              mutId)
    colnames(mutSigRef) <- c("ID", "Sample", 
                               "chr", "pos", 
                               "pos_end", "ref", 
                               "alt", "Hugo_Symbol", 
                               "mut_id")
    
    ## get branch infomation
    branch <- branch[order(nchar(branch), branch)]
    branches <- strsplit(branch, split='∩')
    
    
    ## output collection
    mutSigsOutput <- data.frame()
    mutBranches <- data.frame()
    mutBranchesOutput <- list()
    listBranchName <- c()
    ## generate mutational signautres for different branches
    for (branchCounter in length(branches):1){
        ## generate a single branch
        branch <- Filter(Negate(is.na), 
                         branches[[branchCounter]])
        mutBranch <- mutSigRef[which(mutSigRef$Sample %in% branch), ]
        
        for (tsb in branch){
            ## generate the intersection of the branch
            mutTsb <- mutSigRef[which(mutSigRef$Sample %in% tsb), ]
            mutBranch <- match_df(mutBranch, mutTsb, on=c("chr", "pos", 
                                                             "pos_end", "ref", 
                                                             "alt"))
        }
        
        ## generate the branch name
        branchName <- paste(branch, collapse="∩")
        
        if (length(mutBranch[,1]) == 0){
            ## fix the problem when no intersection found
            next()
        } else{
            ## label the intersection(set) of the branch
            mutSigRef[which(
                mutSigRef[,1] %in% mutBranch[,1]), 2] <- branchName
            
            ## duplicate the same mutation
            mutBranchIntersection <- mutSigRef[which(
                mutSigRef$Sample == branchName & 
                    (!duplicated(mutSigRef$mut_id))), ]
            mutBranches <- rbind(mutBranches, mutBranchIntersection)
            mutBranchesOutput[[branchCounter]] <- subset(
                mutBranchIntersection,select=-c(Sample))
            listBranchName <- c(branchName, listBranchName)
            ## get the mutational signature of the branch
            mutSigsOutput <- .branchMutationalSig(mutBranches, mutSigsOutput, 
                                                      branch, branchName, 
                                                      patientID, driver_genes, 
                                                      driverGenesDir, mutThreshold, 
                                                      refBuild)
        }
        
    }
    names(mutBranchesOutput) <- listBranchName
    ## return the data frame of mutational signature for all branches
    return(list(mutSigsOutput, mutBranchesOutput))
}


## Weight mutational Signature of each branch
.branchMutationalSig <- function(mutSigRef, mutSigsOutput, 
                                   branch, branchName, 
                                   patientID, driverGenes, 
                                   driverGenesDir, mutThreshold, 
                                   refBuild){
    if (length(mutSigRef[which(
        mutSigRef$Sample == branchName), 1]) < mutThreshold){
        sigsMaxName <- "No.Signature"
        sigsMaxProb <- 0
    }else{
        ## deconstructSigs
        sigsInput <- suppressWarnings(mut.to.sigs.input(mut.ref=mutSigRef, 
                                                         sample.id="Sample", 
                                                         chr="chr", 
                                                         pos="pos", 
                                                         ref="ref", 
                                                         alt="alt",
                                                         bsg=get(refBuild)))
        sigsWhich <- whichSignatures(tumor.ref=sigsInput, 
                                      signatures.ref=signatures.cosmic, 
                                      sample.id=branchName,
                                      contexts.needed=TRUE)
        ## get mutational signature with max weight
        sigsMax <- sigsWhich[["weights"]][which.max(sigsWhich[["weights"]])]
        sigsMaxName <- colnames(sigsMax)
        sigsMaxProb <- sigsMax[,1]
    }
    
    ## vectorize branch name
    # branch <- gsub(paste(patientID,"-",sep=""), "", branch)
    
    ## figure out putative driver genes
    if (typeof(driverGenesDir) == "character"){
        ## read putative driver genes' list
        driverGenes <- as.character(read.table(
            driverGenesDir, quote="")[,1])
        ## filter potative driver genes of each branch
        pdgMut <- mutSigRef[which(
            mutSigRef$Sample == branchName &
                as.character(mutSigRef$Hugo_Symbol) %in% driverGenes),]
        pdgBranch <- as.character(pdg.mut$Hugo_Symbol)
        ## collect branches' mutataional signature and potative driver genes
        mutSigsBranch <- data.frame(
            branch=I(list(branch)), 
            sig=sigsMaxName, 
            mut.num=length(
                mutSigRef[which(
                    mutSigRef$Sample == branchName), 1]), 
            sig.prob=sigsMaxProb, 
            putative_driver_genes=I(list(pdg.branch)))
    } else{
        mutSigsBranch <- data.frame(
            branch=I(list(branch)), 
            sig=sigsMaxName, 
            mut.num=length(mutSigRef[which(
                mutSigRef$Sample == branchName), 1]), 
            sig.prob=sigsMaxProb)
    }
    ## collect branches' mutataional signature information
    rbind(mutSigsOutput, mutSigsBranch)
}
