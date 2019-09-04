#' @import ggplot2
#' @export mutStackPlot





oncogeneListFile <- "/home/ninomoriaty/R_Project/MesKit/inst/oncogene.list.txt"
tsgListFile <- "/home/ninomoriaty/R_Project/MesKit/inst/TSG.list.txt"

mutStackPlot <- function(maf, oncogeneListFile, tsgListFile) {
    ## prepare maf data with mut.id
    patientID <- maf@patientID
    mafData <- maf@data
    mutId <- dplyr::select(tidyr::unite(mafData, "mut.id", 
                                        Hugo_Symbol, Chromosome, 
                                        Start_Position, End_Position, 
                                        Reference_Allele, Tumor_Seq_Allele2, 
                                        sep=":"), mut.id)
    mafData <- cbind(mafData, mutId)
    tsbLs <- as.character(unique(mafData$Tumor_Sample_Barcode))
    
    ## read oncogeneList and tsgList
    oncogeneLs <- as.character(read.table(oncogeneListFile, header = TRUE, quote="", sep="\t")[,1])
    tsgLs <- as.character(read.table(tsgListFile, header = TRUE, quote="", sep="\t")[,1])
    
    ## generate mutID-labeled mafData 
    oncogeneData <- .stackDataFilter(mafData, oncogeneLs, geneType="Oncogene")
    tsgData <- .stackDataFilter(mafData, tsgLs, geneType="TSG")
    plotData <- rbind(oncogeneData, tsgData)
    
    ## generate stack plot
    ggplot(plotData, aes(x=GeneType, y=VariantNum,fill=VariantCategory)) + 
        geom_bar(stat="identity",position="stack", width = 0.55) + 
        ggtitle(paste(patientID, "'s Stack Plot", sep=""))+ 
        theme(axis.ticks.length=unit(0.5,'cm')) + 
        guides(fill=guide_legend(title=NULL)) # silence the legend title
    
}

.stackDataFilter <- function(mafData, geneLs, geneType) {
    ## data preparation
    mafDataGene <- mafData[which(as.character(mafData$Hugo_Symbol) %in% geneLs),]
    mafDataGene$tsbID <- rep("", nrow(mafDataGene))
    mutIdLs <- as.character(unique(mafDataGene$mut.id))
    
    ## generate mutid
    for (mutid in mutIdLs){
        mafDataGene[which(mafDataGene$mut.id == mutid), ]$tsbID <- 
            paste(unique(mafDataGene[which(mafDataGene$mut.id == mutid), ]$Tumor_Sample_Barcode), collapse=",")
    }
    
    ## stack: private
    GenePrivate <- mafDataGene[which(lapply(strsplit(mafDataGene$tsbID, split=','), length) == 1), ]
    ## stack: shared
    GeneShared <- mafDataGene[which(lapply(strsplit(mafDataGene$tsbID, split=','), length) == length(tsbLs)), ]
    ## stack: partial-shared
    GeneParitalShared <- mafDataGene[which(lapply(strsplit(mafDataGene$tsbID, split=','), length) > 1 & 
                                               lapply(strsplit(mafDataGene$tsbID, split=','), length) < length(tsbLs)), ]
    
    ## generate plot data for ggplot
    VariantNum <- c(nrow(GenePrivate), nrow(GeneShared), nrow(GeneParitalShared))
    VariantCategory <- c("Private", "Shared", "ParitalShared")
    stackData <- data.frame(VariantCategory, VariantNum)
    stackData$GeneType <- rep(geneType, nrow(stackData))
    return(stackData)
}
