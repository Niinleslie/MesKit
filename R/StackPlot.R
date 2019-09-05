#' @import ggplot2 ggsci
#' @export mutStackPlot



maf.File <- system.file("extdata/multi_lesion/maf", "311252.maf", package = "Meskit")
sampleInfo.File <- system.file("extdata/multi_lesion", "sample_info.txt", package = "Meskit")
maf <- readMaf(mafFile=maf.File, sampleInfoFile=sampleInfo.File, refBuild="hg19")
oncogeneListFile <- system.file("extdata/", "oncogene.list.txt", package = "Meskit")
tsgListFile <- system.file("extdata/", "TSG.list.txt", package = "Meskit")
mutStackPlot(maf, oncogeneListFile, tsgListFile, themeOption="aaas")

mutStackPlot <- function(maf, oncogeneListFile, tsgListFile, themeOption) {
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
    oncogeneData <- .stackDataFilter(mafData, oncogeneLs, tsbLs, geneType="Oncogene")
    tsgData <- .stackDataFilter(mafData, tsgLs, tsbLs, geneType="TSG")
    plotData <- rbind(oncogeneData, tsgData)
    
    

    ## generate stack plot
    ggsciFillPalette <- eval(parse(text="scale_fill_aaas()"))
    ggplot(plotData, aes(x=GeneType, y=VariantNum,fill=VariantCategory)) + 
        theme_bw() +
        geom_bar(stat="identity",position="stack", width = 0.55) + 
        theme(plot.title=element_text(size=15, face='bold'), 
              title=element_text(size=16), 
              text=element_text(size=18),  
              panel.grid=element_blank(), 
              panel.border=element_blank(), 
              axis.line=element_line(size=0.25), 
              axis.title.x = element_blank(), 
              axis.ticks.length=unit(0.5,'cm'), 
              legend.title=element_blank()) + 
        ggtitle(paste("Variants of oncogenes and TSGs in patient ", patientID, sep="")) + 
        labs(y = "Variant Number") + 
        ggsciFillPalette
}

.stackDataFilter <- function(mafData, geneLs, tsbLs, geneType) {
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
    VariantCategory <- c("Private", "Shared", "Parital-Shared")
    stackData <- data.frame(VariantCategory, VariantNum)
    stackData$GeneType <- rep(geneType, nrow(stackData))
    stackData$VariantCategory <- factor(stackData$VariantCategory, levels = c("Private", "Parital-Shared", "Shared"))
    return(stackData)
}
