#' get the sorted binary| ccf mutation matrix
## read mutation matraix, trans to binary matrix. 0 represents mutation absent???
## add a new column(0) representing normal sample
getMutMatrix <- function(mafData, use.ccf = FALSE){
    M <- mafData %>% 
    tidyr::unite("mutation_id", c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), sep = ":", remove = FALSE) %>% 
    dplyr::select(mutation_id, Tumor_Sample_Barcode)

    M$mutation <- 1
    mutBinary <- suppressMessages(tidyr::spread(M, Tumor_Sample_Barcode, mutation)) %>%
      dplyr::select(-mutation_id) 
    mutBinary[!is.na(mutBinary)] <- 1
    mutBinary[is.na(mutBinary)] <- 0
    mutBinary$NORMAL <- 0
    mutBinary <- t(apply(mutBinary, 2, as.numeric))

    sampleOrder <- sort(rowSums(mutBinary), decreasing=TRUE, index.return=TRUE)$ix
    scoreCol <- function(x) {
      score <- 0;
      for(i in 1:length(x)) {
        if(x[i]) {
          score <- score + 2^(length(x)-i)
        }
      }
      return(score)
    }
    scores <- apply(mutBinary[sampleOrder, ], 2, scoreCol)    
    geneOrder <- sort(scores, decreasing = TRUE, index.return = TRUE)$ix
    
    mutMatrix <- t(mutBinary[sampleOrder, geneOrder])
  
  ## get CCF matrix
  if(use.ccf & !"CCF" %in% colnames(mafData)){
    stop("No CCF data was found when during Maf object generation")
  }

  if(use.ccf & "CCF" %in% colnames(mafData)){
    M$CCF <- mafData$CCF
    M[is.na(M$CCF),'CCF'] <- 2

    mutCCF <- select(M, -mutation) %>%
      tidyr::spread(key=Tumor_Sample_Barcode, value=CCF) %>%
      dplyr::select(-mutation_id)

    mutCCF[is.na(mutCCF)] <- 0
    mutCCF$NORMAL <- 0
    mutCCF <- apply(mutCCF, 2, as.numeric)
    mutMatrix <- mutCCF[geneOrder, sampleOrder]
  }

  return(mutMatrix)
}
