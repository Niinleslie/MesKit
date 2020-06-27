## read mutation matraix, trans to binary matrix. 0 represents mutation absent???
## add a new column(0) representing normal sample
getMutMatrix <- function(maf_data, use.ccf = FALSE){
    # maf_data <- maf$`38`@data
    M <- maf_data %>% 
    tidyr::unite("mutation_id", c("Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"), sep = ":", remove = FALSE) %>% 
    dplyr::select(mutation_id, Tumor_Sample_Barcode)
    
    # print(M[450:451,])
    M$mutation <- 1
    mutBinary <- suppressMessages(tidyr::spread(M, Tumor_Sample_Barcode, mutation))
    mut.id <- mutBinary$mutation_id
    mutBinary <- dplyr::select(mutBinary, -mutation_id)
  

    mutBinary[!is.na(mutBinary)] <- 1
    mutBinary[is.na(mutBinary)] <- 0
    if(!"NORMAL" %in% colnames(mutBinary)){
        mutBinary$NORMAL <- 0
    }
    # mutBinary$NORMAL <- 0

    mutBinary <- apply(mutBinary, 2, as.numeric)
    
    if(is(mutBinary, "numeric")){
        mutBinary <- t(as.matrix(mutBinary))
    }
    
    mutBinary.trans <- t(mutBinary)

    sampleOrder <- sort(rowSums(mutBinary.trans), decreasing=TRUE, index.return=TRUE)$ix
    scoreCol <- function(x) {
      score <- 0;
      for(i in seq_len(length(x))) {
        if(x[i]) {
          score <- score + 2^(length(x)-i)
        }
      }
      return(score)
    }
    
    if(ncol(mutBinary.trans) == 1){
        sampleOrder <- sort(mutBinary,decreasing = TRUE,index.return=TRUE)$ix
        mutBinary <- t(mutBinary[,sampleOrder])
        rownames(mutBinary) <- mut.id
        mutMatrix <- mutBinary
    }else{
        scores <- apply(mutBinary.trans[sampleOrder, ], 2, scoreCol) 
        geneOrder <- sort(scores, decreasing = TRUE, index.return = TRUE)$ix
        rownames(mutBinary) <- mut.id
        mutMatrix <- mutBinary[geneOrder, sampleOrder]
    }

  
  ## get CCF matrix
  if(use.ccf & !"CCF" %in% colnames(maf_data)){
    stop("Error: no CCF data was found from Maf object")
  }

  if(use.ccf & "CCF" %in% colnames(maf_data)){
    M$CCF <- maf_data$CCF
    M[is.na(M$CCF),'CCF'] <- 2

    mutCCF <- dplyr::select(M, -mutation) %>%
      tidyr::spread(key=Tumor_Sample_Barcode, value=CCF)

    mut.id <- mutCCF$mutation_id
    mutCCF <- dplyr::select(mutCCF, -mutation_id)

    mutCCF[is.na(mutCCF)] <- 0
    mutCCF$NORMAL <- 0
    mutCCF <- apply(mutCCF, 2, as.numeric)
    rownames(mutCCF) <- mut.id
    mutMatrix <- mutCCF[geneOrder, sampleOrder]
  }

  return(mutMatrix)
}
