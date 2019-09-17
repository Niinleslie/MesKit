#' get the sorted binary| ccf mutation matrix
maf_preprocess <- function(maf.dat, use.indel = F, use.ccf = F, ccf, ccf.mutation.id, ccf.mutation.sep){
  if(!use.indel){
    maf.dat <- maf.dat[which(maf.dat$Variant_Type == "SNP"),]
  }
  mut.id <- tidyr::unite(maf.dat, "mut.id", Hugo_Symbol, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":")$mut.id
  M <- data.frame(mut.id = mut.id, sample = maf.dat$Tumor_Sample_Barcode)

  if(use.ccf){
    if(is.null(ccf)){
      stop("Missing ccf file. Check whether maf@ccf.loci is NULL")
    }
    ccf.mut.id <- as.vector(select(tidyr::unite(maf.dat, "ccf.mut.id", ccf.mutation.id, sep = ccf.mutation.sep), ccf.mut.id))    
    M <- cbind(M, ccf.mut.id=ccf.mut.id)
    M <- select(merge(x = M, y = ccf, by.x = c("ccf.mut.id","sample"), by.y = c("mutation_id", "sample_id"), all.x =T), mut.id, sample, cellular_prevalence)
    colnames(M) <- c("mut.id", "sample", "CCF")
  }

  return(M)
}


##read mutation matraix, trans to 0-1 binary matrix???0 represents mutation absent???
##add a new column(0) representing normal sample
##ccf return from NJtree
mut_ccf_sort <- function(maf.dat, ccf, use.indel, ccf.mutation.id , ccf.mutation.sep){
  index_row_col <- mut_binary_sort(maf.dat, use.indel = use.indel, returnOrder = TRUE)
  M <- maf_preprocess(maf.dat, use.indel = use.indel, use.ccf = TRUE, ccf, ccf.mutation.id = ccf.mutation.id, ccf.mutation.sep = ccf.mutation.sep)
  M[is.na(M$CCF),'CCF'] <- 2

  mut_samples <- suppressMessages(dcast(M, mut.id~sample, value.var = "CCF")[,-1])
  #mut_samples[mut_samples == 0 ] <- 2
  mut_samples[is.na(mut_samples)] <- 0
  mut_samples$NORMAL <- 0

  mut_samples<- apply(mut_samples, 2, as.numeric)
  mut_sort <- mut_samples[index_row_col[[1]], index_row_col[[2]]]
  return(mut_sort)
}
##sort the matrix; return row order and coloum order
mut_binary_sort <- function(maf.dat, use.indel = FALSE, returnOrder = FALSE){
  M <- maf_preprocess(maf.dat, use.indel = use.indel)
  mut_samples <- suppressMessages(dcast(M, mut.id~sample)[,-1])
  mut_samples[!is.na(mut_samples)] <- 1
  mut_samples[is.na(mut_samples)] <- 0
  mut_samples$NORMAL <- 0
  mut_binary<- apply(mut_samples, 2, as.numeric)

  mut_binary <- t(mut_binary)
  sampleOrder <- sort(rowSums(mut_binary), decreasing = TRUE, index.return = TRUE)$ix
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i)
      }
    }
    return(score)
  }
  scores <- apply(mut_binary[sampleOrder, ], 2, scoreCol)
  geneOrder <- sort(scores, decreasing = TRUE, index.return = TRUE)$ix
  mut_sort <- t(mut_binary[sampleOrder,geneOrder])
  if(returnOrder){
    return(list(geneOrder, sampleOrder))
  }
  else{
    return(mut_sort)
  }
}



