#' Check Mutational Signature for each branch of phylogenetic tree
#' @description Read maf file as data.frame. Define branches' set relationship by re-labeling their tumor sample 
#' barcode from the smallest set. Calcualte each branch's mutational signature weight according to cosmic reference
#' and pick the maxium. Return a data frame of each set/branch's mutational signature.
#' 
#' @import reshape2 BSgenome BSgenome.Hsapiens.UCSC.hg19 GenomeInfoDb grDevices graphics utils deconstructSigs
#' @import plyr
#' 
#' @param maf_file specify a maf document/directory as the input of the function
#' @param branch_file specify a txt document/directory as the input of the branches (needed to be refined)
#' @return data frame of each set/branch's mutational signature.
#'
#' @examples
#' \dontrun{
#' Mutational_sigs_tree(maf_file, branch_file)
#' Mutational_sigs_tree(maf_file, branch_file, driver_genes_dir)
#' Mutational_sigs_tree(maf_file, branch_file, driver_genes_dir, mut.threshold = 30)
#'}


# import pkgs
# dependencies of deconstructSigs
library(reshape2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomeInfoDb)
library(grDevices)
library(graphics)
library(utils)
library(deconstructSigs)
# data frame needed
library(plyr)

# main function
# Usage: Mutational_Sigs_branch(maf_file, samples_vector)
Mutational_sigs_tree <- function(maf.dat, branch, driver_genes_dir = FALSE, mut.threshold = 50){
  maf_input <- maf.dat
  # get mutationalSigs-related  infomation
  dat.sample <- data.frame(as.character(maf_input$Tumor_Sample_Barcode), stringsAsFactors=FALSE)
  dat.chr <- data.frame(as.character(maf_input$Chromosome), stringsAsFactors=FALSE)
  dat.chr[,1] <- paste("chr", dat.chr[,1], sep="")
  dat.pos.start <- maf_input$Start_Position
  dat.pos.end <- maf_input$End_Position
  dat.ref <- maf_input$Reference_Allele
  dat.alt <- maf_input$Tumor_Seq_Allele2
  dat.num <- 1:length(dat.alt)
  dat.mutgene <-  maf_input$Hugo_Symbol
  mut.sig.ref <- data.frame(dat.num, dat.sample, dat.chr, dat.pos.start, dat.pos.end, dat.ref, dat.alt, dat.mutgene)
  colnames(mut.sig.ref) <- c("ID", "Sample", "chr", "pos", "pos_end", "ref", "alt", "Hugo_Symbol")
  
  # Synchronize sample name with mut.sig.ref
  patientID = strsplit(as.character(maf_input$Tumor_Sample_Barcode[1]), "-")[[1]][1]
  ID_prefix = paste(" ", patientID, "-", sep = "")
  
  # get branch infomation
  branches <- strsplit(branch,split='∩')
  
  # output collection
  mut.sigs.output <- data.frame()
  mut.branches <- data.frame()
  # generate mutational signautres for different branches
  for (branch_counter in length(branches):1){
    # generate a single branch
    branch <- Filter(Negate(is.na), branches[[branch_counter]])
    mut.branch <- mut.sig.ref[which(mut.sig.ref$Sample %in% branch), ]
    
    for (tsb in branch){
      # generate the intersection(set) of the branch(different from duplication)
      mut.tsb <- mut.sig.ref[which(mut.sig.ref$Sample %in% tsb), ]
      mut.branch <- match_df(mut.branch, mut.tsb, on = c("chr", "pos", "pos_end", "ref", "alt"))
    }
    
    # generate the branch name
    branch_name <- paste(branch, collapse = "+")
    
    if (length(mut.branch[,1]) == 0){
      # fix the problem when no intersection found
      next()
    } else{
      # label the intersection(set) of the branch
      mut.sig.ref[which(mut.sig.ref[,1] %in% mut.branch[,1]), 2] <- branch_name
      # duplicate the same mutation
      mut.branch.intersection <- mut.sig.ref[which(mut.sig.ref$Sample == branch_name 
                                                   & (!duplicated(mut.sig.ref$chr) 
                                                      | !duplicated(mut.sig.ref$pos)
                                                      | !duplicated(mut.sig.ref$pos_end)
                                                      | !duplicated(mut.sig.ref$ref)
                                                      | !duplicated(mut.sig.ref$ alt))),]
      mut.branches <- rbind(mut.branches, mut.branch.intersection)
      # get the mutational signature of the branch
      ### However, this part could be optimized as sigs.input should be just calculated once. ###
      mut.sigs.output <- Mutational_sigs_branch(mut.branches, mut.sigs.output, branch, branch_name, patientID, driver_genes, driver_genes_dir, mut.threshold)
    }
    
  }
  # return the data frame of mutational signature for all branches
  mut.sigs.output
}


# Weight mutational Signature of each branch
Mutational_sigs_branch <- function(mut.sig.ref, mut.sigs.output, branch, branch_name, patientID, driver_genes, driver_genes_dir, mut.threshold){
  if (length(mut.sig.ref[which(mut.sig.ref$Sample == branch_name), 1]) < mut.threshold){
    sigs.max.name <- "No.Signature"
    sigs.max.prob <- 0
  }else{
    # deconstructSigs
    sigs.input <- suppressWarnings(mut.to.sigs.input(mut.ref = mut.sig.ref, 
                                                     sample.id = "Sample", 
                                                     chr = "chr", 
                                                     pos = "pos", 
                                                     ref = "ref", 
                                                     alt = "alt"))
    sigs.which <- whichSignatures(tumor.ref = sigs.input, 
                                  signatures.ref = signatures.cosmic, 
                                  sample.id = branch_name,
                                  contexts.needed = TRUE)
    # get mutational signature with max weight
    sigs.max <- sigs.which[["weights"]][which.max(sigs.which[["weights"]])]
    sigs.max.name <- colnames(sigs.max)
    sigs.max.prob <- sigs.max[,1]
  }
  
  # # vectorize branch name
  # branch <- gsub(paste(patientID,"-",sep=""), "", branch)
  
  # figure out putative driver genes
  if (typeof(driver_genes_dir) == "character"){
    # read putative driver genes' list
    driver_genes <- as.character(read.table(driver_genes_dir, quote = "")[,1])
    # filter potative driver genes of each branch
    pdg.mut <- mut.sig.ref[which(mut.sig.ref$Sample == branch_name &
                                   as.character(mut.sig.ref$Hugo_Symbol) %in% driver_genes),]
    pdg.branch <- as.character(pdg.mut$Hugo_Symbol)
    # collect branches' mutataional signature and potative driver genes information
    mut.sigs.branch <- data.frame(branch = I(list(branch)), sig = sigs.max.name, mut.num = length(mut.sig.ref[which(mut.sig.ref$Sample == branch_name), 1]), sig.prob = sigs.max.prob, putative_driver_genes = I(list(pdg.branch)))
  } else{
    mut.sigs.branch <- data.frame(branch = I(list(branch)), sig = sigs.max.name, mut.num = length(mut.sig.ref[which(mut.sig.ref$Sample == branch_name), 1]), sig.prob = sigs.max.prob)
  }
  # collect branches' mutataional signature information
  rbind(mut.sigs.output, mut.sigs.branch)
}
