#' Check Mutational Signature for each branch of phylogenetic tree
#' @description Read maf file as data.frame. Define branches' set relationship by re-labeling their tumor sample 
#' barcode from the smallest set. Calcualte each branch's mutational signature weight according to cosmic reference
#' and pick the maxium. Return a data frame of each set/branch's mutational signature.
#' 
#' @param maf_file specify a maf document/directory as the input of the function
#' @param branch_file specify a txt document/directory as the input of the branches (needed to be refined)
#' @return data frame of each set/branch's mutational signature.
#'
#' @examples
#' \dontrun{
#' Mutational_sigs_tree(maf_file, branch_file)
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
Mutational_sigs_tree <- function(maf_file, branch_file, driver_genes_dir = FALSE, mut.threshold = 50){
  # read .maf file
  maf_input <- read.table(maf_file, quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # get mutationalSigs-related  infomation
  dat.sample <- data.frame(as.character(maf_input[,ncol(maf_input)]), stringsAsFactors=FALSE)
  dat.chr <- data.frame(as.character(maf_input[,2]), stringsAsFactors=FALSE)
  dat.chr[,1] <- paste("chr", dat.chr[,1], sep="")
  dat.pos.start <- maf_input[,3]
  dat.pos.end <- maf_input[,4]
  dat.ref <- maf_input[,7]
  dat.alt <- maf_input[,9]
  dat.num <- 1:length(dat.alt)
  dat.mutgene <-  maf_input$Hugo_Symbol
  mut.sig.ref <- data.frame(dat.num, dat.sample, dat.chr, dat.pos.start, dat.pos.end, dat.ref, dat.alt, dat.mutgene)
  colnames(mut.sig.ref) <- c("ID", "Sample", "chr", "pos", "pos_end", "ref", "alt", "Hugo_Symbol")
  
  # Synchronize sample name with mut.sig.ref
  patientID = strsplit(as.character(maf_input$Tumor_Sample_Barcode[1]), "-")[[1]][1]
  ID_prefix = paste(" ", patientID, "-", sep = "")
  
  # get branch infomation
  branch_input <- gsub("∩", ID_prefix, readLines(branch_file))
  branches <- strsplit(as.character(paste(patientID, "-", branch_input, sep = "")), split=" ")
  
  # output collection
  mut.sigs.output <- data.frame()
  
  # generate mutational signautres for different branches
  for (branch_counter in length(branches):1){
    # generate a single branch
    branch <- Filter(Negate(is.na), branches[[branch_counter]])
    mut.branch <- mut.sig.ref[which(mut.sig.ref$Sample %in% branch), ]
    
    for (tsb in branch){
      # generate the intersection(set) of the branch
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
      # get the mutational signature of the branch
      ### However, this part could be optimized as sigs.input should be just calculated once. ###
      mut.sigs.output <- Mutational_sigs_branch(mut.sig.ref, mut.sigs.output, branch, branch_name, patientID, driver_genes, driver_genes_dir, mut.threshold)
    }
    
  }
  # return the data frame of mutational signature for all branches
  mut.sigs.output
}

  
# Weight mutational Signature of each branch
Mutational_sigs_branch <- function(mut.sig.ref, mut.sigs.output, branch, branch_name, patientID, driver_genes, driver_genes_dir, mut.threshold){
  if (length(mut.sig.ref[which(mut.sig.ref$Sample == branch_name), 1]) < mut.threshold){
    sigs.max <- "No Signature"
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
    sigs.max <- colnames(sigs.which[["weights"]][which.max(sigs.which[["weights"]])])
  }
  
  # vectorize branch name
  branch <- gsub(paste(patientID,"-",sep=""), "", branch)
  
  # figure out putative driver genes
  if (typeof(driver_genes_dir) == "character"){
    # read putative driver genes' list
    driver_genes <- as.character(read.table(driver_genes_dir, quote = "")[,1])
    # filter potative driver genes of each branch
    pdg.mut <- mut.sig.ref[which(mut.sig.ref$Sample == branch_name &
                                   as.character(mut.sig.ref$Hugo_Symbol) %in% driver_genes),]
    pdg.branch <- as.character(pdg.mut$Hugo_Symbol)
    # collect branches' mutataional signature and potative driver genes information
    mut.sigs.branch <- data.frame(branch = I(list(branch)), mut.sig = sigs.max, putative_driver_genes = I(list(pdg.branch)))
  } else{
    mut.sigs.branch <- data.frame(branch = I(list(branch)), mut.sig = sigs.max)
  }
  # collect branches' mutataional signature information
  rbind(mut.sigs.output, mut.sigs.branch)
}


###### putative driver genes ######
driver_genes_dir <- "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/putative_driver_genes.txt"

###### output test ######
# setwd("/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/results")
# maf_file_ls = c(maf_file1, maf_file2, maf_file3, maf_file4, maf_file5, maf_file6, maf_file7)
# branch_file_ls = c(branch_file1, branch_file2, branch_file3, branch_file4, branch_file5, branch_file6, branch_file7)
# for (counter in 1:length(branch_file_ls)){
#   result <- Mutational_sigs_tree(maf_file_ls[counter], branch_file_ls[counter])
#   write.table(result, file=paste("output", counter,".txt", sep = ""))
# }

# sigs <- read.table("/home/ninomoriaty/Nutstore Files/Nutstore/VAF_plot_beta/Mutation_sign/output2.txt", stringsAsFactors=F, quote = "", header = TRUE, fill = TRUE, sep = ' ')
# 
# maf_file1 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/311252.snv_indel.imputed.maf"
# branch_file1 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/311252.NJtree.edges"
# 
# maf_file2 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/313544.snv_indel.imputed.maf"
# branch_file2 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/313544.NJtree.edges"
# 
# maf_file3 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/313935.snv_indel.imputed.maf"
# branch_file3 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/313935.NJtree.edges"
# 
maf_file4 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/313953.snv_indel.imputed.maf"
branch_file4 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/313953.NJtree.edges"
# 
# maf_file5 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/314007.snv_indel.imputed.maf"
# branch_file5 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/314007.NJtree.edges"
# 
# maf_file6 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/314069.snv_indel.imputed.maf"
# branch_file6 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/314069.NJtree.edges"
# 
# maf_file7 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/314155.snv_indel.imputed.maf"
# branch_file7 = "/home/ninomoriaty/Nutstore Files/Nutstore/edges_mafs/314155.NJtree.edges"

# Confirm sets of mutation
# Mutation_sets <- function(mut.sig.ref, branches){
#   # generate branch name
#   for (branch_counter in 1:length(branches)){
#     branch <- Filter(Negate(is.na), branches[,branch_counter])
#     mut.branch <- mut.sig.ref[which(mut.sig.ref$Sample %in% branch), ]
#     for (tsb in branch){
#       mut.tsb <- mut.sig.ref[which(mut.sig.ref$Sample %in% tsb), ]
#       mut.branch <- match_df(mut.branch, mut.tsb, on = c("chr", "pos", "pos_end", "ref", "alt"))
#     }
#     branch_name <- paste(branch, collapse = "+")
#     Mutational_sigs_branch(mut.branch, branch_name)
#     mut.sig.ref[which(mut.sig.ref[,1] %in% mut.branch[,1]), 2] <- branch_name
#   }
#   mut.sig.ref
# }

