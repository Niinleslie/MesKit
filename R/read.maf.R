<<<<<<< HEAD
#' Output an new maf document with information of sample_info 
#' @description Add sample_info to the original maf file to get a new maf.The new maf file adds three pieces of information:lesion,patient and time
#'
#' @import maftools ggplot2
#'
#' @param patientID patient/sample name
#' @param dat.dir specify a data directory as the input of the function
#' @param use.ccf = Default FALSE. TRUE if ccf data is provided
#' @param plot.mafSummary Default TRUE. FALSE if the summary figure is unnecessary
#' @return a Maf object/class includes information of sample_info and mut.id and summary figure of it
#'
#' @examples
#' read.maf(patientID, dat.dir = "./data", use.ccf = FALSE, plot.mafSummary = TRUE) 

library(maftools)
library(ggplot2)

# Maf class
Maf <- setClass(Class = "Maf", contains = "MAF", slots =  c(ccf.cluster = 'data.table', ccf.loci = 'data.table', patientID='character'))

# read.maf main function
read.Maf<- function(patientID, dat.dir = "./data", use.ccf = FALSE, plot.mafSummary = TRUE){
  # read maf file
  maf_input <- read.table(paste(dat.dir,'/maf/',patientID,'.maf',sep = ""), quot = "", header = TRUE, fill = TRUE, sep = '\t')
  # read info file
  sample_info_input <-  read.table(paste(dat.dir, '/sample_info.txt',sep = ""), quot = "", header = TRUE, fill = TRUE, sep = '')
  # read ccf file
  if (use.ccf) {
    ccf.cluster.tsv_input <- read.table(paste(dat.dir, '/ccf/', patientID, '.cluster.tsv', sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
    ccf.loci.tsv_input <- read.table(paste(dat.dir, '/ccf/', patientID, '.loci.tsv', sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
  } else {
    ccf.cluster.tsv_input <- NULL
    ccf.loci.tsv_input <- NULL
  }
  # Generate patient,lesion and time information in the last three columns of the maf file
  maf_input$patient = ""
  maf_input$lesion = ""
  maf_input$time = ""
  
  # combine sample_info_input with maf_input
  tsb_SampleInfo <- unique(sample_info_input$sample)
  tsb_ls <- tsb_SampleInfo[which(tsb_SampleInfo %in% maf_input$Tumor_Sample_Barcode)]
  for (tsb in tsb_ls) {
    patient <- sample_info_input[which(sample_info_input$sample == tsb),]$patient
    lesion <- as.character(sample_info_input[which(sample_info_input$sample == tsb),]$lesion)
    time <- as.character(sample_info_input[which(sample_info_input$sample == tsb),]$time)
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$patient <- patient
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$lesion <- lesion
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$time <- time
  }
  
  
  # fix: Error in setattr(x, "row.names", rn) : row names must be 'character' or 'integer', not 'integer'
  maf_input$Hugo_Symbol <- as.character(maf_input$Hugo_Symbol)
  
  # transform data.frame to data.table
  maf.data <- data.table::setDT(maf_input)
  ccf.cluster.tsv <- data.table::setDT(ccf.loci.tsv_input)
  ccf.loci.tsv <- data.table::setDT(ccf.loci.tsv_input)
  
  # generate maf.silent and filter maf.data
  vc.nonSilent =  c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                    "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                    "In_Frame_Ins", "Missense_Mutation")
  
  maf.silent = maf.data[!Variant_Classification %in% vc.nonSilent] #Silent variants
  if(nrow(maf.silent) > 0){
    maf.silent.vc = maf.silent[,.N, .(Tumor_Sample_Barcode, Variant_Classification)]
    maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N') #why dcast is not returning it as data.table ?
    summary.silent = data.table::data.table(ID = c('Samples',colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]),
                                            N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,2:ncol(maf.silent.vc.cast), with = FALSE])))
    
    #maf.data = maf.data[Variant_Classification %in% vc.nonSilent] #Choose only non-silent variants from main table
  }
  
  # summarize sample_info and mut.id with summarizeMaf
  maf.summary <- maftools:::summarizeMaf(maf = maf.data, chatty = TRUE)
  maf <- Maf(data = maf.data, variants.per.sample = maf.summary$variants.per.sample, variant.type.summary = maf.summary$variant.type.summary,
              variant.classification.summary = maf.summary$variant.classification.summary, gene.summary = maf.summary$gene.summary,
              summary = maf.summary$summary, maf.silent = maf.silent, clinical.data = maf.summary$sample.anno, ccf.cluster = ccf.cluster.tsv, ccf.loci = ccf.loci.tsv,patientID = patientID)
  
  # print the summary plot
  if (plot.mafSummary) {
    ggsave(plotmafSummary(maf=maf, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE), filename = paste(patientID, ".png", sep=""), width = 12, height = 9, dpi = 300)
  }
  
  return(maf)
}




# # Error in setattr(x, "row.names", rn) : row names must be 'character' or 'integer', not 'integer'
# gs.dat = getGeneSummary(m2)
# data.table::setDF(gs.dat)
# rownames(gs.dat) = as.character(gs.dat$Hugo_Symbol)
# 
# 
# # maftools comparasion
# maftools.result <- read.maf(paste(maf.dir,'/',patientID,'.maf',sep = ""))
# plotmafSummary(maf=maftools.result, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE)







=======
#' Output an new maf document with information of sample_info 
#' @description Add sample_info to the original maf file to get a new maf.The new maf file adds three pieces of information:lesion,patient and time
#'
#' @import maftools ggplot2
#'
#' @param patientID patient/sample name
#' @param dat.dir specify a data directory as the input of the function
#' @param use.ccf = Default FALSE. TRUE if ccf data is provided
#' @param plot.mafSummary Default TRUE. FALSE if the summary figure is unnecessary
#' @return a Maf object/class includes information of sample_info and mut.id and summary figure of it
#'
#' @examples
#' read.maf(patientID, dat.dir = "./data", use.ccf = FALSE, plot.mafSummary = TRUE) 

library(maftools)
library(ggplot2)

# Maf class
Maf <- setClass(Class = "Maf", contains = "MAF", slots =  c(ccf.cluster.tsv = 'data.table', ccf.loci.tsv = 'data.table'))

# directories
setwd("/home/ninomoriaty/R_Project")
patientID = "311252"

# read.maf main function
read.Maf<- function(patientID, dat.dir = "./data", use.ccf = FALSE, plot.mafSummary = TRUE){
  # read maf file
  maf_input <- read.table(paste(dat.dir,'/maf/',patientID,'.maf',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # read info file
  sample_info_input <-  read.table(paste(dat.dir, '/sample_info.txt',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '')
  # read ccf file
  if (use.ccf) {
    ccf.cluster.tsv_input <- read.table(paste(dat.dir,'/ccf/',patientID,'.cluster.tsv',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
    ccf.loci.tsv_input <- read.table(paste(dat.dir,'/ccf/',patientID,'.loci.tsv',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
  } else {
    ccf.cluster.tsv_input <- NULL
    ccf.loci.tsv_input <- NULL
  }
  # Generate patient,lesion and time information in the last three columns of the maf file
  maf_input$patient = ""
  maf_input$lesion = ""
  maf_input$time = ""
  
  # combine sample_info_input with maf_input
  tsb_SampleInfo <- unique(sample_info_input$sample)
  tsb_ls <- tsb_SampleInfo[which(tsb_SampleInfo %in% maf_input$Tumor_Sample_Barcode)]
  for (tsb in tsb_ls) {
    patient <- sample_info_input[which(sample_info_input$sample == tsb),]$patient
    lesion <- as.character(sample_info_input[which(sample_info_input$sample == tsb),]$lesion)
    time <- as.character(sample_info_input[which(sample_info_input$sample == tsb),]$time)
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$patient <- patient
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$lesion <- lesion
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$time <- time
  }
  
  
  # fix: Error in setattr(x, "row.names", rn) : row names must be 'character' or 'integer', not 'integer'
  maf_input$Hugo_Symbol <- as.character(maf_input$Hugo_Symbol)
  
  # transform data.frame to data.table
  maf.data <- data.table::setDT(maf_input)
  ccf.cluster.tsv <- data.table::setDT(ccf.loci.tsv_input)
  ccf.loci.tsv <- data.table::setDT(ccf.loci.tsv_input)
  
  # generate maf.silent and filter maf.data
  vc.nonSilent =  c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                    "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                    "In_Frame_Ins", "Missense_Mutation")
  
  maf.silent = maf.data[!Variant_Classification %in% vc.nonSilent] #Silent variants
  if(nrow(maf.silent) > 0){
    maf.silent.vc = maf.silent[,.N, .(Tumor_Sample_Barcode, Variant_Classification)]
    maf.silent.vc.cast = data.table::dcast(data = maf.silent.vc, formula = Tumor_Sample_Barcode ~ Variant_Classification, fill = 0, value.var = 'N') #why dcast is not returning it as data.table ?
    summary.silent = data.table::data.table(ID = c('Samples',colnames(maf.silent.vc.cast)[2:ncol(maf.silent.vc.cast)]),
                                            N = c(nrow(maf.silent.vc.cast), colSums(maf.silent.vc.cast[,2:ncol(maf.silent.vc.cast), with = FALSE])))
    
    maf.data = maf.data[Variant_Classification %in% vc.nonSilent] #Choose only non-silent variants from main table
  }
  
  # summarize sample_info and mut.id with summarizeMaf
  maf.summary <- maftools:::summarizeMaf(maf = maf.data, chatty = TRUE)
  maf <- Maf(data = maf.data, variants.per.sample = maf.summary$variants.per.sample, variant.type.summary = maf.summary$variant.type.summary,
              variant.classification.summary = maf.summary$variant.classification.summary, gene.summary = maf.summary$gene.summary,
              summary = maf.summary$summary, maf.silent = maf.silent, clinical.data = maf.summary$sample.anno, ccf.cluster.tsv = ccf.cluster.tsv, ccf.loci.tsv = ccf.loci.tsv)
  
  # print the summary plot
  if (plot.mafSummary) {
    ggsave(plotmafSummary(maf=maf, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE), filename = paste(patientID, ".png", sep=""), width = 12, height = 9, dpi = 300)
  }
  
  return(maf)
}


# # Error in setattr(x, "row.names", rn) : row names must be 'character' or 'integer', not 'integer'
# gs.dat = getGeneSummary(m2)
# data.table::setDF(gs.dat)
# rownames(gs.dat) = as.character(gs.dat$Hugo_Symbol)
# 
# 
# # maftools comparasion
# maftools.result <- read.maf(paste(maf.dir,'/',patientID,'.maf',sep = ""))
# plotmafSummary(maf=maftools.result, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE)







>>>>>>> 4efd87d3fb7554c11df72b0bed97d2640207d2e3
