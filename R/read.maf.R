#' Output an new maf document with information of sample_info 
#' @description Add sample_info to the original maf file to get a new maf.The new maf file adds three pieces of information:lesion,patient and time
#'
#' @import maftools
#'
#' @param patientID patient/sample name
#' @param maf.dir specify a maf document/directory as the input of the function
#' @param sample_info_file specify a txt document/directory as the input of the sample_info (Created by yourself)
#' @return the new maf file includes information of sample_info and mut.id
#'
#' @examples
#' read.maf(patientID, maf.dir = "./data/maf", SampleInfo.dir = "./data/") 

library(maftools)

MAF2 <- setClass(Class = "MAF2", contains="MAF", slots =  c(ccf.cluster.tsv = 'data.table', ccf.loci.tsv = 'data.table'))

# directories
maf.dir = "/home/ninomoriaty/R_Project/MesKit/inst/extdata/multi_lesion/maf"
ccf.dir = "/home/ninomoriaty/R_Project/MesKit/inst/extdata/multi_lesion/ccf"
patientID = "311252"
SampleInfo.dir = "/home/ninomoriaty/R_Project/MesKit/inst/extdata/multi_lesion/sample_info.txt"

# read.maf main function
read.maf2 <- function(patientID, maf.dir, ccf.dir, SampleInfo.dir){
  # read maf file
  maf_input <- read.table(paste(maf.dir,'/',patientID,'.maf',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # read info file
  sample_info_input <-  read.table(SampleInfo.dir, quote = "", header = TRUE, fill = TRUE, sep = '')
  # read ccf file
  ccf.cluster.tsv_input <- read.table(paste(ccf.dir,'/',patientID,'.cluster.tsv',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
  ccf.loci.tsv_input <- read.table(paste(ccf.dir,'/',patientID,'.loci.tsv',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
  # Generate patient,lesion and time information in the last three columns of the maf file
  maf_input$patient = ""
  maf_input$lesion = ""
  maf_input$time = ""
  
  # combine sample_info_input with maf_input
  tsb_ls_total <- unique(sample_info_input$sample)
  tsb_ls <- tsb_ls_total[which(tsb_ls_total %in% maf_input$Tumor_Sample_Barcode)]
  for (tsb in tsb_ls) {
    patient <- sample_info_input[which(sample_info_input$sample == tsb),]$patient
    lesion <- as.character(sample_info_input[which(sample_info_input$sample == tsb),]$lesion)
    time <- as.character(sample_info_input[which(sample_info_input$sample == tsb),]$time)
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$patient <- patient
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$lesion <- lesion
    maf_input[which(maf_input$Tumor_Sample_Barcode == tsb),]$time <- time
  }
  
  maf_input$Hugo_Symbol <- as.character(maf_input$Hugo_Symbol)
  maf.data <- data.table::setDT(maf_input)
  ccf.cluster.tsv <- data.table::setDT(ccf.loci.tsv_input)
  ccf.loci.tsv <- data.table::setDT(ccf.loci.tsv_input)
  
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
  
  maf.summary <- maftools:::summarizeMaf(maf = maf.data, chatty = TRUE)
  
  
  
  m2 <- MAF2(data = maf.data, variants.per.sample = maf.summary$variants.per.sample, variant.type.summary = maf.summary$variant.type.summary,
             variant.classification.summary = maf.summary$variant.classification.summary, gene.summary = maf.summary$gene.summary,
             summary = maf.summary$summary, maf.silent = maf.silent, clinical.data = maf.summary$sample.anno, ccf.cluster.tsv = ccf.cluster.tsv, ccf.loci.tsv = ccf.loci.tsv)
  
  return(m2)
}

plotmafSummary(maf=m2, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE)


# # Error in setattr(x, "row.names", rn) : row names must be 'character' or 'integer', not 'integer'
# gs.dat = getGeneSummary(m2)
# data.table::setDF(gs.dat)
# rownames(gs.dat) = as.character(gs.dat$Hugo_Symbol)
# 
# 
# # maftools comparasion
# maftools.result <- read.maf(paste(maf.dir,'/',patientID,'.maf',sep = ""))
# plotmafSummary(maf=maftools.result, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE)







