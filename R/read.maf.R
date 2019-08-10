#' Output an new maf document with information of sample_info 
#' @description Add sample_info to the original maf file to get a new maf.The new maf file adds three pieces of information:lesion,patient and time
#'
#' @import maftools ggplot2
#'
#' @param patientID patient/sample name
#' @param dat.dir specify a data directory as the input of the function
#' @param use.ccf = Default FALSE. TRUE if ccf data is provided
#' @param plot.mafSummary Default TRUE. FALSE if the summary figure is unnecessary
#' @param ref.build Default hg19. specify a referential genome including BSgenome.Hsapiens.UCSC.hg19 and BSgenome.Hsapiens.UCSC.hg38 
#' @return a Maf object/class includes information of sample_info and mut.id and summary figure of it
#'
#' @examples
#' read.maf(patientID, dat.dir = "./data", use.ccf = FALSE, plot.mafSummary = TRUE) 

# Maf class
Maf <- setClass(Class = "Maf", contains = "MAF", slots =  c(ccf.cluster = 'data.table', ccf.loci = 'data.table', patientID = 'character', ref.build='character'))

# read.maf main function
read.Maf <- function(patientID, maf.dir, sample_info.dir, ccf.cluster.dir = NULL, ccf.loci.dir = NULL, plot.mafSummary = TRUE, ref.build = "hg19"){
  # read maf file
  maf_input <- read.table(maf.dir, quot = "", header = TRUE, fill = TRUE, sep = '\t')
  # read info file
  sample_info_input <-  read.table(sample_info.dir, quot = "", header = TRUE, fill = TRUE, sep = '', stringsAsFactors = F)
  # read ccf file
  if (!is.null(ccf.cluster.dir)&!is.null(ccf.loci.dir)) {
    ccf.cluster.tsv_input <- read.table(ccf.cluster.dir, quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
    ccf.loci.tsv_input <- read.table(ccf.loci.dir, quote = "", header = TRUE, fill = TRUE, sep = '\t', stringsAsFactors=F)
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
    patient <- as.character(sample_info_input[which(sample_info_input$sample == tsb),]$patient) 
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
  ccf.cluster.tsv <- data.table::setDT(ccf.cluster.tsv_input)
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
  maf.summary <- suppressMessages(maftools:::summarizeMaf(maf = maf.data, chatty = TRUE))
  maf <- Maf(data = maf.data, variants.per.sample = maf.summary$variants.per.sample, variant.type.summary = maf.summary$variant.type.summary,
              variant.classification.summary = maf.summary$variant.classification.summary, gene.summary = maf.summary$gene.summary,
              summary = maf.summary$summary, maf.silent = maf.silent, clinical.data = maf.summary$sample.anno, ccf.cluster = ccf.cluster.tsv, ccf.loci = ccf.loci.tsv, patientID = patientID, ref.build = ref.build)
  
  # print the summary plot
  if (plot.mafSummary) {
    pic <- plotmafSummary(maf=maf, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE)
    #ggsave(plotmafSummary(maf=maf, rmOutlier = TRUE, addStat='median', dashboard = T, titvRaw = FALSE), filename = paste(patientID, ".VariantSummary.png", sep=""), width = 12, height = 9, dpi = 800, path = "./output")
  }
  
  return(maf)
}


