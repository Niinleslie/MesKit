#' Output an new maf document with information of sample_info 
#' @description Add sample_info to the original maf file to get a new maf.The new maf file adds three pieces of information:lesion,patient and time
#' 
#' @param patientID patient/sample name
#' @param maf.dir specify a maf document/directory as the input of the function
#' @param sample_info_file specify a txt document/directory as the input of the sample_info (Created by yourself)
#' @return the new maf file includes information of sample_info and mut.id
#'
#' @examples
#' read.maf(patientID, maf.dir = "./data/maf", SampleInfo.dir = "./data/") 

library(maftools)

## MAF object
MAF <- setClass(Class = 'MAF', slots =  c(data = 'data.table', variants.per.sample = 'data.table', variant.type.summary = 'data.table',
                                          variant.classification.summary = 'data.table', gene.summary = 'data.table',
                                          summary = 'data.table', maf.silent = 'data.table'))

setMethod(f = 'show', signature = 'MAF', definition = function(object){
  cat(paste('An object of class ', class(object), "\n"))
  print(object@summary)
})


# directories
maf.dir = "/home/ninomoriaty/R_Project/MesKit/inst/extdata/multi_lesion/maf"
patientID = "311252"
SampleInfo.dir = "/home/ninomoriaty/R_Project/MesKit/inst/extdata/multi_lesion/sample_info.txt"


# maftools comparasion
maftools.result <- read.maf(paste(maf.dir,'/',patientID,'.maf',sep = ""))
# vcs = getSampleSummary(maftools.result)
maftools.result@summary
maftools.result@variant.classification.summary
m <- MAF(data = maftools.result@data, variants.per.sample = maftools.result@variants.per.sample, variant.type.summary = maftools.result@variant.type.summary,
         variant.classification.summary = maftools.result@variant.classification.summary, gene.summary = maftools.result@gene.summary,
         summary = maftools.result@summary, maf.silent = maftools.result@maf.silent)

plotmafSummary(maf = m, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

# read.maf main function
read.maf <- function(patientID, maf.dir, SampleInfo.dir){
  # read maf file
  maf_input <- read.table(paste(maf.dir,'/',patientID,'.maf',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t')
  # read info file
  sample_info_input <-  read.table(SampleInfo.dir, quote = "", header = TRUE, fill = TRUE, sep = '')
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
  
  # 
  
}










