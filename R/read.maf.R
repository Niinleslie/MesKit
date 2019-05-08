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
  
  maf_input$Variant_Classification = as.factor(as.character(maf_input$Variant_Classification))
  maf_input$Variant_Type = as.factor(as.character(maf_input$Variant_Type))
  
  nGenes = length(unique(maf_input$Hugo_Symbol))
  maf.tsbs = levels(maf_input$Tumor_Sample_Barcode)
  nSamples = length(levels(maf_input$Tumor_Sample_Barcode))
  
  maf.data <- data.table::setDT(maf_input)
  
  # MAF:data
  tsb = maf_input[,.N, Tumor_Sample_Barcode]
  colnames(tsb)[2] = 'Variants'
  tsb = tsb[order(tsb$Variants, decreasing = TRUE),]
  
  # MAF:
  
  # mafSummary = summarizeMaf(maf = maf_input, anno = clinicalData, chatty = verbose)
  # 
  # m = MAF(data = maf, variants.per.sample = mafSummary$variants.per.sample, variant.type.summary = mafSummary$variant.type.summary,
  #         variant.classification.summary = mafSummary$variant.classification.summary, gene.summary = mafSummary$gene.summary,
  #         summary = mafSummary$summary, maf.silent = maf.silent, clinical.data = mafSummary$sample.anno)
  

  h <- function(x){
    Hugo_Symbol=x[1]
    Chromosome=x[2]
    Start_Position=x[3]
    End_Position=x[4]
    Tumor_Seq_Allele1=x[8]
    Tumor_Seq_Allele2=x[9]
    mut.id=paste(Hugo_Symbol,Chromosome,Start_Position,End_Position,Tumor_Seq_Allele1,Tumor_Seq_Allele2,sep = '_')
    x[19]=mut.id
    #Find the location of sample at sample_info
    row <- which(x[15]==sample_info_input$sample)
    x[16]=as.character(sample_info_input$patient[row])
    x[17]=as.character(sample_info_input$lesion[row])
    x[18]=as.character(sample_info_input$time[row])
    return(x)
  }
  maf_input=t(apply(maf_input,1,h))
  return(maf_input)
}










