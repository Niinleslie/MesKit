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


read.maf <- function(patientID, maf.dir, SampleInfo.dir){
  #read maf file
  maf_input <- read.table(paste(maf.dir,'/',patientID,'.maf',sep = ""), quote = "", header = TRUE, fill = TRUE, sep = '\t')
  #read info file
  sample_info_input <-  read.table(SampleInfo.dir, quote = "", header = TRUE, fill = TRUE, sep = '')
  #Generate patient,lesion and time information in the last three columns of the maf file
  maf_input$patient=''
  maf_input$lesion=''
  maf_input$time=''
  maf_input$mut.id=''
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










