#--- filter mutations in CNV regions 

copyNumberFilter <- function(maf_data, seg, use.tumorSampleLabel = FALSE){
  ## combine data frame
  if(is(seg, "list")){
    seg <- dplyr::bind_rows(seg) %>%
      dplyr::filter(Patient_ID == unique(maf_data$Patient_ID)) %>% 
      as.data.table()
  }
  
  # if("LOH" %in% colnames(seg)){
  #   seg <- seg[seg$LOH == FALSE,]
  #   # message("Remove segment with LOH" )
  # }
  
  if(use.tumorSampleLabel){
    seg$Tumor_Sample_Barcode <- seg$Tumor_Sample_Label
  }
  
  seg <- seg[!seg$Chromosome %in% c("X","Y")]
  maf_data$mut_id <-  dplyr::select(tidyr::unite(maf_data, "mut_id", 
                                          "Hugo_Symbol", "Chromosome", 
                                          "Start_Position", "End_Position",
                                          "Reference_Allele", "Tumor_Seq_Allele2",
                                          "Tumor_Sample_Barcode", 
                                          sep=":"), "mut_id")
  seg$Chromosome <- as.character(seg$Chromosome)
  data.table::setkey(x = seg, "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position")
  sampleNames <- unique(seg$Tumor_Sample_Barcode)
  sampleDat <- maf_data[maf_data$Tumor_Sample_Barcode %in% sampleNames,]
  resID <- maf_data[!maf_data$mut_id %in% sampleDat$mut_id]$mut_id
  overlapsDat <- data.table::foverlaps(x = sampleDat, y = seg, 
                                       by.x = c('Tumor_Sample_Barcode','Chromosome',
                                                'Start_Position', 'End_Position'))
  
  col_keep <- c(
    "Hugo_Symbol", 
    "Chromosome", 
    "i.Start_Position",
    "i.End_Position",
    "Tumor_Sample_Barcode",
    "VAF", 
    "Start_Position",
    "End_Position",
    "CopyNumber", 
    "Type", 
    "mut_id"
  )
  if("LOH" %in% colnames(seg)){
    col_keep <- c(col_keep, "LOH")
  }
  
  overlapsDat <-  overlapsDat %>% 
    dplyr::select(
      all_of(col_keep)
    ) %>% 
    dplyr::rename(
      "Start_Position" = "i.Start_Position", 
      "End_Position" = "i.End_Position",
      "Segment_Start" = "Start_Position",
      "Segment_End"= "End_Position"
    )
  num_all <- nrow(overlapsDat)
  
  # remove mutations within copy-number altered regions 
  overlapsDat <-  overlapsDat[c(overlapsDat$Type == "Neutral"|is.na(overlapsDat$Type)), ]
  
  if("LOH" %in% colnames(seg)){
    mes <- 'within copy-number altered or LOH regions.'
    # remove mutations within copy-number altered regions 
    overlapsDat <- overlapsDat[overlapsDat$LOH == "FALSE"|is.na(overlapsDat$Type),]
  }else{
    mes <- 'within copy-number altered regions.'
  }
  num_filt <- num_all - nrow(overlapsDat)
  
  message(paste('Removed ', num_filt,' ', mes, sep = ''))
  
  maf_data <- maf_data[maf_data$mut_id %in% c(overlapsDat$mut_id, resID)] %>% 
    dplyr::select(
      -"mut_id"
    )
  
  return(maf_data)
}