#--- filter mutations in CNV regions 

copyNumberFilter <- function(maf, seg){
  ## combine data frame
  if(is(seg, "list")){
    seg <- plyr::rbind.fill(seg)
  }
  seg <- seg[!Chromosome %in% c("X","Y")]
  maf_data <- getMafData(maf)
  maf_data$ID <-  dplyr::select(tidyr::unite(maf_data, "ID", 
                                          Hugo_Symbol, Chromosome, 
                                          Start_Position, End_Position,
                                          Reference_Allele, Tumor_Seq_Allele2,Tumor_Sample_Barcode, 
                                          sep=":"), ID)
  seg$Chromosome <- as.character(seg$Chromosome)
  data.table::setkey(x = seg, Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position)
  sampleNames <- unique(seg[, Tumor_Sample_Barcode])
  sampleDat <- maf_data[Tumor_Sample_Barcode %in% sampleNames,]
  resID <- maf_data[!ID %in% sampleDat$ID]$ID
  overlapsDat <- data.table::foverlaps(x = sampleDat, y = seg, 
                                       by.x = c('Tumor_Sample_Barcode','Chromosome',
                                                'Start_Position', 'End_Position'))
  overlapsDat <-  overlapsDat[,.(Hugo_Symbol, Chromosome, i.Start_Position, i.End_Position,
                       Tumor_Sample_Barcode, VAF, Start_Position, End_Position, CopyNumber, Type, ID)] %>% 
                        dplyr::rename(Start_Position = i.Start_Position, End_Position = i.End_Position,
                                Segment_Start = Start_Position , Segment_End= End_Position)
  # colnames(overlapsDat)[c(3:4, 7:8)] <-  c('Start_Position', 'End_Position', 'Segment_Start', 'Segment_End')
  if(nrow(overlapsDat[is.na(overlapsDat$CopyNumber)]) > 0){
    message(paste('Removed ', nrow(overlapsDat[is.na(overlapsDat$CopyNumber)]), ' variants with no copy number data.', sep = ''))
    # print(overlapsDat[is.na(overlapsDat$CopyNumber)])
    overlapsDat = overlapsDat[!is.na(overlapsDat$CopyNumber)]
  }
  message(paste('Removed ', nrow(overlapsDat[overlapsDat$Type != "Neutral",]), ' variants in copy number altered regions.', sep = ''))
  # print(overlapsDat[overlapsDat$Type != "Neutral",])
  overlapsDat <-  overlapsDat[overlapsDat$Type == "Neutral",]
  maf_data <- maf_data[ID %in% c(overlapsDat$ID, resID)] %>% subset(select = -c(ID))
  maf@data <- maf_data
  return(maf)
}