
# extract recurrent CNA genes from gistic file
readGisticGene <- function(gisticGenesFile = NULL, Gistic.type = NULL){
   
   if(!is.null(gisticGenesFile)){
      gisticGenes <- data.table::fread(input = gisticGenesFile, stringsAsFactors = FALSE, header = TRUE)
      if(is.na(unique(gisticGenes[, ncol(gisticGenes), with = FALSE])[1])){     
         gisticGenes <- gisticGenes[,-ncol(gisticGenes), with = FALSE]
      }
      
      #wide peak boundaries
      wpb <- as.character(gisticGenes[3]) 
      wpb <- wpb[2:length(wpb)]
      
      # keep the information of peak boundaries and genes 
      gisticGenes <- gisticGenes[c(4:nrow(gisticGenes)),]
      gisticGenes$cytoband <- gsub(pattern = ' ', replacement = '_', x = gisticGenes$cytoband)
      colnames(gisticGenes) <- c( 'cytoband', paste0( colnames(gisticGenes)[2:length(colnames(gisticGenes))], '_',wpb))
      
      gisticGenes <- suppressWarnings(data.table::melt(gisticGenes, id.vars = 'cytoband'))
      gisticGenes <- gisticGenes[!value %in% '']
      
      #gisticGenes <- gisticGenes[!grep(pattern = '|', x = gisticGenes$value, fixed = TRUE)] #remove genes with ambiguous annotation
      gisticGenes <- gisticGenes[,.(variable, value)] %>%
         tidyr::separate(col = variable, into = c("Cytoband", "Wide_Peak_Boundaries"), 
                         sep = "_", remove = TRUE) %>%
         dplyr::mutate(Gene=sub("\\|.+", "", value), Gistic.type = Gistic.type) %>%
         dplyr::select(-value)
      
      return(gisticGenes)
   }
}

## read GISTIC all_lesions files
readGisticAllLesions <- function(gisticAllLesionsFile = NULL, verbose = TRUE){
   if(verbose){
      cat(paste0('--Processing ', basename(gisticAllLesionsFile), '\n'))
   }
   
   gisticLesions <- data.table::fread(input = gisticAllLesionsFile, stringsAsFactors = FALSE, header = TRUE)
   if(is.na(unique(gisticLesions[, ncol(gisticLesions), with = FALSE])[1])){     
      gisticLesions <- gisticLesions[,-ncol(gisticLesions), with = FALSE]
   }
   
   gisticLesions <- gisticLesions %>%
      dplyr::select(1,2,3,6) %>%
      `colnames<-` (c("PeakID", "Cytoband", "WPB", "Qvalue")) %>%
      dplyr::filter(!grepl("values", PeakID)) %>%
      dplyr::mutate(WPB = sapply(strsplit(.[, "WPB"], split = "(", fixed =  TRUE), '[', 1)) %>%
      dplyr::mutate(Gistic.type = substr(PeakID, start = 1, stop = 3), PeakID = NULL)
   
   return(gisticLesions)
   
}