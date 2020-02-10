#' readMaf
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param mafFile tab delimited MAF file (plain text or *.gz compressed). The prefix to the filename of mafFile is considered as "patientID", which will be used for downstream analysis and visualization. Required. 
#' @param mutType Select proper variant classification you need. Default "All".Option: "nonSilent". 
#' @param mutNonSilent variant classifications which are considered as non-silent. Default NULL. Option: "Default". 
#' @param chrSilent Select chromosomes needed to be dismissed. Default NULL. 
#' @param use.indel logic. whether to use INDELs besides somatic SNVs. Default FALSE. 
#' @param ccfFile CCF file of SNVs. Default NULL. 
#' @param refBuild human reference genome versions of hg18, hg19 or hg38 by UCSC. Default "hg19". 
#' 
#' 
#' @examples
#' maf.File <- system.file("extdata/maf", "HCC6046.maf", package = "MesKit")
#' ccf.File <- system.file("extdata/", "HCC6046.CCF.txt", package = "MesKit")
#' maf <- readMaf(mafFile=maf.File, refBuild="hg19")
#' maf <- readMaf(mafFile=maf.File, ccfFile=ccf.File, refBuild="hg19")
#' @return an object of class Maf. 
#' 
#' @exportClass classMaf
#' @export readMaf


## read.maf main function
readMaf <- function(
    ## maf parameters
    mafFile,
    ## filter selection
    mutType="All", mutNonSilent="Default", chrSilent=NULL, use.indel=FALSE, 
    ## ccf parameters             
    ccfFile=NULL, 
    ## supplyment
    refBuild="hg19"
){
    
    ## read maf file from .maf or .gz file
    if (.substrRight(mafFile, 3) == ".gz"){
        mafData <- read.table(mafGz <- gzfile(mafFile, "r"), quote="", 
                               header=TRUE, fill=TRUE, 
                               sep='\t', stringsAsFactors=FALSE)
        
        close(mafGz)
    } else {
        mafData <- read.table(mafFile, quote="", 
                               header=TRUE, fill=TRUE, 
                               sep='\t', stringsAsFactors=FALSE)
    }
    
    ## if the filename is exactly the patientID
    fileName <- unlist(strsplit(mafFile, "/"))[length(unlist(strsplit(mafFile, "/")))]
    patientID <- strsplit(as.character(fileName), ".maf")[[1]][1]
    
                                  
    ## read ccf files
    if (!is.null(ccfFile)) {
        ccfInput <- read.table(ccfFile, quote="", 
                                      header=TRUE, fill=TRUE, 
                                      sep='\t', stringsAsFactors=FALSE)


        mafData <- uniteCCF(mafData, ccf)
    } 

    
    ## filter variant classification
    if (mutType == "nonSilent"){
        if (mutNonSilent == "Default"){
            nonSilent <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                           "Translation_Start_Site", "Nonsense_Mutation", 
                           "Nonstop_Mutation", "In_Frame_Del",
                           "In_Frame_Ins", "Missense_Mutation")
        } else {
            nonSilent <- mutNonSilent 
        }
        mafData <- mafData[which(mafData$Variant_Classification %in% nonSilent), ]
    } else if (mutType == "All"){
        # message("All variant classification submitted")
    } else {
        error("parameter `mut.type` error. 
              The mut.type should be either 'All' or 'nonSilent'. 
              You could further settle the filter by parameter 'mutNonSilent'.")
    }
    
    ## use.indel filter
    if(!use.indel){
        mafData <- mafData[which(mafData$Variant_Type == "SNP"),]
    }
    
    ## chromosome filter 
    if (!is.null(chrSilent)){
        mafData <- mafData[which(!mafData$Chromosome %in% chrSilent), ]
    }

    
    ## generate classMaf
    maf <- classMaf(data=data.table::setDT(mafData),  
                    patientID=patientID, 
                    ref.build=refBuild)
    
    # ## for parameter vafColumn="VAF", select particular VAF column
    # colnames(maf@data)[colnames(maf@data) == vafColumn] <- "VAF"
       
    return(maf)
}

.substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
}


##--- combine CCF into maf object
uniteCCF <- function(maf, ccf){
    maf <- tidyr::unite(
        maf, "mutID", c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Variant_Type"), 
        sep=":", remove = FALSE
        )
    ccf <- mutate(
        ccfInput, Variant_Type="SNP") %>%
        tidyr::unite("mutID", c("Sample_ID", "Chromosome", "Start_Position", "Variant_Type"), sep = ":", remove = FALSE) %>%
        dplyr::select(mutID, CCF, CCF_std)

    maf_merge_ccf <- merge(maf, ccf, by="mutID", all.x=TRUE)%>%
        dplyr::select(-mutID) 
}


##--- classMaf class
classMaf <- setClass(Class="classMaf", 
                     slots= c(data='data.table', patientID='character', 
                              ref.build='character'))
