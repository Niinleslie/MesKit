#' Subset Maf object
#' @description Read tab delimited MAF (can be plain text or *.gz compressed) file along with sample information file.
#'
#' @param geneList subset by geneList.Default:NULL
#' @param patient.id Select the specific patients. Default: NULL, all patients are included.
#' @param chrSilent Select chromosomes needed to be dismissed. Default NULL.
#' @param mutType select Proper variant classification you need. Default "All". Option: "nonSyn".
#' @param use.indel Logical value. Whether to use INDELs besides somatic SNVs. Default: TRUE.
#' @param min.vaf The minimum VAF for filtering variants. Default: 0.
#' @param max.vaf The maximum VAF for filtering variants. Default: 1.
#' @param min.average.vaf The minimum tumor average VAF for filtering variants. Default: 0.
#' @param min.average.adj.vaf The minimum tumor average ajust VAF for filtering variants. Default: 0.
#' @param min.ccf The minimum CCF for filtering variants. Default: 0.
#' @param min.ref.depth The minimum reference allele depth for filtering variants. Default: 0.
#' @param min.alt.depth The minimum alteratation allele depth for filtering variants. Default: 0.
#' @param min.total.depth The minimum total allele depth for filtering variants. Default: 0.
#' @param clonalStatus subset by clonal status.Default: NULL.Option: "Clonal","Subclonal".
#' @param use.adjVAF let VAF = VAF_adj, Tumor_Average_VAF = Tumor_Average_VAF_adj.Default: FALSE. 
#'
#'
#' @examples
#' maf <- readMaf(mafFile=maf.File, refBuild="hg19")
#' maf <- subsetMaf(maf)
#' @return Maf object.
#' @export

subsetMaf <- function(maf,
                      geneList = NULL,
                      patient.id = NULL,
                      chrSilent = NULL,
                      mutType = "All",
                      use.indel = TRUE,
                      min.vaf = 0,
                      max.vaf = 1,
                      min.average.vaf = 0,
                      min.average.adj.vaf = 0,
                      min.ccf = 0,
                      min.ref.depth = 0,
                      min.alt.depth = 0,
                      min.total.depth = 0,
                      clonalStatus = NULL,
                      use.adjVAF = FALSE){
   
   mafData <- maf@data
   nonSyn.vc <- maf@nonSyn.vc
   
   ## merge VAF,VAF_adj,CCF by tumor id
   ## calculate average VAF
   mafData <- mafData %>% 
      dplyr::group_by(Patient_ID,Tumor_ID,Chromosome,Start_Position,Reference_Allele,Tumor_Seq_Allele2) %>%
      dplyr::mutate(Total_allele_depth = Ref_allele_depth + Alt_allele_depth) %>% 
      dplyr::mutate(Tumor_Average_VAF = round(sum(VAF * Total_allele_depth)/sum(Total_allele_depth),3)) %>% 
      dplyr::filter(Tumor_Average_VAF >=min.average.vaf)
   
   ## calculate average adjust VAF
   if("VAF_adj" %in% colnames(mafData)){
      mafData <- mafData %>%
         dplyr::mutate(Tumor_Average_VAF_adj = round(sum(VAF_adj * Total_allele_depth)/sum(Total_allele_depth),3)) %>%
         dplyr::filter(Tumor_Average_VAF_adj>=min.average.adj.vaf)
   }
   
   if(use.adjVAF){
       if(!"VAF_adj" %in% colnames(mafData)){
           stop("Adjusted VAF was not found in maf data.")
       }
       else{
           mafData$VAF <- mafData$VAF_adj
           mafData$Tumor_Average_VAF <- mafData$Tumor_Average_VAF_adj
       }
   }
   
   ## calculate average adjust CCF
   if("CCF" %in% colnames(mafData)){
      mafData <- mafData %>% 
         dplyr::mutate(Tumor_Average_CCF = round(sum(CCF * Total_allele_depth)/sum(Total_allele_depth),3))
   }
   
   mafData <- mafData %>% 
      dplyr::ungroup() %>%
      dplyr::select(-Total_allele_depth) %>%  
      as.data.table()
   
   ## patient filter
   if(!is.null(patient.id)){
      patient.setdiff <- setdiff(patient.id, unique(mafData$Patient_ID))
      if(length(patient.setdiff) > 0){
         stop(paste0("Patient ", patient.setdiff, " can not be found in your data"))
      }
      mafData <- mafData[Patient_ID %in% patient.id]
   }
   
   ## chromosome filter
   if (!is.null(chrSilent)) {
      mafData <- mafData[mafData$Chromosome %in% chrSilent]
   }
   
   ## filter variant classification
   mutType <- match.arg(mutType, choices = c("All","nonSyn"), several.ok = FALSE)
   if (mutType == "nonSyn") {
      mafData <-  mafData[Variant_Classification %in% nonSyn.vc]
   }
   
   ## use.indel filter
   if (!use.indel) {
      mafData <- mafData[Variant_Type == "SNP"]
   }
   
   ## Select mutations in selected genes
   if(!is.null(geneList)){
      mafData <- mafData[Hugo_Symbol %in% geneList]
   }
   
   ## allele depth filter
   mafData <- mafData[Ref_allele_depth>=min.ref.depth&
                      Alt_allele_depth>=min.alt.depth &
                      Ref_allele_depth + Alt_allele_depth>=min.total.depth]
   mafData$Ref_allele_depth[is.na(mafData$Ref_allele_depth)] <- 0
   mafData$Alt_allele_depth[is.na(mafData$Alt_allele_depth)] <- 0
   
   ## vaf filter
   mafData <- mafData[VAF >= min.vaf & VAF <= max.vaf]
   
   ## ccf filter
   if("CCF" %in% colnames(mafData)){
      mafData <- mafData[CCF >= min.ccf]
   }
   
   if(!is.null(clonalStatus)){
      if(!"Clonal_Status" %in% colnames(mafData)){
         stop("Missing information about clonal status in maf data")
      }
      clonalStatus <- match.arg(clonalStatus, choices = c("Clonal", "Subclonal"))
      mafData <- mafData[Clonal_Status == clonalStatus]
   }
   
   
   maf@data <- mafData
   
   return(maf) 
}