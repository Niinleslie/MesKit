# sub-functions to classify mutations

# type = "SP"(shared pattern)
mutType_SP <- function(
  #Tumor_Type, 
  unique_barcode_count, 
  total_barcode_count,
  type_barcode_count,
  type_barcode_ID, 
  TypeCount,
  classByType) {

  if(classByType){
    dplyr::case_when(
        unique_barcode_count == 1 ~ "Private",
        #(TypeCount == "Single") & (unique_barcode_count == 1) ~ "Private",
        #(TypeCount == "Multi") & (unique_barcode_count == 1) ~ paste0("Private_", Tumor_Type),
        (unique_barcode_count == total_barcode_count) ~ "Public",        
        (TypeCount == "Single") & (1 < type_barcode_count) &
            (type_barcode_count < total_barcode_count) ~ "Shared",
        (TypeCount == "Multi") & (1 < type_barcode_count) &
        (type_barcode_count < total_barcode_count) ~ paste0("Shared_", type_barcode_ID)   
    )
  }
  else{
    dplyr::case_when(
        unique_barcode_count == 1 ~ "Private",
        unique_barcode_count == total_barcode_count ~ "Public",
        (1 < unique_barcode_count) &
            (unique_barcode_count < total_barcode_count) ~ "Shared"
    )
  }
}

# type = "SPCS"(shared pattern & clonal/subclonal)
mutType_SPCS <- function(
  #Tumor_Type, 
  unique_barcode_count, 
  total_barcode_count,
  type_barcode_count,
  type_barcode_ID, 
  Clonal_Status, 
  TypeCount,
  classByType) {

  if(classByType){
    dplyr::case_when(
      #TypeCount == "Single" & unique_barcode_count == 1 ~ paste0("Private_", Clonal_Status),
      #TypeCount == "Multi" & unique_barcode_count == 1 ~ paste0("Private_", "_", Clonal_Status),      
      unique_barcode_count == 1 ~ paste0("Private_", Clonal_Status),
      unique_barcode_count == total_barcode_count ~ "Public",        
      (TypeCount == "Single") & (1 < type_barcode_count) &
          (type_barcode_count < total_barcode_count) ~ paste0("Shared_", Clonal_Status),
      TypeCount == "Multi" & (1 < type_barcode_count) &
      (unique_barcode_count < type_barcode_count)
        ~ paste0("Shared_", type_barcode_ID, Clonal_Status)   
    )
  } else{
    dplyr::case_when(
        unique_barcode_count == 1 ~ paste0("Private_", Clonal_Status),
        unique_barcode_count == total_barcode_count ~ paste0("Shared_", Clonal_Status),
        (1 < unique_barcode_count) &
            (unique_barcode_count < total_barcode_count) ~ paste0("P_shared_", Clonal_Status)
    )       
  }   
}


multiHits <- function(x) {
    if (length(x) > 1) {
        return(paste0(paste0(x, collapse = ";"), ";Multi_hits"))
    }
    else {
        return(x)
    }
}
