# sub-functions to classify mutations

# type = "SP"(shared pattern)
mutType_SP <- function(
    #Tumor_ID, 
    unique_barcode_count, 
    total_barcode_count,
    unique_tumor_count,
    tumor_barcode_ID, 
    total_tumor_count,
    TypeCount,
    classByTumor) {
    
    if(classByTumor){
        dplyr::case_when(
            #unique_barcode_count == 1 ~ "Private",
            (unique_tumor_count == total_tumor_count) ~ "Public",  
            (TypeCount == "Single") & (unique_tumor_count == 1) ~ "Private",
            (TypeCount == "Multi") & (unique_tumor_count == 1) ~ paste0("Private_", tumor_barcode_ID),
            (TypeCount == "Single") & (1 < unique_tumor_count) &
                (unique_tumor_count < total_tumor_count) ~ "Shared",
            (TypeCount == "Multi") & (1 < unique_tumor_count) &
                (unique_tumor_count < total_tumor_count) ~ paste0("Shared_", tumor_barcode_ID)   
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
    #Tumor_ID, 
    unique_barcode_count, 
    total_barcode_count,
    unique_tumor_count,
    tumor_barcode_ID, 
    total_tumor_count,
    Clonal_Status, 
    TypeCount,
    classByTumor) {
    
    if(classByTumor){
        dplyr::case_when(
            unique_tumor_count == total_tumor_count ~ paste0("Public_", Clonal_Status),        
            (TypeCount == "Single") & (unique_tumor_count == 1)~ paste0("Private_", Clonal_Status),
            (TypeCount == "Multi") & (unique_tumor_count == 1)~ paste0("Private_", tumor_barcode_ID, "_", Clonal_Status),
            (TypeCount == "Single") & (1 < unique_tumor_count) &
                (unique_tumor_count < total_tumor_count) ~ paste0("Shared_", Clonal_Status),
            (TypeCount == "Multi") & (1 < unique_tumor_count) & (unique_tumor_count < total_tumor_count)
            ~ paste0("Shared_", tumor_barcode_ID,"_", Clonal_Status)   
        )
    } else{
        dplyr::case_when(
            unique_barcode_count == 1 ~ paste0("Private_", Clonal_Status),
            unique_barcode_count == total_barcode_count ~ paste0("Public_", Clonal_Status),
            (1 < unique_barcode_count) &
                (unique_barcode_count < total_barcode_count) ~ paste0("Shared_", Clonal_Status)
        )       
    }   
}


multiHits <- function(x) {
    if (length(x) > 1) {
        if (length(unique(x)) == 1){
            l <- paste0(x, collapse = ";")
        }else{
            l <- paste0(paste0(x, collapse = ";"), ";Multi_hits")
        }
        return(l)
    }
    else {
        return(x)
    }
}


do.classify <- function(
    maf_data,
    class = "SP",
    classByTumor = FALSE) {
    
    class.options = c('SP', 'CS', 'SPCS')
    if(!class %in% class.options){
        stop("Class can only be either 'SP', 'CS' or 'SPCS'")
    }
    
    maf_data <- maf_data %>%
        tidyr::unite(
            "mutation_id",
            c("Chromosome",
              "Start_Position",
              "Reference_Allele",
              "Tumor_Seq_Allele2"
            ),
            sep = ":",
            remove = FALSE
        )
    
    mutation_barcode_count <-
        maf_data %>%
        dplyr::group_by(.data$Patient_ID, .data$mutation_id) %>%
        #dplyr::group_by(., Patient_ID) %>%
        #{if(classByTumor) dplyr::group_by(., Patient_ID, Tumor_ID, mutation_id)
        #else dplyr::group_by(., Patient_ID, mutation_id)
        #} %>%
        dplyr::summarise(unique_barcode_count = dplyr::n_distinct(.data$Tumor_Sample_Barcode))
    
    patient_barcode_count <-
        maf_data %>%
        dplyr::group_by(.data$Patient_ID) %>%
        #{if(classByTumor) dplyr::group_by(., Patient_ID, Tumor_ID)
        #else dplyr::group_by(., Patient_ID)
        #} %>%
        dplyr::summarise(total_barcode_count = dplyr::n_distinct(.data$Tumor_Sample_Barcode),
                         total_tumor_count = dplyr::n_distinct(.data$Tumor_ID))
    
    
    
    unique_tumor_count <-maf_data %>%
        #dplyr::group_by(., Patient_ID, Tumor_ID, mutation_id) %>%
        #dplyr::summarise(unique_tumor_count = dplyr::n_distinct(Tumor_Sample_Barcode))
        dplyr::group_by(.data$Patient_ID, .data$mutation_id) %>%
        dplyr::summarise(
            unique_tumor_count = dplyr::n_distinct(.data$Tumor_ID), 
            tumor_barcode_ID = paste(unique(.data$Tumor_ID), collapse = "_")
        )

    maf_data <-
        suppressMessages(maf_data %>%
                             dplyr::left_join(mutation_barcode_count) %>%
                             dplyr::left_join(patient_barcode_count) %>%
                             dplyr::left_join(unique_tumor_count)
        )
    
    if(length(unique(maf_data$Tumor_ID)) == 1){
        TypeCount = "Single"
    }else{TypeCount = "Multi"}
    
    if(class == "SP"){
        maf_data <- maf_data %>%
            dplyr::mutate(
                Mutation_Type = mutType_SP(
                    #Tumor_ID, 
                    unique_barcode_count =  .data$unique_barcode_count, 
                    total_barcode_count =  .data$total_barcode_count,
                    unique_tumor_count =  .data$unique_tumor_count,
                    tumor_barcode_ID =  .data$tumor_barcode_ID, 
                    total_tumor_count =  .data$total_tumor_count,
                    TypeCount = TypeCount,
                    classByTumor)
            )
    }else{
        if(! "Clonal_Status" %in% colnames(maf_data)){
            stop(paste0("Clonal status could not be identified without CCF and CCF_Std data!"))
        }
        else{
            if(class == "CS"){
                maf_data <- maf_data %>%
                    dplyr::mutate(
                        Mutation_Type = .data$Clonal_Status
                    )
            }
            else if(class == "SPCS"){
                maf_data <- maf_data %>%
                    dplyr::filter(!is.na(.data$Clonal_Status)) %>%
                    dplyr::mutate(Mutation_Type =
                                      mutType_SPCS(
                                          #Tumor_ID, 
                                          unique_barcode_count =  .data$unique_barcode_count, 
                                          total_barcode_count =  .data$total_barcode_count,
                                          unique_tumor_count =  .data$unique_tumor_count,
                                          tumor_barcode_ID =  .data$tumor_barcode_ID,                                          
                                          total_tumor_count =  .data$total_tumor_count,
                                          Clonal_Status =  .data$Clonal_Status,
                                          TypeCount =  TypeCount,
                                          classByTumor =  classByTumor)
                    )
            }
        }
    }
    
    maf_data <- maf_data %>%
        dplyr::select("Hugo_Symbol", "Chromosome", 
                      "Start_Position", "End_Position",
                      "Reference_Allele", "Tumor_Seq_Allele2", 
                      "Tumor_Sample_Barcode", "Mutation_Type",
                      "unique_barcode_count", "Patient_ID"
                      #type_barcode_ID,unique_tumor_count,total_barcode_count
        )
    
    
    return(maf_data)
}
