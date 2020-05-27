# sub-functions to classify mutations

# type = "SP"(shared pattern)
mutType_SP <- function(
    #Tumor_ID, 
    unique_barcode_count, 
    total_barcode_count,
    type_barcode_count,
    tumor_barcode_ID, 
    TypeCount,
    classByTumor) {
    
    if(classByTumor){
        dplyr::case_when(
            #unique_barcode_count == 1 ~ "Private",
            (TypeCount == "Single") & (unique_barcode_count == 1) ~ "Private",
            (TypeCount == "Multi") & (unique_barcode_count == 1) ~ paste0("Private_", tumor_barcode_ID),
            (unique_barcode_count == total_barcode_count) ~ "Public",        
            (TypeCount == "Single") & (1 < type_barcode_count) &
                (type_barcode_count < total_barcode_count) ~ "Shared",
            (TypeCount == "Multi") & (1 < type_barcode_count) &
                (type_barcode_count < total_barcode_count) ~ paste0("Shared_", tumor_barcode_ID)   
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
    type_barcode_count,
    tumor_barcode_ID, 
    Clonal_Status, 
    TypeCount,
    classByTumor) {
    
    if(classByTumor){
        dplyr::case_when(
            #TypeCount == "Single" & unique_barcode_count == 1 ~ paste0("Private_", Clonal_Status),
            #TypeCount == "Multi" & unique_barcode_count == 1 ~ paste0("Private_", "_", Clonal_Status),      
            unique_barcode_count == 1 ~ paste0("Private_", Clonal_Status),
            unique_barcode_count == total_barcode_count ~ paste0("Public_", Clonal_Status),        
            (TypeCount == "Single") & (1 < type_barcode_count) &
                (type_barcode_count < total_barcode_count) ~ paste0("Shared_", Clonal_Status),
            TypeCount == "Multi" & (1 < type_barcode_count) &
                (unique_barcode_count < type_barcode_count)
            ~ paste0("Shared_", tumor_barcode_ID, Clonal_Status)   
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
    maf,
    class = "SP",
    patient.id = NULL,
    classByTumor = FALSE) {
    
    class.options = c('SP', 'CS', 'SPCS')
    if(!class %in% class.options){
        stop("Error:class can only be either 'SP', 'CS' or 'SPCS'")
    }
    
    
    if(is.null(patient.id)){
        patient.id = unique(maf@data$Patient_ID)
    }else{
        patient.setdiff <- setdiff(patient.id, unique(maf@data$Patient_ID))
        if(length(patient.setdiff) > 0){
            stop(paste0(patient.setdiff, " can not be found in your data"))
        }
    }
    
    maf_data <- maf@data %>%
        tidyr::unite(
            "mutation_id",
            c("Chromosome",
              "Start_Position",
              "Reference_Allele",
              "Tumor_Seq_Allele2"
            ),
            sep = ":",
            remove = FALSE
        ) %>%
        dplyr::filter(Patient_ID %in% patient.id)
    
    mutation_barcode_count <-
        maf_data %>%
        dplyr::group_by(., Patient_ID, mutation_id) %>%
        #dplyr::group_by(., Patient_ID) %>%
        #{if(classByTumor) dplyr::group_by(., Patient_ID, Tumor_ID, mutation_id)
        #else dplyr::group_by(., Patient_ID, mutation_id)
        #} %>%
        dplyr::summarise(unique_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode))
    
    patient_barcode_count <-
        maf_data %>%
        dplyr::group_by(., Patient_ID) %>%
        #{if(classByTumor) dplyr::group_by(., Patient_ID, Tumor_ID)
        #else dplyr::group_by(., Patient_ID)
        #} %>%
        dplyr::summarise(total_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode))
    
    
    
    type_barcode_count <-maf_data %>%
        #dplyr::group_by(., Patient_ID, Tumor_ID, mutation_id) %>%
        #dplyr::summarise(type_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode))
        dplyr::group_by(., Patient_ID, mutation_id) %>%
        dplyr::summarise(
            type_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode), 
            tumor_barcode_ID = paste(unique(Tumor_ID), collapse = "_")
        )
    
    maf_data <-
        suppressMessages(maf_data %>%
                             dplyr::left_join(mutation_barcode_count) %>%
                             dplyr::left_join(patient_barcode_count) %>%
                             dplyr::left_join(type_barcode_count)
        )
    
    if(length(unique(maf_data$Tumor_ID)) == 1){
        TypeCount = "Single"
    }else{TypeCount = "Multi"}
    
    if(class == "SP"){
        maf_data <- maf_data %>%
            dplyr::mutate(
                Mutation_Type = mutType_SP(
                    #Tumor_ID, 
                    unique_barcode_count, 
                    total_barcode_count,
                    type_barcode_count,
                    tumor_barcode_ID, 
                    TypeCount,
                    classByTumor)
            )
    }else{
        if(! "Clonal_Status" %in% colnames(maf_data)){
            stop(paste0("Clonal status could not be identified without CCF data!"))
        }
        else{
            if(class == "CS"){
                maf_data <- maf_data %>%
                    dplyr::mutate(
                        Mutation_Type = Clonal_Status
                    )
            }
            else if(class == "SPCS"){
                maf_data <- maf_data %>%
                    dplyr::filter(!is.na(Clonal_Status)) %>%
                    dplyr::mutate(Mutation_Type =
                                      mutType_SPCS(
                                          #Tumor_ID, 
                                          unique_barcode_count, 
                                          total_barcode_count,
                                          type_barcode_count,
                                          tumor_barcode_ID, 
                                          Clonal_Status,
                                          TypeCount,
                                          classByTumor)
                    )
            }
        }
    }
    
    maf_data <- maf_data %>%
        dplyr::select(Hugo_Symbol, Chromosome, 
                      Start_Position, End_Position,
                      Reference_Allele, Tumor_Seq_Allele2, 
                      Tumor_Sample_Barcode, Mutation_Type,
                      unique_barcode_count, Patient_ID
                      #type_barcode_ID,type_barcode_count,total_barcode_count
        )
    
    
    return(maf_data)
}
