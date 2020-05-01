do.classify <- function(
                   maf,
                   class = "SP",
                   patient.id = NULL,
                   classByType = FALSE) {

  class.options = c('SP', 'CS', 'SPCS')
  if(!class %in% class.options){
    stop("class can only be either 'SP', 'CS' or 'SPCS'")
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
    #{if(classByType) dplyr::group_by(., Patient_ID, Tumor_Type, mutation_id)
    #else dplyr::group_by(., Patient_ID, mutation_id)
    #} %>%
    dplyr::summarise(unique_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode))
  
  patient_barcode_count <-
    maf_data %>%
    dplyr::group_by(., Patient_ID) %>%
    #{if(classByType) dplyr::group_by(., Patient_ID, Tumor_Type)
    #else dplyr::group_by(., Patient_ID)
    #} %>%
    dplyr::summarise(total_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode))
  
  
  
  type_barcode_count <-maf_data %>%
    #dplyr::group_by(., Patient_ID, Tumor_Type, mutation_id) %>%
    #dplyr::summarise(type_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode))
    dplyr::group_by(., Patient_ID, mutation_id) %>%
    dplyr::summarise(
      type_barcode_count = dplyr::n_distinct(Tumor_Sample_Barcode), 
      type_barcode_ID = paste(unique(Tumor_Type), collapse = "_")
    )
  
  maf_data <-
    suppressMessages(maf_data %>%
                       dplyr::left_join(mutation_barcode_count) %>%
                       dplyr::left_join(patient_barcode_count) %>%
                       dplyr::left_join(type_barcode_count)
    )
  
  if(length(unique(maf_data$Tumor_Type)) == 1){
    TypeCount = "Single"
  }else{TypeCount = "Multi"}
  
  if(class == "SP"){
    maf_data <- maf_data %>%
      dplyr::mutate(
        Mutation_Type = mutType_SP(
          #Tumor_Type, 
          unique_barcode_count, 
          total_barcode_count,
          type_barcode_count,
          type_barcode_ID, 
          TypeCount,
          classByType)
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
                            #Tumor_Type, 
                            unique_barcode_count, 
                            total_barcode_count,
                            type_barcode_count,
                            type_barcode_ID, 
                            Clonal_Status,
                            TypeCount,
                            classByType)
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


genHeatmapPlotMat <- function(
        maf,
        patient.id = NULL,
        class = "SP",
        classByType = FALSE, 
        topGenesCount = 15, 
        geneList = NULL) {
     
    maf_data <- do.classify(maf, classByType = classByType, patient.id = patient.id, class = class)
  
    if (!is.null(geneList)) {  
  
      geneSelect <- geneList
  
      maf_data <- maf_data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(Selected_Mut = dplyr::if_else(
          any(Hugo_Symbol %in% geneSelect),
          TRUE,
          FALSE))
    
      maf_data <- maf_data %>%
          dplyr::filter(Selected_Mut)
    }
 
    col_labels <- dplyr::select(maf_data, Patient_ID, Tumor_Sample_Barcode)%>%
                    distinct(.)
    col_labels <- as.vector(col_labels$Tumor_Sample_Barcode)

    # split by patient
    patient.split <- maf_data %>%
        dplyr::select(Patient_ID, Tumor_Sample_Barcode) %>%
        dplyr::distinct() %>%
        dplyr::select(Patient_ID) %>%
        as.matrix() %>%
        as.vector() %>%
        as.character()

    if(length(unique(patient.split)) == 1){
        patient.split = NULL
    }

    # long -> wider
    #col_labels <- unique(maf$Tumor_Sample_Barcode)
    mat <- maf_data %>%
        dplyr::ungroup() %>%
        dplyr::group_by(Hugo_Symbol) %>%
        dplyr::mutate(
            total_barcode_count = sum(unique_barcode_count)
            ) %>%
        dplyr::select(Hugo_Symbol,
                      Patient_ID,
                      Tumor_Sample_Barcode,
                      Mutation_Type,
                      total_barcode_count
                      ) %>%
        tidyr::pivot_wider(
            #names_from = Tumor_Sample_Barcode,
            names_from = c(Patient_ID, Tumor_Sample_Barcode),
            values_from = Mutation_Type,
            values_fn = list(Mutation_Type = multiHits)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(total_barcode_count)) %>%       
        #dplyr::select_if(function(x) {!all(is.na(x))}) %>%
        dplyr::slice(1:topGenesCount) %>%
        tibble::column_to_rownames(., "Hugo_Symbol") %>% 
        dplyr::select(-total_barcode_count) %>%
        as.matrix()
    
    
    return(list(mat, patient.split, col_labels))
}
