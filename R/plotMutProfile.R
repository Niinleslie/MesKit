#' @import ComplexHeatmap

genHeatmapPlotMatrix <- function(
        maf_data, 
        topGenesCount = 15) {

    if ("Selected_Mut" %in% colnames(maf_data)) {
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

plotMutProfile <- function(maf_data,
                           class = "SP",
                           topGenesCount = 15,
                           classByType = FALSE,
                           bgCol = "#f0f0f0",
                           patientsCol = NULL,
                           remove_empty_columns = TRUE,
                           remove_empty_rows = TRUE, 
                           showColnames = TRUE) {
        
    maf.plot <- genHeatmapPlotMatrix(maf_data, topGenesCount = topGenesCount)  
    mat <- maf.plot[[1]]

    #col_labels <- dplyr::select(maf_data, Patient_ID, Tumor_Sample_Barcode)%>%
                    #distinct(.)
    #col_labels <- as.vector(col_labels$Tumor_Sample_Barcode)
    patient.split <- maf.plot[[2]]
    col_labels <- maf.plot[[3]]

    # get the order or rows
    stat <- rep(0, topGenesCount)
    for(i in 1:nrow(mat)){
      stat[i] <- sum(!is.na(mat[i, ])) / ncol(mat)
    }
    
    rowOrderFrame <- data.frame(Genes = rownames(mat), freq = stat)
    rowOrder <- as.numeric(rownames(rowOrderFrame[order(rowOrderFrame$freq, decreasing = TRUE), ]))
    
    # View(mat)

    #patient_id_cols <-
        #RColorBrewer::brewer.pal(length(unique(patient.split)), "Set")
    #names(patient_id_cols) <- unique(patient.split)
#
    #patient_barcode <- maf_data %>%
        #dplyr::select(Patient_ID, Tumor_Sample_Barcode) %>%
        #dplyr::distinct() %>%
        #dplyr::mutate(color = patient_id_cols[patient.split]) 
#
    #sample_barcode <- patient_barcode$color
    #names(sample_barcode) <- patient_barcode$Tumor_Sample_Barcode


    multi_hit_exist = FALSE
    for (i in 1:nrow(mat)) {
        for (j in 1:ncol(mat)) {
            if (length(temp <- grep("Multi_hits", mat[i, j]))) {
                multi_hit_exist = TRUE
                break
            }
        }
    }
        
    if (!(classByType)) {
      col_type <- function(class) {
          if (class == "SP") {
              cols <- c("#3C5488FF", "#00A087FF", "#F39B7fFF")
              names(cols) <- c("Public","Shared", "Private")
          } else if (class == "CS") {
              cols <- c("#00A087FF", "#3C5488FF")
              names(cols) <- c("Clonal", "Subclonal")
          } else if (class == "SPCS") {
              cols <-
                  c(
                      "#00A087FF",
                      "#3C5488FF",
                      "#8491B4FF",
                      "#F39B7FFF", 
                      "#E64B35FF",                    
                      "#4DBBD5FF"                    
                  )
              names(cols) <-
                  c(
                      "Public_Clonal",
                      "Public_Subclonal",
                      "Shared_Clonal",
                      "Shared_Subclonal",
                      "Private_Clonal",
                      "Private_Subclonal"                    
                  )
          }
        
          return(cols)
      }

      alter_fun <- function(class){
          if(class == "SP"){
              l <- list(
              background = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = bgCol, col = NA)),
              Private = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SP")["Private"], col = NA)),
              Public = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SP")["Public"], col = NA)),
              Shared = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SP")["Shared"], col = NA)),
              Multi_hits = function(x, y, w, h)
                  grid::grid.points(x, y, pch = 16, size = grid::unit(0.5, "char")
              ))
          }else if(class == "CS"){
              l <- list(
              background = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = bgCol, col = NA)),
              Clonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("CS")["Clonal"], col = NA)),
              Subclonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("CS")["Subclonal"], col = NA)),
              Multi_hits = function(x, y, w, h)
                  grid::grid.points(x, y, pch = 16, size = grid::unit(0.5, "char") 
                  ))
          }else if(class == "SPCS" ){
              l <- list(
              background = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = bgCol, col = NA)),
              Private_Clonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SPCS")["Private_Clonal"], col = NA)),
              Private_Subclonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SPCS")["Private_Subclonal"], col = NA)),
              Shared_Clonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SPCS")["Public_Clonal"], col = NA)),
              Shared_Subclonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SPCS")["Public_Subclonal"], col = NA)),
              P_shared_Clonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SPCS")["Shared_Clonal"], col = NA)),
              P_shared_Subclonal = function(x, y, w, h)
                  grid::grid.rect(x, y, w * 0.9, h * 0.9,
                                  gp = grid::gpar(fill = col_type("SPCS")["Shared_Subclonal"], col = NA)),
              Multi_hits = function(x, y, w, h)
                  grid::grid.points(x, y, pch = 16, size = grid::unit(0.5, "char") 
              ))
          }
          return(l)            
      }
     
    } else{
      
      # set certain colors
      ColorScale <- c("#E41A1C","#377EB8","#7F0000",
                      "#35978f","#FC8D62","#2166ac",
                      "#E78AC3","#A6D854","#FFD92F",
                      "#E5C494","#8DD3C7", "#6E016B" ,
                      "#BEBADA", "#e08214", "#80B1D3",
                      "#d6604d","#ffff99","#FCCDE5",
                      "#FF6A5A","#BC80BD","#CCEBC5" ,
                      "#fb9a99","#B6646A", "#9F994E", 
                      "#7570B3" ,"#c51b7d" ,"#66A61E" ,
                      "#E6AB02" ,"#003c30", "#666666")
      
      mutationTypes <- na.omit(unique(maf_data$Mutation_Type))
      
      # filter mutation types
      filteredTypes <- c()
      for (i in 1:length(mutationTypes)){
          if (length(grep(mutationTypes[i], mat)) != 0) {
              filteredTypes <- c(filteredTypes, mutationTypes[i])
          }
      }
      mutationTypes <- filteredTypes
      
      
      # sort types in legend
      if (class == "SP" | class == "SPCS") {
        sortType <- function(types) {
            publicType <- types[grep("Public", types)]
            sharedType <- types[grep("Shared", types)]
            privateType <- types[grep("Private", types)]
            return(c(publicType, sharedType, privateType))
        }
      
        mutationTypes <- sortType(mutationTypes)
      }
      
      col_type <- function(class) {
        set.seed(123)
        cols <- sample(ColorScale, size = length(mutationTypes), replace = FALSE)
        names(cols) <- mutationTypes
        return(cols)
        
      }
      
      # prepare functions for assignment
      alter_fun_functions <- list()
      for (type_num in 1:length(mutationTypes)){
        alter_fun_function <- paste0("function(x, y, w, h) grid::grid.rect(x, y, w * 0.9, h * 0.9,
                              gp = grid::gpar(fill = col_type(\'",class,"\')[\'",mutationTypes[type_num], "\'],col = NA))")
        alter_fun_functions <- c(alter_fun_functions, eval(parse(text = alter_fun_function)))
      }
      names(alter_fun_functions) <- mutationTypes
      
      alter_fun <- function(class){
        l <- c(alter_fun_functions, Multi_hits = function(x, y, w, h)
          grid::grid.points(x, y, pch = 16, size = grid::unit(0.5, "char") 
          ), background = function(x, y, w, h)
            grid::grid.rect(x, y, w * 0.9, h * 0.9,
                            gp = grid::gpar(fill = bgCol, col = NA)))
        return(l)            
      }
    }
    # prepare legends
    
    ## type legend
    
    heatmapLegend <- ComplexHeatmap::Legend(title = "Type", 
                            title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                            at = names(col_type(class)),
                            labels = sub("_", "-", names(col_type(class))),
                            labels_gp = grid::gpar(fontsize = 12),
                            grid_width = unit(4, "mm"),
                            grid_height = unit(4, "mm"), legend_gp = grid::gpar(fill = col_type(class)))
    
    ## patient legend
    patient.id <- unique(patient.split)
    
    set.seed(1234)
    patientsCol <- sample(colors(), length(patient.id), replace = FALSE)    
    names(patientsCol) <- patient.id
    
    patientLegend <-  ComplexHeatmap::Legend(
                                 labels = patient.id, 
                                 legend_gp = grid::gpar(fill = patientsCol), 
                                 title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                                 labels_gp = grid::gpar(fontsize = 12),
                                 grid_width = unit(4, "mm"),
                                 grid_height = unit(4, "mm"), title = "Patient")
    
    ## multi-hits legend
    multiLegend <- ComplexHeatmap::Legend(
                          labels = "Multi_hits",
                          labels_gp = grid::gpar(fontsize = 12),
                          type = "points",
                          pch = 16,
                          grid_width = unit(4, "mm"),
                          grid_height = unit(4, "mm")
                   )
    
    
    ## type-multi legend
    hm <- ComplexHeatmap::packLegend(heatmapLegend, multiLegend, direction = "vertical", gap = unit(0.3, "mm"))
    
    ## type-multi-patient legend
    hmp <- ComplexHeatmap::packLegend(hm, patientLegend, direction = "vertical", gap = unit(0.3, "cm"))
    
    ## type-patient legend
    hp <- ComplexHeatmap::packLegend(heatmapLegend, patientLegend, direction = "vertical", gap = unit(1.2, "cm"))
    
    
    
    ht <- suppressMessages(
        ComplexHeatmap::oncoPrint(
            mat,
            alter_fun = alter_fun(class),
            col = col_type(class),
            column_title = "Mutational profile",
            column_title_gp = grid::gpar(fontsize = 16, fontface = "bold", col = "black"),
            row_title_gp = grid::gpar(fontsize = 11, fontface = "plain", col = "black"),
            #heatmap_legend_param = heatmap_legend(class),
            show_heatmap_legend = FALSE,
            remove_empty_columns = remove_empty_columns,
            remove_empty_rows = remove_empty_rows,
            row_order = rowOrder,
            row_names_gp = grid::gpar(fontsize = 11, fontface = "italic", col = "black"),
            column_names_gp = grid::gpar(fontsize = 11, fontface = "plain", col = "black"),
            pct_digits = 2,
            pct_side = "right",
            row_names_side = "left", 
            column_split = patient.split,
            column_order = colnames(mat),
            column_labels = col_labels,
            show_column_names = showColnames,
            bottom_annotation = if(
                is.null(patient.split)) NULL else{
                ComplexHeatmap::HeatmapAnnotation(
                #df = data.frame(patient = colnames(mat)),
                df = data.frame(Patient = patient.split),
                show_annotation_name = FALSE,
                col = list(Patient = patientsCol),
                simple_anno_size = unit(0.2, "cm"),
                show_legend = FALSE,
                annotation_legend_param = list(title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                                               labels_gp = grid::gpar(fontsize = 10),
                                               grid_width = unit(3.5, "mm"),
                                               grid_height = unit(3.5, "mm")
                                               #plot = FALSE
                                               )
                
                )}                
        )
    )



    if (multi_hit_exist) {
      ComplexHeatmap::draw(ht, heatmap_legend_list = hmp,
                           padding = unit(c(3, 3, 3, 3), "mm"))
            #annotation_legend_list = if(
                #is.null(patient.split)) NULL else{
                #list(ComplexHeatmap::Legend(
                #title = "patient",
                #title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                #at = patient.split,
                ##labels = names(patient_id_cols),
                ##legend_gp = grid::gpar(fill = patient_id_cols),
                #labels_gp = grid::gpar(fontsize = 10)
                #))}
            #)
    } else {
      ComplexHeatmap::draw(ht, heatmap_legend_list = hp,
                           padding = unit(c(3, 3, 3, 3), "mm"))
            #annotation_legend_list = if(
                #is.null(patient.split)) NULL else{
                #list(ComplexHeatmap::Legend(
                #title = "patient",
                #title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
                #at = patient.split,
                ##labels = names(patient_id_cols),
                ##legend_gp = grid::gpar(fill = patient_id_cols),
                #labels_gp = grid::gpar(fontsize = 10)
                #))}
  }


}
