suppressMessages(library(shiny))
suppressMessages(library(ggplot2))



# Define server logic required to plot various variables against mpg
shinyServer(function(input, output, session){
  observeEvent(input$help, {
    updateTabItems(session, "sidername", "home")
  })
  observeEvent(input$contact, {
    updateTabItems(session, "sidername", "home")
  })
  phylotree.type <- reactive({
    return(input$phylotTreeType)
  })
  width1 <- reactive({
    return(input$width1)
  })
  width2 <- reactive({
    return(input$width2)
  })
  width3 <- reactive({
    return(input$width3)
  })
  width4 <- reactive({
    return(input$width4)
  })
  width5 <- reactive({
    return(input$width5)
  })
  width6 <- reactive({
    return(input$width6)
  })
  height6 <- reactive({
    return(input$height6)
  })
  width7 <- reactive({
    return(input$width7)
  })
  height7 <- reactive({
    return(input$height7)
  })
  widthsig1 <- reactive({
    return(input$widthsig1)
  })
  heightsig1 <- reactive({
    return(input$heightsig1)
  })
  widthsig2 <- reactive({
    return(input$widthsig2)
  })
  heightsig2 <- reactive({
    return(input$heightsig2)
  })
  width11 <- reactive({
    return(input$width11)
  })
  height11 <- reactive({
    return(input$height11)
  })
  widthccfDen <- reactive({
    return(input$widthccfden)
  })
  
  
  
  mafName <- reactive({
    name <- input$mafFile$name
    patientID <- strsplit(name,"\\.")[[1]][1]
    return(patientID)
  })
  inputSilent <- observe({
    if (!is.null(input$mafFile)) {
      mafFile <- input$mafFile$datapath
      .substrRight <- function(x, n){
        substr(x, nchar(x)-n+1, nchar(x))
      }
      ## read maf file
      if (.substrRight(mafFile, 3) == ".gz") {
        mafInput <- read.table(mafGz <- gzfile(mafFile, "r"), quote="",
                               header=TRUE, fill=TRUE,
                               sep='\t')
        close(mafGz)
      } else {
        mafInput <- read.table(mafFile, quote="",
                               header=TRUE, fill=TRUE,
                               sep='\t')
      }
    } else {
      mafFile <- system.file("extdata", "HCC_LDC.maf", package = "MesKit")
      ccfFile <- system.file("extdata", "HCC_LDC.ccf.tsv", package = "MesKit")
      mafInput <- read.table(mafFile, quote="",
                             header=TRUE, fill=TRUE,
                             sep='\t')
    }
    colMt <- unique(mafInput$Variant_Classification)
    updateSelectInput(session, "mutNonSilent", 
                      choices=colMt, 
                      selected = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                                   "Translation_Start_Site", "Nonsense_Mutation", 
                                   "Nonstop_Mutation", "In_Frame_Del",
                                   "In_Frame_Ins", "Missense_Mutation"))
    return(colMt)
  })
  
  inputData <- eventReactive(input$submit1, {
    if(input$submit1){
      
      if (!is.null(input$mutNonSilent)){
        ls.mutNonSilent <- strsplit(input$mutNonSilent, ",")
      } else {
        ls.mutNonSilent <- "Default"
        updateSelectInput(session, "mutNonSilent", 
                          selected = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", 
                                       "Translation_Start_Site", "Nonsense_Mutation", 
                                       "Nonstop_Mutation", "In_Frame_Del",
                                       "In_Frame_Ins", "Missense_Mutation"))
      }
      
      # if (!is.null(input$chrSilent)){
      #   ls.chrSilent <- strsplit(input$chrSilent, ",")
      # } else {
      #   ls.chrSilent <- NULL
      # }
      
      if(is.null(input$mafFile)){
        mafFile <- system.file("extdata", "HCC6046.maf", package = "MesKit")
        ccfFile <- system.file("extdata", "HCC6046.ccf.tsv", package = "MesKit")
        maf <- readMaf(mafFile = mafFile,ccfFile = ccfFile)
      } else {
        if(!is.null(input$mafFile) & !is.null(input$ccfFile)){
          maf <- readMaf(mafFile = input$mafFile$datapath,
                         ccfFile =  input$ccfFile$datapath)          
          # if(maf@patientID == "0"){
          #     name <- as.character(mafName()) 
          #     maf@patientID <- name
          # }
        }
        else{
          maf <-  readMaf(mafFile = input$mafFile$datapath)
          # if(maf@patientID == "0"){
          #     name <- as.character(mafName()) 
          #     maf@patientID <- name
          # }
        }
      }
      return(maf)
    }
  })
  
  varsLs <- reactiveValues()
  
  observeEvent(input$submit1,{
    withProgress(min = 0, max = 2, value = 0, {
      ## Rshiny: progress bar
      setProgress(message = 'Input data: Generating ', detail = paste("MAF Class", sep="")) 
      varsLs[['maf']] <- inputData()
      incProgress(amount=1)
      
      ## Rshiny: progress bar
      setProgress(message = 'Input data: Generating ', detail = paste("phyloTree from MAF ",sep="")) 
      colNames <- colnames(inputData()@data)
      standardCol <- c("Hugo_Symbol","Chromosome","Start_Position","End_Position",
                       "Variant_Classification", "Variant_Type", "Reference_Allele",
                       "Tumor_Seq_Allele2", "VAF", "Tumor_Sample_Barcode")
      is <- intersect(colNames,standardCol)
      if(length(is)== 10){
          varsLs[['phyloTree']] <-  phyloTree <- getPhyloTree(isolate(varsLs$maf),method = input$method)
      }
      incProgress(amount=1)
      
      setProgress(message = paste("Input data: MAF and phyloTree Generation Done!", sep=""), detail = "") 
      Sys.sleep(1)
      
    })
  })
  
  buttonValue <- reactiveValues(a = 0, b = 0, c = 0, d = 0)
  observeEvent(input$iecontrol01,{
    buttonValue$a <- buttonValue$a + 1
    if(buttonValue$a == 2){
      buttonValue$a <- 0
    }
    buttonValue$b <- 0
    buttonValue$c <- 0
    buttonValue$d <- 0
  })
  observeEvent(input$iecontrol02,{
    buttonValue$a <- 0
    buttonValue$b <- buttonValue$b + 1
    if(buttonValue$b == 2){
      buttonValue$b <- 0
    }
    buttonValue$c <- 0
    buttonValue$d <- 0
  })
  observeEvent(input$iecontrol03,{
    buttonValue$a <- 0
    buttonValue$b <- 0
    buttonValue$c <- buttonValue$c + 1
    if(buttonValue$c == 2){
      buttonValue$c <- 0
    }
    buttonValue$d <- 0
  })
  observeEvent(input$iecontrol04,{
    buttonValue$a <- 0
    buttonValue$b <- 0
    buttonValue$c <- 0
    buttonValue$d <- buttonValue$d + 1
    if(buttonValue$d == 2){
      buttonValue$d <- 0
    }
  })
  ## output Introduction of maf datatable
  output$ie1 <- renderUI({
    if(buttonValue$a == 1){
      tagList(
        box(
          width = NULL,
          div(
            h3(strong("The MAF files")),
            p("MAF files contain many fields of information about chromosome and gene mutations and their annotations. The following fields are highly recommended to be contained in the MAF files.",
              style = "font-size:16px; font-weight:500;line-height:30px;"),
            h4(strong("Mandatory fields:")),
            p("Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, VAF, Tumor_Sample_Barcode.",
              style = "font-size:16px; font-weight:500;line-height:30px;"),
            h4(strong("Example MAF file"))
          ),
          DT::dataTableOutput("ied1",width = "100%"),
          br()
        )
      )
    }
  })
  output$ied1 <- renderDataTable({
    if(input$iecontrol01){
      maftable <- read.table('dom/maf.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(maftable, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
    }
  })
  ## output Introduction of sampleinfo datatable
  # output$ie2 <- renderUI({
  #   if(buttonValue$b == 1){
  #     box(
  #       width = NULL,
  #       tagList(
  #         div(
  #           h3(strong("Information of samples")),
  #           p("Below is an example of the first four rows of sample_info.txt. It should contain the sampleID, patientID, lesion and sampling time. The input files are located under the '/inst/extdata/' folder.",
  #             style = "font-size:16px; font-weight:500;line-height:30px;"),
  #           style = "width : 800px"
  #         ),
  #         tags$li(strong("tumors sampling across multiple spatially-distinct regions"),
  #                 style = "width : 800px;font-size:16px; font-weight:500;"),
  #         br(),
  #         DT::dataTableOutput("ied2_1",width = "70%"),
  #         br(),
  #         tags$li(strong("tumors sampling across multiple time points"),
  #                 style = "width : 800px;font-size:16px; font-weight:500;"),
  #         br(),
  #         DT::dataTableOutput("ied2_2",width = "70%"),
  #         br(),
  #         div(tags$span("Note:",style = "font-size:16px; font-weight:700;"),
  #             tags$span("'-' represents sampling at the same time or sampling from the same site.",
  #                       style = "font-size:16px; font-weight:500;")),
  #         br()
  #       )
  #     )
  #   }
  # })
  # output$ied2_1 <- renderDataTable({
  #   if(input$iecontrol02){
  #     spd1 <- read.table('dom/sampleinfo1.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
  #     datatable(spd1, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
  #   }
  # })
  # output$ied2_2 <- renderDataTable({
  #   if(input$iecontrol02){
  #     spd1 <- read.table('dom/sampleinfo2.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
  #     datatable(spd1, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
  #   }
  # })
  ## output Introduction of CCF
  output$ie3 <- renderUI({
    if(buttonValue$c == 1){
      box(
        width = NULL,
        tagList(
          h3(strong("The CCF file")),
          p("CCF files contain cancer cell fraction of each mutation.",
            style = "font-size:16px; font-weight:500;line-height:30px;"),
          h4(strong("Mandatory fields:")),
          p("Chromosome,Start,Sample,CCF",
            style = "font-size:16px; font-weight:500;line-height:30px;"),
          h4(strong("Example CCF file:")),
          DT::dataTableOutput("ied3"),
          br()
        )
      )
    }
  })
  output$ied3 <- renderDataTable({
    if(input$iecontrol03){
      spd3 <- read.table('dom/ccf.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(spd3, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
    }
  })
  ## output Introduction of ccf.loci
  # output$ie4 <- renderUI({
  #   if(buttonValue$d == 1){
  #     box(
  #       width = NULL,
  #       tagList(
  #         h3(strong("Example ccf.loci file")),
  #         DT::dataTableOutput("ied4"),
  #         br()
  #       )
  #     )
  #   }
  # })
  # output$ied4 <- renderDataTable({
  #   if(input$iecontrol04){
  #     spd4 <- read.table('dom/ccf.loci.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
  #     datatable(spd4, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
  #   }
  # })
  output$datapreview <- renderUI({
      if(input$submit1){
          colNames <- colnames(inputData()@data)
          standardCol <- c("Hugo_Symbol","Chromosome","Start_Position","End_Position",
                           "Variant_Classification", "Variant_Type", "Reference_Allele",
                           "Tumor_Seq_Allele2", "VAF", "Tumor_Sample_Barcode")
          is <- intersect(colNames,standardCol)
          if(length(is) == 10){
              DT::dataTableOutput('maftable', width = '100%')
          }
          else{
              uiOutput("warningMessage00")
          }
      }
  })
  output$maftable <- DT::renderDataTable({
    if(input$submit1){
        datatable(inputData()@data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)  
    }
  })
  stopButtonValue2 <- reactiveValues(a = 0)
  observeEvent(input$stop2,{
    stopButtonValue2$a <- 1
  })
  observeEvent(input$submit3,{
    stopButtonValue2$a <- 0
  })
  stopButtonValue2 <- reactiveValues(a = 0)
  observeEvent(input$stop2,{
    stopButtonValue2$a <- 1
  })
  observeEvent(input$submit2,{
    stopButtonValue2$a <- 0
  })
  ms <- eventReactive(input$submit2, {
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'MATH Score: Calculation in progress',
                 detail = 'This may take a while...')
    
    if(input$submit2 & stopButtonValue2$a != 1){
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(maf)), "")
      )
      mathScore(maf,min.vaf = input$minvaf)
    }
  })
  output$mathScore <- DT::renderDataTable({
    ms()
  })
  # ms2 <- eventReactive(input$submit0, {
  #   progress <- Progress$new(session, min=1, max=15)
  #   on.exit(progress$close())
  #   progress$set(message = 'TMB: Calculation in progress',
  #                detail = 'This may take a while...')
  #   
  #   if(input$submit0){
  #     maf <- isolate(varsLs$maf)
  #     validate(
  #         need(!(is.null(maf)), "")
  #     )
  #     getTMB(maf,tsb = NULL,
  #                       minvaf = input$minvaf, 
  #                       maxvaf = input$maxvaf)$sampleLevel
  #   }
  # })
  # output$mathScoreTMB <- DT::renderDataTable({
  #   ms2()
  # })
  
  output$msdb <- renderUI({
    if(!is.null(ms())){
      fluidRow(
        column(
          width = 9
        ),
        column(
          width = 3,
          downloadBttn('DownloadMathScore', 'Download')
        )
      )
    }
  })
  
  # output$msdbtmb <- renderUI({
  #   if(!is.null(ms2())){
  #     fluidRow(
  #       column(
  #         width = 9
  #       ),
  #       column(
  #         width = 3,
  #         downloadBttn('DownloadTMB', 'Download')
  #       )
  #     )
  #   }
  # })
  
  stopButtonValue3 <- reactiveValues(a = 0)
  observeEvent(input$stop3,{
    stopButtonValue3$a <- 1
  })
  observeEvent(input$submit3,{
    stopButtonValue3$a <- 0
  })
  sg <- eventReactive(input$useseg,{
      if(!is.null(input$segFile)){
          seg <- readSegment(segCN.file = input$segFile$datapath)
      }
      else{
          seg <- NULL
      }
      return(seg)
  })
  vc <- eventReactive(input$submit3, {
    if(input$submit3 & stopButtonValue3$a != 1){
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(maf)), "")
      )
      pmax <- length(unique(maf@data$Patient_ID))
      plot.list <- NULL
      
      withProgress(min = 0, max = pmax+1, value = 0, {
        setProgress(message = 'vafCluster: Calculation in progress',
                    detail = 'This may take a while...')
          
        plot.list <- vafCluster(maf,
                                plotOption = input$plotOption, 
                                showMATH = input$showMATH,
                                segCN.file = input$segFile$datapath)
        return(plot.list)
     })
    }
  })
  

  
  getpatient.vafcluster <- eventReactive(input$vc.pl, {
      return(input$vc.pl)
  })
  
  output$vafcluster.patientlist <- renderUI({
      if(!is.null(vc())){
          plot.list <- vc()
          names <- names(plot.list)
          selectInput("vc.pl", "Patient",
                      choices = names, width = 600)
      }
  })
  
  output$vaf <- renderPlot({
      if(!is.null(vc())){
          plot.list <- vc()
          plot.list[[getpatient.vafcluster()]]
      }
  },width = width1,
  height = 560,
  res = 100)
  
  output$vcdb <- renderUI({
    if(!is.null(vc())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons(inputId = 'DownloadVafPlotCheck', 
                       label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),
                       choiceValues = c("png", "pdf"), 
                       inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadVafPlot', 'Download')
        )
      )
    }
  })
  
  ccfauc <- eventReactive(input$submit_ccfauc,{
      withProgress(min = 0, max = 2, value = 0, {
          setProgress(message = 'ccfAUC: Calculation in progress',
                      detail = 'This may take a while...')
          maf <- varsLs$maf
          validate(
              need(!(is.null(maf)), "")
          )
          cc <- ccfAUC(maf, min.ccf = input$minccf_ccfauc, withinType = input$withintype_ccfauc)
          incProgress(amount = 1)
          setProgress(message = 'ccfAUC: Calculation done!')
      })
      return(cc)
  })
  
  output$ccfauc.patientlist <- renderUI({
      if(!is.null(ccfauc())){
          plot.list <- ccfauc()$CCF.density.plot
          names <- names(plot.list)
          tagList(
              selectInput("auc.pl", "Patient",
                          choices = names, width = 600) 
          )
      }
  })
  
  getpatient.ccfauc <- eventReactive(input$auc.pl,{
      return(input$auc.pl)
  })
  
  output$ccfauc_plot <- renderPlot({
      if(!is.null(ccfauc())){
          return(ccfauc()$CCF.density.plot[[getpatient.ccfauc()]])
      }
  },  
  width = 560,
  height = 560,
  res = 100)
  
  output$ccfauc_db_ui <- renderUI({
      if(!is.null(vc())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_ccfauc_plot_check', 
                               label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                               choiceNames = list(
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                               ),
                               choiceValues = c("png", "pdf"), 
                               inline = T)
              ),
              column(
                  width = 3,
                  downloadBttn('Download_ccfauc_plot', 'Download')
              )
          )
      }
  })
  
  output$ccfauc_table <- DT::renderDataTable({
      if(!is.null(ccfauc())){
          t <- ccfauc()$AUC.value
          rows <- which(t$Patient_ID == getpatient.ccfauc())
          dt <- datatable(t[rows,][,1:3],
                           options = list(searching = TRUE,
                                          pageLength = 10, 
                                          scrollX = TRUE,
                                          dom = "t",
                                          fixedHeader = TRUE),
                           rownames = F) 
          return(dt)
      }
  })
  
  
  output$ccfauc_table_ui <- renderUI({
      if(!is.null(ccfauc())){
          tagList(
              h4(strong('AUC value')),
              br(),
              DT::dataTableOutput('ccfauc_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_ccfauc_table', 'Download')
                  )
              )
          )
      }
  })
  
  ## calfst sever
  calfst <- eventReactive(input$submit_calfst,{
      withProgress(min = 0, max = 2, value = 0, {
          setProgress(message = 'ccfAUC: Calculation in progress',
                      detail = 'This may take a while...')
          maf <- varsLs$maf
          validate(
              need(!(is.null(maf)), "")
          )
          cc <- calFst(maf, min.vaf = input$minvaf_calfst, withinType = input$withintype_calfst)
          incProgress(amount = 1)
          setProgress(message = 'ccfAUC: Calculation done!')
      })
      return(cc)
  })
  
  output$calfst.patientlist <- renderUI({
      if(!is.null(calfst())){
          plot.list <- calfst()$Fst.plot
          names <- names(plot.list)
          tagList(
              selectInput("calfst.pl", "Patient",
                          choices = names, width = 600) 
          )
      }
  })
  
  getpatient.calfst <- eventReactive(input$calfst.pl,{
      return(input$calfst.pl)
  })
  
  output$calfst_plot <- renderPlot({
      if(!is.null(calfst())){
          return(calfst()$Fst.plot[[getpatient.calfst()]])
      }
  },  
  width = 560,
  height = 560,
  res = 100)
  
  output$calfst_db_ui <- renderUI({
      if(!is.null(calfst())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_calfst_plot_check', 
                               label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                               choiceNames = list(
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                               ),
                               choiceValues = c("png", "pdf"), 
                               inline = T)
              ),
              column(
                  width = 3,
                  downloadBttn('Download_calfst_plot', 'Download')
              )
          )
      }
  })
  
  output$calfst_avg_table <- DT::renderDataTable({
      if(!is.null(calfst())){
          t <- calfst()$Fst.avg
          rows <- which(t$Patient_ID == getpatient.calfst())
          dt <- datatable(t[rows,],
                          options = list(searching = TRUE,
                                         pageLength = 10, 
                                         scrollX = TRUE,
                                         dom = "t",
                                         fixedHeader = TRUE),
                          rownames = F) 
          return(dt)
      }
  })
  
  
  output$calfst_avg_table_ui <- renderUI({
      if(!is.null(calfst())){
          tagList(
              h4(strong('Fst average value')),
              br(),
              DT::dataTableOutput('calfst_avg_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_calfst_avg_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$calfst_pair_table <- DT::renderDataTable({
      if(!is.null(calfst())){
          m <- calfst()$Fst.pair
          rownames(m) <- colnames(m)
          m <- as.data.frame(m)
          colnames(m) <- rownames(m)
          dt <- datatable(m,options = list(searching = TRUE,
                                         pageLength = 10, 
                                         scrollX = TRUE,
                                         dom = "t",
                                         fixedHeader = TRUE),
                          rownames = T) 
          return(dt)
      }
  })
  
  
  output$calfst_pair_table_ui <- renderUI({
      if(!is.null(calfst())){
          tagList(
              h4(strong('Fst pair')),
              br(),
              DT::dataTableOutput('calfst_pair_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_calfst_pair_table', 'Download')
                  )
              )
          )
      }
  })
  
  
  ## neidist sever
  calneidist <- eventReactive(input$submit_calneidist,{
      withProgress(min = 0, max = 2, value = 0, {
          setProgress(message = 'calNeiDist: Calculation in progress',
                      detail = 'This may take a while...')
          maf <- varsLs$maf
          validate(
              need(!(is.null(maf)), "")
          )
          cc <- calNeiDist(maf, min.vaf = input$minvaf_calneidist, withinType = input$withintype_calneidist)
          incProgress(amount = 1)
          setProgress(message = 'ccfAUC: Calculation done!')
      })
      return(cc)
  })
  
  output$calneidist.patientlist <- renderUI({
      if(!is.null(calneidist())){
          plot.list <- calneidist()$Nei.plot
          names <- names(plot.list)
          tagList(
              selectInput("calneidist.pl", "Patient",
                          choices = names, width = 600) 
          )
      }
  })
  
  getpatient.calneidist <- eventReactive(input$calneidist.pl,{
      return(input$calneidist.pl)
  })
  
  output$calneidist_plot <- renderPlot({
      if(!is.null(calneidist())){
          return(calneidist()$Nei.plot[[getpatient.calneidist()]])
      }
  },  
  width = 560,
  height = 560,
  res = 100)
  
  output$calneidist_db_ui <- renderUI({
      if(!is.null(calneidist())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_calneidist_plot_check', 
                               label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                               choiceNames = list(
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                               ),
                               choiceValues = c("png", "pdf"), 
                               inline = T)
              ),
              column(
                  width = 3,
                  downloadBttn('Download_calneidist_plot', 'Download')
              )
          )
      }
  })
  
  output$calneidist_avg_table <- DT::renderDataTable({
      if(!is.null(calneidist())){
          t <- calneidist()$Nei.dist.avg
          rows <- which(t$Patient_ID == getpatient.calneidist())
          dt <- datatable(t[rows,],
                          options = list(searching = TRUE,
                                         pageLength = 10, 
                                         scrollX = TRUE,
                                         dom = "t",
                                         fixedHeader = TRUE),
                          rownames = F) 
          return(dt)
      }
  })
  
  
  output$calneidist_avg_table_ui <- renderUI({
      if(!is.null(calneidist())){
          tagList(
              h4(strong('Nei dist average value')),
              br(),
              DT::dataTableOutput('calneidist_avg_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_calneidist_avg_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$calneidist_pair_table <- DT::renderDataTable({
      if(!is.null(calneidist())){
          m <- calneidist()$Nei.dist
          rownames(m) <- colnames(m)
          m <- as.data.frame(m)
          colnames(m) <- rownames(m)
          dt <- datatable(m,options = list(searching = TRUE,
                                           pageLength = 10, 
                                           scrollX = TRUE,
                                           dom = "t",
                                           fixedHeader = TRUE),
                          rownames = T) 
          return(dt)
      }
  })
  
  
  output$calneidist_pair_table_ui <- renderUI({
      if(!is.null(calneidist())){
          tagList(
              h4(strong('Nei dist pair')),
              br(),
              DT::dataTableOutput('calneidist_pair_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_calneidist_pair_table', 'Download')
                  )
              )
          )
      }
  })
   
  ## comparejsi sever
  comparejsi <- eventReactive(input$submit_comparejsi,{
      withProgress(min = 0, max = 2, value = 0, {
          setProgress(message = 'compareJSI : Calculation in progress',
                      detail = 'This may take a while...')
          maf <- varsLs$maf
          validate(
              need(!(is.null(maf)), "")
          )
          cc <- compareJSI(maf, min.vaf = input$minvaf_comparejsi, pairByType = input$pairbytype_comparejsi)
          incProgress(amount = 1)
          setProgress(message = 'ccfAUC: Calculation done!')
      })
      return(cc)
  })
  
  output$comparejsi.patientlist <- renderUI({
      if(!is.null(comparejsi())){
          plot.list <- comparejsi()$JSI.plot
          names <- names(plot.list)
          tagList(
              selectInput("comparejsi.pl", "Patient",
                          choices = names, width = 600) 
          )
      }
  })
  
  getpatient.comparejsi <- eventReactive(input$comparejsi.pl,{
      return(input$comparejsi.pl)
  })
  
  output$comparejsi_plot <- renderPlot({
      if(!is.null(comparejsi())){
          return(comparejsi()$JSI.plot[[getpatient.comparejsi()]])
      }
  },  
  width = 560,
  height = 560,
  res = 100)
  
  output$comparejsi_db_ui <- renderUI({
      if(!is.null(comparejsi())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_comparejsi_plot_check', 
                               label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                               choiceNames = list(
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                               ),
                               choiceValues = c("png", "pdf"), 
                               inline = T)
              ),
              column(
                  width = 3,
                  downloadBttn('Download_comparejsi_plot', 'Download')
              )
          )
      }
  })
  
  output$comparejsi_avg_table <- DT::renderDataTable({
      if(!is.null(comparejsi())){
          t <- comparejsi()$JSI.multi
          rows <- which(t$Patient_ID == getpatient.comparejsi())
          dt <- datatable(t[rows,],
                          options = list(searching = TRUE,
                                         pageLength = 10, 
                                         scrollX = TRUE,
                                         dom = "t",
                                         fixedHeader = TRUE),
                          rownames = F) 
          return(dt)
      }
  })
  
  
  output$comparejsi_avg_table_ui <- renderUI({
      if(!is.null(comparejsi())){
          tagList(
              h4(strong('JSI average value')),
              br(),
              DT::dataTableOutput('comparejsi_avg_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_comparejsi_avg_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$comparejsi_pair_table <- DT::renderDataTable({
      if(!is.null(comparejsi())){
          m <- comparejsi()$JSI.pair
          m <- as.data.frame(m)
          dt <- datatable(m,options = list(searching = TRUE,
                                           pageLength = 10, 
                                           scrollX = TRUE,
                                           dom = "t",
                                           fixedHeader = TRUE),
                          rownames = T) 
          return(dt)
      }
  })
  
  
  output$comparejsi_pair_table_ui <- renderUI({
      if(!is.null(comparejsi())){
          tagList(
              h4(strong('JSI pair')),
              br(),
              DT::dataTableOutput('comparejsi_pair_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_comparejsi_pair_table', 'Download')
                  )
              )
          )
      }
  })
  
  stopButtonValue4 <- reactiveValues(a = 0)
  observeEvent(input$stop4,{
    stopButtonValue4$a <- 1
  })
  observeEvent(input$submit4,{
    stopButtonValue4$a <- 0
  })
  msp <- eventReactive(input$submit4, {
    if(input$submit4 & stopButtonValue4$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Mutsharedprivateplot: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(maf)), "")
      )
      return(mutPrivateShared(maf, show.num = input$showNum1))
    }
  })
  # 
  # output$mutSharedPrivatePlot <- renderPlot({
  #   if (input$submit4 & input$plotChoiceSpp == "sharedPrivatePlot"){
  #     msp()
  #   } else if (input$submit5 & input$plotChoiceSpp == "stackPlot") {
  #     stk()
  #   }
  # },
  # width = width2,
  # height = 560,
  # res = 100
  # )
  output$mutSharedPrivatePlot <- renderPlot({
    msp()
  },
  width = width2,
  height = 560,
  res = 100
  )
  output$mspdb <- renderUI({
    if(!is.null(msp())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons('DownloadSharedPlotCheck',
                       label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),
                       choiceValues = c("png", "pdf"), 
                       inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadSharedPlot', 'Download')
        )
      )
    }
    # else if(!is.null(stk())){
    #     fluidRow(
    #       column(
    #         width = 7
    #       ),
    #       column(
    #         width = 2,
    #         radioButtons('DownloadStackPlottCheck',
    #                      label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
    #                      choiceNames = list(
    #                        tags$span(style = "font-size:14.5px; font-weight:400; ", "png"),
    #                        tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
    #                      ),
    #                      choiceValues = c("png", "pdf"),
    #                      inline = T)
    #       ),
    #       column(
    #         width = 3,
    #         downloadBttn('DownloadStackPlot', 'Download')
    #       )
    #     )
    #   }
  })
  stopButtonValue5 <- reactiveValues(a = 0)
  observeEvent(input$stop5,{
    stopButtonValue5$a <- 1
  })
  observeEvent(input$submit5,{
    stopButtonValue5$a <- 0
  })
  stk <- eventReactive(input$submit5, {
    if(input$submit5 & stopButtonValue5$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'mutOncoTSG: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      if(is.null(input$oncogeneListFile$datapath)){
        oncogeneListFile <- system.file("extdata/", "oncogene.list.txt", package = "MesKit")
      }
      else{
        oncogeneListFile <- input$oncogeneListFile$datapath
      }
      if(is.null(input$tsgListFile$datapath)){
        tsgListFile <- system.file("extdata/", "TSG.list.txt", package = "MesKit")
      }
      else{
        tsgListFile <- input$tsgListFile$datapath
      }
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(maf)), "")
      )
      mutOncoTSG(maf, oncogeneListFile = oncogeneListFile,
                   tsgListFile = tsgListFile, 
                   show.percentage = input$showPercentage)
    }
  })
  output$mutoncotsg <- renderPlot({
    stk()
  },
  width = width3,
  height = 560,
  res = 100
  )
  output$stkdb <- renderUI({
    if(!is.null(stk())){
      fluidRow(
        column(
          width = 6
        ),
        column(
          width = 2,
          radioButtons('DownloadStackPlotCheck',
                       label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"),
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),
                       choiceValues = c("png", "pdf"),
                       inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadStackPlot', 'Download')
        )
      )
    }
  })
  stopButtonValue6 <- reactiveValues(a = 0)
  observeEvent(input$stop6,{
    stopButtonValue6$a <- 1
  })
  observeEvent(input$submit6,{
    stopButtonValue6$a <- 0
  })
  ji <- eventReactive(input$submit6, {
    if(stopButtonValue6$a != 1&input$submit6){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Paired-samples similarity: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(maf)), "")
      )
      return(MesKit::JaccardIndex(maf,type = input$JItype))
    }
    # progress <- Progress$new(session, min=1, max=15)
    # on.exit(progress$close())
    # progress$set(message = 'Calculation in progress',
    #              detail = 'This may take a while...')
    # 
    # for (i in 1:15) {
    #   progress$set(value = i)
    #   Sys.sleep(0.01)
    # }
    # maf <- isolate(varsLs$maf)
    # return(MesKit::JaccardIndex(maf,type = input$JItype))
  })
  output$JaccardIndex <- renderPlot({
    ji()
  },
  width = width4,
  height = width4,
  res = 100
  )
  output$jidb <- renderUI({
    if(!is.null(ji())){
      fluidRow(
        column(
          width = 5
        ),
        column(
          width = 2,
          radioButtons('DownloadJaccardIndexCheck',label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),
                       choiceValues = c("png", "pdf"), 
                       inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadJaccardIndex', 'Download')
        )
      )
    }
  })
  stopButtonValue7 <- reactiveValues(a = 0)
  observeEvent(input$stop7,{
    stopButtonValue7$a <- 1
  })
  observeEvent(input$submit7,{
    stopButtonValue7$a <- 0
  })
  clp <- eventReactive(input$submit7, {
    if(input$submit7 & stopButtonValue7$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Subclonal plot: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      if(!is.null(input$mafFile) & !is.null(input$sampleInfo)){
        validate(
          need(!(is.null(input$ccf.cluster$datapath)), "click the button 'use ccf',Upload ccf.cluster in Session 'Input Data' ")
        )
        validate(
          need(!(is.null(input$ccf.loci$datapath)), "Upload ccf.loci in Session 'Input Data'")
        )
        maf <- isolate(varsLs$maf)
        validate(
            need(!(is.null(maf)), "")
        )
        p <- MesKit::tumorClonesPlot(maf)
        return(p)
      }
      else{
        maf <- isolate(varsLs$maf)
        p <- MesKit::tumorClonesPlot(maf)
        return(p)
      }
    }
  })
  output$cloneplot <- renderPlot({
    clp()
  },
  width = width5,
  height = 560,
  res = 100
  )
  output$clpdb <- renderUI({
    if(!is.null(clp())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons('DownloadClonePlotCheck',label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),
                       choiceValues = c("png", "pdf"), 
                       inline = T)
        ),
        column(
          width = 3,
          div(downloadBttn('DownloadClonePlot', 'Download'))
        )
      )
    }
  })
  
  ccfden <- eventReactive(input$submitccfden, {
    if(input$submitccfden){
      progress <- Progress$new(session, min=0, max=1)
      on.exit(progress$close())
      progress$set(message = 'CCF Density: Calculation in progress',
                   detail = 'This may take a while...')
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(maf)), "")
      )
      cd <- ccfDensity(maf, show.density = input$showdensity)
      progress$set(value = 1)
      return(cd)
    }
  })
  output$ccfdenplot <- renderPlot({
    ccfden()
  },
  width = widthccfDen,
  height = 560,
  res = 100
  )
  output$ccfdendb <- renderUI({
    if(!is.null(ccfden())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons('DownloadCCFDensityCheck',label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),
                       choiceValues = c("png", "pdf"), 
                       inline = T)
        ),
        column(
          width = 3,
          div(downloadBttn('DownloadCCFDensity', 'Download'))
        )
      )
    }
  })
  
  stopButtonValue8 <- reactiveValues(a = 0)
  observeEvent(input$stop8,{
    stopButtonValue8$a <- 1
  })
  observeEvent(input$submit8,{
    stopButtonValue8$a <- 0
  })
  GO <- eventReactive(input$submit8, {
    if(input$submit8 & stopButtonValue8$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'GO analysis: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      phyloTree <- isolate(varsLs$phyloTree)
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(phyloTree)), "")
      )
      if (input$driverGenesMapping1) {
          if(is.null(input$driverGenesFile1$datapath)){
              driverGenesFile <- system.file("extdata", "putative_driver_genes.txt", package = "MesKit")
          } else{
              driverGenesFile <- input$driverGenesFile1$datapath
          }
      }
      else{
          driverGenesFile <- NULL
      }
      treeGO(phyloTree, 
                        GO.type = input$GOtype, 
                        plotType = input$plotType, 
                        pAdjustMethod=input$pAdjustMethod, 
                        qval = as.numeric(input$qval1), 
                        pval = as.numeric(input$pval1), 
                        showCategory = input$showCategory,
                        driverGenesFile = driverGenesFile)
    }
  })
  # Datatable under GO plot
  output$gotui <- renderUI({
    if(!is.null(GO())){
      tagList(
        h4(strong('Result list')),
        br(),
        DT::dataTableOutput('gotable'),
        fluidRow(
          column(
            width = 9
          ),
          column(
            width = 3,
            br(),
            downloadBttn('DownloadGOTable', 'Download')
          )
        )
      )
    }
  })
  output$gotable <- renderDataTable({
    data <- GO()[[1]][[which(names(GO()[[1]]) == input$gl)]]
    # targets <- which(colnames(data) %in% c("geneID"))-1
    datatable(data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
  })
  output$chooselist1 <- renderUI({
    if(!is.null(GO())){
      if("All" %in% names(GO()[[2]])){
        names <- names(GO()[[2]])
        names <- names[-(which(names == "All"))]
        names <- append(names,"All",after = 0)
      }
      else{
        names <- names(GO()[[2]])
      }
      selectInput("gl","Branch",
                  choices = names ,selected = names[1],width = 600)
    }
  })
  output$GOplot <- renderPlot({
    return(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
  },
  width = width6,
  height = height6
  )
  output$GOdb <- renderUI({
    if(!is.null(GO())){
      fluidRow(
        column(
          width = 6
        ),
        column(
          width = 2,
          radioButtons('DownloadGOPlotCheck',label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),choiceValues = c("png", "pdf"), inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadGOPlot', 'Download')
        )
      )
    }
  })
  stopButtonValue9 <- reactiveValues(a = 0)
  observeEvent(input$stop9,{
    stopButtonValue9$a <- 1
  })
  observeEvent(input$submit9,{
    stopButtonValue9$a <- 0
  })
  Path <- eventReactive(input$submit9, {
    if(input$submit9 != 0 & stopButtonValue9$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Pathway analysis: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      phyloTree <- isolate(varsLs$phyloTree)
      maf <- isolate(varsLs$maf)
      validate(
          need(!(is.null(phyloTree)), "")
      )
      if (input$driverGenesMapping2) {
          if(is.null(input$driverGenesFile2$datapath)){
              driverGenesFile <- system.file("extdata", "putative_driver_genes.txt", package = "MesKit")
          } else{
              driverGenesFile <- input$driverGenesFile2$datapath
          }
      }
      else{
          driverGenesFile <- NULL
      }
      list <- treePathway(phyloTree, driverGenesFile = driverGenesFile,
                                     pathway.type=input$pathwaytype, 
                                     plotType = input$pathplotType, 
                                     pAdjustMethod=input$pathpAdjustMethod, 
                                     qval = as.numeric(input$qval2), 
                                     pval = as.numeric(input$pval2),
                                     showCategory = input$pathshowCategory
                                     
      )
      return(list)
    }
  })
  output$Pathdb <- renderUI({
    if(!is.null(Path())){
      fluidRow(
        column(
          width = 6
        ),
        column(
          width = 2,
          radioButtons('DownloadPathPlotCheck',label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),choiceValues = c("png", "pdf"), inline = T)
        ),
        column(
          width = 3,
          br(),
          downloadBttn('DownloadPathPlot', 'Download')
        )
      )
      
    }
  })
  output$chooselist2 <- renderUI({
    if(!is.null(Path())){
      if("All" %in% names(Path()[[2]])){
        names <- names(Path()[[2]])
        names <- names[-(which(names == "All"))]
        names <- append(names,"All",after = 0)
      }
      else{
        names <- names(Path()[[2]])
      }
      selectInput("pl","Branch",
                  choices = names ,selected = names[1],width = 600)
    }
  })
  # Datatable under Pathway plot
  output$patht <- renderUI({
    if(!is.null(Path())){
      tagList(
        h4(strong('Result list')),
        br(),
        DT::dataTableOutput('pathtable'),
        fluidRow(
          column(
            width = 9
          ),
          column(
            width = 3,
            downloadBttn('DownloadPathTable', 'Download')
          )
        )
      )
    }
  })
  output$pathtable <- renderDataTable({
    data <- Path()[[1]][[which(names(Path()[[1]]) == input$pl)]]
    return(DT::datatable(data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)) 
    
  })
  output$Pathwayplot <- renderPlot({
    return(Path()[[2]][[which(names(Path()[[2]]) == input$pl)]]) 
  },
  width = width7,
  height = height7,
  res = 100
  )
  
  tms <- reactiveValues()
  
  stopButtonValueSig <- reactiveValues(a = 0)
  observeEvent(input$stopSig,{
    stopButtonValueSig$a <- 1
  })
  observeEvent(input$submitSig,{
    stopButtonValueSig$a <- 0
  })
  
  observeEvent(input$submitSig,{
      if(is.null(tree.mutSig[["value"]])){
          phyloTree <- isolate(varsLs$phyloTree)
          withProgress(min = 0, max = length(names(phyloTree))+1, value = 0,{
              setProgress(message = 'treeMutSig: Calculation in progress',
                          detail = 'This may take a while...')
              
              tms[["value"]] <- treeMutSig(phyloTree, 
                                                   geneList = NULL, 
                                                   min.mut.count = input$mutThreshold,
                                                   signaturesRef=input$signaturesRef)
             Sys.sleep(1)
          })
      }
  })
  sigOFA <- eventReactive(input$submitSig, {
    if (input$oncogeneMapping){
      if(is.null(input$driverGenesFile$datapath)){
        driverGenesFile <- system.file("extdata", "putative_driver_genes.txt", package = "MesKit")
      } else{
        driverGenesFile <- input$driverGenesFile$datapath
      }
        driverGene <- as.character(read.table(driverGene)$V1)
    }
    else{
        driverGene <- NULL 
    }
      
   phyloTree <- isolate(varsLs$phyloTree) 
   
   validate(
       need(!(is.null(phyloTree)), "")
   )
   
  if(input$submitSig & stopButtonValueSig$a != 1){
         
         mutSig.summary <- NULL
         

         mutSig.summary <- isolate(tms$value)$mutSig.summary
         
         output$sigsummary.patientlist <- renderUI({
             if(!is.null(mutSig.summary)){
                 names <- names(mutSig.summary)
                 tagList(
                     selectInput("ss.pl", "Patient",
                                 choices = names,
                                 selected = names[1],
                                 width = 600) 
                 )
             }
         })
         
         if(length(names(mutSig.summary))  == 1){
             n <- names(mutSig.summary)[1]
         }
         else{
             n <- getpatient.sigsummary()
             
         }
         
         s <- mutSig.summary[[n]]
         print(s)
         output$sigOFAt <- DT::renderDataTable({
             return(datatable(s, 
                              options = list(searching = TRUE, 
                                             pageLength = 10, 
                                             lengthMenu = c(5, 10, 15, 18), 
                                             scrollX = TRUE), 
                              rownames = FALSE))
         })
         # print(input$ss.pl)
         # idxp <- which(names(mutSig.summary) == "HCC6046")
         return(s)
  }
 })
  
  getpatient.sigsummary <- eventReactive(input$ss.pl,{
      return(input$ss.pl)
  })
  
  # observeEvent(input$submitSig4,{
  #     if(is.null(mtb[["value"]])){
  #         
  #         if (input$oncogeneMapping4){
  #             if(is.null(input$driverGenesFile$datapath)){
  #                 driverGenesFile <- system.file("extdata", "putative_driver_genes.txt", package = "MesKit")
  #             } else{
  #                 driverGenesFile <- input$driverGenesFile4$datapath
  #             }
  #             driverGene <- as.character(read.table(driverGenesFile)$V1)
  #         }
  #         else{
  #             driverGene <- NULL 
  #         }
  #         
  #         phyloTree <- isolate(varsLs$phyloTree)
  #         withProgress(min = 0, max = length(names(phyloTree))+1, value = 0,{
  #             setProgress(message = 'treeMutSig: Calculation in progress',
  #                         detail = 'This may take a while...')
  #             
  #             x <- mutTrunkBranch(phyloTree, 
  #                                 geneList = driverGene, 
  #                                 use.shiny = TRUE)
  #             
  #             mtb$summary <- plyr::rbind.fill(x$mutTrunkBranch.res)
  #             mtb$plot <- x$mutTrunkBranch.plot
  #             
  #             Sys.sleep(1)
  #         })
  #     }
  # })
  
  sigBT <- eventReactive(input$submitSig4, {
      if (input$oncogeneMapping4){
          if(is.null(input$driverGenesFile$datapath)){
              driverGenesFile <- system.file("extdata", "putative_driver_genes.txt", package = "MesKit")
          } else{
              driverGenesFile <- input$driverGenesFile4$datapath
          }
          driverGene <- as.character(read.table(driverGenesFile)$V1)
      }
      else{
          driverGene <- NULL 
      }
      
      phyloTree <- isolate(varsLs$phyloTree)
      withProgress(min = 0, max = length(names(phyloTree))+1, value = 0,{
          setProgress(message = 'mutTrunkBranch: Calculation in progress',
                      detail = 'This may take a while...')
          
          
          mtb <- mutTrunkBranch(phyloTree, 
                               geneList = driverGene)
          
          Sys.sleep(1)
      })
      
        return(plyr::rbind.fill(mtb$mutTrunkBranch.res))
  })
 
  output$sigBTt <- DT::renderDataTable({
    return(datatable(sigBT(), 
                     options = list(searching = TRUE, 
                                    pageLength = 10, 
                                    lengthMenu = c(5, 10, 15, 18), 
                                    scrollX = TRUE), 
                     rownames = FALSE))
  })
  
  
  stopButtonValueSig1 <- reactiveValues(a = 0)
  observeEvent(input$stopSig1,{
    stopButtonValueSig1$a <- 1
  })
  observeEvent(input$submitSig1,{
    stopButtonValueSig1$a <- 0
  })
  sigOFA1 <- eventReactive(input$submitSig1, {
      
      phyloTree <- isolate(varsLs$phyloTree)
      withProgress(min = 0, max = length(names(phyloTree))+1, value = 0,{
          setProgress(message = 'treeMutSig: Calculation in progress',
                      detail = 'This may take a while...')
          
          tms <- treeMutSig(phyloTree, 
                            min.mut.count = input$mutThreshold, 
                            signaturesRef=input$signaturesRef1)
          
          Sys.sleep(1)
      })
      
      df.signature <- tms$mutSig.summary
      plot.list <- plotMutSigProfiler(tms$mutSig.spectrum)
      return(list(plot.list, df.signature))
  })
  
  getpatient.tms <- eventReactive(input$tms.pl,{
      return(input$tms.pl)
  })
  
  output$treemutsig.patientlist <- renderUI({
      if(!is.null(sigOFA1())){
          plot.list <- sigOFA1()[[1]]
          names <- names(plot.list)
          selectInput("tms.pl", "Patient_branches",
                      choices = names, width = 600)
      }
  })
  
  output$sigOFATableUI1 <- renderUI({
    if(!is.null(sigOFA1())){
      tagList(
        h4(strong('Signature summary')),
        br(),
        DT::dataTableOutput('sigOFATable1'),
        br(),
        fluidRow(
          column(
            width = 9
          ),
          column(
            width = 3,
            downloadBttn('DownloadSigOFATable1', 'Download')
          )
        )
      )
    }
  })
  output$sigOFAPlot1 <- renderPlot({
    return(sigOFA1()[[1]][[getpatient.tms()]])
  },
  width = widthsig1,
  height = heightsig1,
  res = 100
  )
  output$sigOFATable1 <- DT::renderDataTable({
    data <- sigOFA1()[[2]][[getpatient.tms()]][,c(1:2)]
    datatable(data, options = list(searching = TRUE, pageLength = 10, 
                                   lengthMenu = c(5, 10, 15, 18), 
                                   scrollX = TRUE, dom = "t",
                                   fixedHeader = TRUE),rownames = F)
  })
  stopButtonValueSig2 <- reactiveValues(a = 0)
  observeEvent(input$stopSig2,{
    stopButtonValueSig2$a <- 1
  })
  observeEvent(input$submitSig2,{
    stopButtonValueSig2$a <- 0
  })
  sigOFA2 <- eventReactive(input$submitSig2, {
      phyloTree <- isolate(varsLs$phyloTree)
      withProgress(min = 0, max = length(names(phyloTree))+1, value = 0,{
          setProgress(message = 'treeMutSig: Calculation in progress',
                      detail = 'This may take a while...')
          
          
          mtb <- mutTrunkBranch(phyloTree, 
                                geneList = NULL)
          
          Sys.sleep(1)
      })
      return(mtb$mutTrunkBranch.plot)
  })
  getpatient.mtb <- eventReactive(input$mtb.pl,{
      return(input$mtb.pl)
  })
  
  output$mtb.patientlist <- renderUI({
      if(!is.null(sigOFA2())){
          plot.list <- sigOFA2()
          names <- names(plot.list)
          selectInput("mtb.pl", "Patient_branches",
                      choices = names, width = 600)
      }
  })
  # output$sigOFATableUI2 <- renderUI({
  #   if(!is.null(sigOFA2()[[2]]))
  #   tagList(
  #     h4(strong('Signature summary')),
  #     br(),
  #     DT::dataTableOutput('sigOFATable2'),
  #     fluidRow(
  #       column(
  #         width = 9
  #       ),
  #       column(
  #         width = 3,
  #         downloadBttn('DownloadSigOFATable2', 'Download')
  #       )
  #     )
  #   )
  # })
  # output$sigOFATable2 <- renderDataTable({
  #   if(!is.null(sigOFA2()[[2]])){
  #     return(datatable(sigOFA2()[[2]], options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T)))
  #   }
  # })
  
  output$sigOFAPlot2 <- renderPlot({
      if(!is.null(sigOFA2())){
          return(sigOFA2()[[getpatient.mtb()]]) 
      }
  },
  width = widthsig2,
  height = heightsig2,
  res = 100
  )
  output$sigpdb <- renderUI({
    if(!is.null(sigOFA())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2
        ),
        column(
          width = 3,
          downloadBttn('DownloadSignatureSummary', 'Download')
        )
      )
    }
  })
  output$sigbtdb <- renderUI({
    if(!is.null(sigBT())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2
        ),
        column(
          width = 3,
          downloadBttn('DownloadSignatureBT', 'Download')
        )
      )
    }
  })
  output$sigpdb1 <- renderUI({
    if(!is.null(sigOFA1())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons('DownloadSignaturePlotCheck1', label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),choiceValues = c("png", "pdf"), inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadSignaturePlot1', 'Download')
        )
      )
    }
  })
  output$sigpdb2 <- renderUI({
    if(!is.null(sigOFA2())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons('DownloadSignaturePlotCheck2', label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),choiceValues = c("png", "pdf"), inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadSignaturePlot2', 'Download')
        )
      )
    }
  })
  stopButtonValue10<- reactiveValues(a = 0)
  observeEvent(input$stop10,{
    stopButtonValue10$a <- 1
  })
  observeEvent(input$submit10,{
    stopButtonValue10$a <- 0
  })
  pht <- eventReactive(input$submit10, {
      
      phyloTree <- isolate(varsLs$phyloTree)
      withProgress(min = 0, max = length(names(phyloTree))+1, value = 0,{
          setProgress(message = 'Phylogenetic tree: Calculation in progress',
                      detail = 'This may take a while...')
          
          
          plot.list <- plotPhyloTree(phyloTree)
          
          Sys.sleep(1)
      })
      
      
      output$tree.patientlist <- renderUI({
          if(!is.null(plot.list)){
              names <- names(plot.list)
              selectInput("t.pl", "Patient",
                          choices = names, width = 600)
          }
      })
      
      if(length(names(plot.list)) == 1){
          n <- names(plot.list)[1]
      }else{
          n <- getpatient.tree()
      }
      
      return(plot.list[[n]])
  })
  
  getpatient.tree <- eventReactive(input$t.pl,{
      return(input$t.pl)
  })
  
  output$phylotree <- renderPlot({
    return(pht()) 
  },
  res = 100
  )
  output$phtdb <- renderUI({
    if(!is.null(pht())){
      br()
      br()
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons('DownloadPhyloTreeCheck', label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),
                       choiceValues = c("png", "pdf"), inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadPhyloTree', 'Download')
        )
      )
    }
  })
  
  ## Download control  
  output$DownloadMathScore <- downloadHandler(
    filename = 'Rtable.csv',
    content = function(file){
      data <- ms()
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadTMB <- downloadHandler(
    filename = "Rtable.csv",
    content = function(file){
      data <- ms2()
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadVafPlot <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadVafPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadVafPlotCheck == "png"){
        png(file,width = input$width1 , height = 560,res = 100)
      }
      else if (input$DownloadVafPlotCheck == "pdf"){
        pdf(file,width = input$width1/100 , height = 6)
      }
      print(vc())
      dev.off()
    },
    contentType = paste('image/',input$DownloadVafPlotCheck,sep="")
  )
  
  output$DownloadStackPlot <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadStackPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadStackPlotCheck == "png"){
        png(file,width = input$width3 , height = 560,res = 100)
      }
      else if (input$DownloadStackPlotCheck == "pdf"){
        pdf(file,width = input$width3/100 , height = 6)
      }
      print(stk())
      dev.off()
    },
    contentType = paste('image/',input$DownloadStackPlotCheck,sep="")
  )
  output$DownloadJaccardIndex <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadJaccardIndexCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadJaccardIndexCheck == "png"){
        png(file,width = input$width4 , height = input$width4,res = 100)
      }
      else if (input$DownloadJaccardIndexCheck == "pdf"){
        pdf(file,width = input$width4/100 , height = input$width4/100)
      }
     maf <- isolate(varsLs$maf)
     MesKit::JaccardIndex(maf,type = input$JItype)
      dev.off()
    },
    contentType = paste('image/',input$DownloadJaccardIndexCheck,sep="")
  )
  output$DownloadSharedPlot <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadSharedPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadSharedPlotCheck == "png"){
        png(file,width = input$width2 , height = 560, res = 100)
        
      }
      else if (input$DownloadSharedPlotCheck == "pdf"){
        pdf(file,width = input$width2/100 , height = 6)
      }
      print(msp())
      dev.off()
    },
    contentType = paste('image/',input$DownloadSharedPlotCheck,sep="")
  )
  output$DownloadCCFDensity <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadClonePlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadCCFDensityCheck == "png"){
        png(file,width = input$widthccfden , height = 560,res = 100)
      }
      else if (input$DownloadCCFDensityCheck == "pdf"){
        pdf(file,width = input$widthccfden/100 , height = 6)
      }
      print(ccfden())
      dev.off()
    },
    contentType = paste('image/',input$DownloadCCFDensityCheck,sep="")
  )
  output$DownloadClonePlot <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadClonePlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadClonePlotCheck == "png"){
        png(file,width = input$width5 , height = 560,res = 100)
      }
      else if (input$DownloadClonePlotCheck == "pdf"){
        pdf(file,width = input$width5/100 , height = 6)
      }
      print(clp())
      dev.off()
    },
    contentType = paste('image/',input$DownloadClonePlotCheck,sep="")
  )
  output$DownloadPhyloTree <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadPhyloTreeCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadPhyloTreeCheck == "png"){
        png(file,width = 1000, height = 650,res = 100)
      }
      else if (input$DownloadPhyloTreeCheck == "pdf"){
        pdf(file,width = 10, height = 6.5)
      }
      print(pht())
      dev.off()
    }
  )
  
  output$DownloadGOPlot <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadGOPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadGOPlotCheck == "png"){
        png(file,width = input$width6, height = input$height6,res = 100)
      }
      else if (input$DownloadGOPlotCheck == "pdf"){
        pdf(file,width = input$width6/100, height = input$height6/100)
      }
      print(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
      dev.off()
    },
    contentType = paste('image/',input$DownloadGOPlotCheck,sep="")
  )
  output$DownloadGOTable <- downloadHandler(
    filename = 'Rtable.csv',
    content = function(file){
      data <- GO()[[1]][[which(names(GO()[[1]]) == input$gl)]]
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  output$DownloadPathPlot <- downloadHandler(
    filename = function() {
      paste("Pathwayplot.",input$DownloadPathPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadPathPlotCheck == "png"){
        png(file,width = input$width7,height = input$height7,res = 100)
      }
      else if (input$DownloadPathPlotCheck == "pdf"){
        pdf(file,width = input$width7/100, height = input$height7/100)
      }
      print(Path()[[2]][[which(names(Path()[[2]]) == input$pl)]])
      dev.off()
    },
    contentType = paste('image/',input$DownloadPathPlotCheck,sep="")
  )
  output$DownloadPathTable <- downloadHandler(
    filename = "Rtable.csv",
    content = function(file){
      data <- Path()[[1]][[which(names(Path()[[1]]) == input$pl)]]
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  output$DownloadSignatureSummary <- downloadHandler(
    filename = "Rtable.csv",
    content = function(file){
      data <- sigOFA()
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  output$DownloadSignatureBT <- downloadHandler(
    filename = "Rtable.csv",
    content = function(file){
      data <- sigBT()
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  output$DownloadSigOFATable1 <- downloadHandler(
    filename = "Rtable.csv",
    content = function(file){
      data <- sigOFA1()[[2]][,c(1:2)]
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadSignaturePlot1 <- downloadHandler(
    filename = function() {
      paste("Rplot.", input$DownloadSignaturePlotCheck1, sep='')
    },
    content = function(file) {
      if (input$DownloadSignaturePlotCheck1 == "png"){
        png(file,width = input$widthsig1, height = input$heightsig1,res = 100)
      }
      else if (input$DownloadSignaturePlotCheck1 == "pdf"){
        pdf(file,width = input$widthsig1/100, height = input$heightsig1/100)
      }
      print(sigOFA1()[[1]])
      dev.off()
    },
    contentType = paste('image/',input$DownloadSignaturePlotCheck1,sep="")
  )
  
  output$DownloadSignaturePlot2 <- downloadHandler(
    filename = function() {
      paste("Rplot.",input$DownloadSignaturePlotCheck2, sep='')
    },
    content = function(file) {
      if (input$DownloadSignaturePlotCheck2 == "png"){
        png(file,width = input$widthsig2, height = input$heightsig2,res = 100)
      }
      else if (input$DownloadSignaturePlotCheck2 == "pdf"){
        pdf(file,width = input$widthsig2/100, height = input$heightsig2/100)
      }
      print(sigOFA2())
      dev.off()
    },
    contentType = paste('image/',input$DownloadSignaturePlotCheck2,sep="")
  )
  output$warningMessage00 <- renderUI({
      if(input$submit1){
          maf <- isolate(varsLs$maf)
          colNames <- colnames(maf@data)
          standardCol <- c("Hugo_Symbol","Chromosome","Start_Position","End_Position",
                           "Variant_Classification", "Variant_Type", "Reference_Allele",
                           "Tumor_Seq_Allele2", "VAF", "Tumor_Sample_Barcode")
          is <- intersect(colNames,standardCol)
          if(length(is) != 10){
              tagList(
                  tags$p("Wrong data format!",br(),"Please click the button",
                         tags$img(src = 'image/button.png',width = "40px",height = "40px"),
                         " to see the example file.",
                         style = "color: red;
                          font-size:27px; 
                          font-weight:500;")) 
          }
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage01 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit0){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage02 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit2){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage03 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit3){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage04 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit4){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage05 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit5){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage06 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit6){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage07 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit7){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage08 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submitccfden){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage09 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit8){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage10 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit9){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage11 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submitSig4){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage12 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submitSig2){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage13 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submitSig){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage14 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submitSig1){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  output$warningMessage15 <- renderUI({
      maf <- isolate(varsLs$maf)
      if(is.null(maf)&input$submit10){
          tagList(
              tags$p("Please upload data in session 'Input Data'!",
                     style = "color: red;
                          font-size:27px; 
                          font-weight:500;")
          )
      }
      else{
          return(NULL)
      }
  })
  
  
  
  
})  