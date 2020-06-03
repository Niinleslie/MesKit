
# Define server logic required to plot various variables against mpg
shinyServer(function(input, output, session){
  observeEvent(input$help, {
    updateTabItems(session, "sidername", "home")
  })
  observeEvent(input$contact, {
    updateTabItems(session, "sidername", "home")
  })
  
  
 
  
  inputData <- eventReactive(input$submit1, {
    if(input$submit1){
      
      
      if(is.null(input$mafFile)){
        mafFile <- system.file("extdata", "HCC6046.maf", package = "MesKit")
        ccfFile <- system.file("extdata", "HCC6046.ccf.tsv", package = "MesKit")
        maf <- readMaf(mafFile = mafFile,ccfFile = ccfFile)
      } else {
        if(!is.null(input$mafFile) & !is.null(input$ccfFile)){
          maf <- readMaf(mafFile = input$mafFile$datapath,
                         ccfFile =  input$ccfFile$datapath,
                         refBuild = input$ref)          
        }
        else{
          maf <-  readMaf(mafFile = input$mafFile$datapath, refBuild = input$ref)
        }
      }
      return(maf)
    }
  })
  
  varsMaf <- reactiveValues()
  
  observeEvent(input$submit1,{
    withProgress(min = 0, max = 1, value = 0, {
      ## Rshiny: progress bar
      setProgress(message = 'Input data: Generating ', detail = paste("Maf/MafList Class", sep="")) 
      varsMaf[['maf']] <- inputData()
      incProgress(amount=1)
      
      
      setProgress(message = paste("Input data: Maf/MafList Generation Done!", sep=""), detail = "") 
      Sys.sleep(1)
      
    })
  })
  
  buttonValue <- reactiveValues(a = 0, b = 0)
  observeEvent(input$iecontrol01,{
    buttonValue$a <- buttonValue$a + 1
    if(buttonValue$a == 2){
      buttonValue$a <- 0
    }
    buttonValue$b <- 0
  })
  observeEvent(input$iecontrol02,{
    buttonValue$a <- 0
    buttonValue$b <- buttonValue$b + 1
    if(buttonValue$b == 2){
      buttonValue$b <- 0
    }
  })
  ## output Introduction of maf datatable
  output$ie1 <- renderUI({
    if(buttonValue$a == 1){
      tagList(
          div(
            h3(strong("The MAF files")),
            p("MAF files contain many fields of information about chromosome and gene mutations and 
              their annotations. The following fields are highly recommended to be contained in
              the MAF files.",
              style = "font-size:16px; font-weight:500;line-height:30px;"),
            p(strong("Mandatory fields"),": Hugo_Symbol,Chromosome,Start_Position,End_Position,
                        Variant_Classification, Variant_Type,Reference_Allele,
                        Tumor_Seq_Allele2,Ref_allele_depth,Alt_allele_depth,
                        VAF,Tumor_Sample_Barcode,Patient_ID,Tumor_ID",
              style = "font-size:16px; font-weight:500;line-height:30px;"),
            h3(strong("Example MAF file"))
          ),
          DT::dataTableOutput("ied1",width = "100%"),
          br()
      )
    }
  })
  output$ied1 <- renderDataTable({
    if(input$iecontrol01){
      mafFile <- system.file("extdata/", "HCC_LDC.maf", package = "MesKit")
      maf_data <- data.table::fread(
          file = mafFile,
          quote = "",
          header = TRUE,
          data.table = TRUE,
          fill = TRUE,
          sep = '\t',
          skip = "Hugo_Symbol",
          stringsAsFactors = FALSE
      )
      d <- datatable(maf_data, options = list(searching = TRUE, pageLength = 5, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
      return(d)
    }
  })
  ## output Introduction of CCF
  output$ie2 <- renderUI({
    if(buttonValue$b == 1){
        tagList(
          h3(strong("The CCF file")),
          p("CCF files contain cancer cell fraction of each mutation.",
            style = "font-size:16px; font-weight:500;line-height:30px;"),
          p(strong("Mandatory fields"),": Patient_ID,Tumor_Sample_Barcode,Chromosome,Start_Position,End_Position",
            style = "font-size:16px; font-weight:500;line-height:30px;"),
          h3(strong("Example CCF file:")),
          DT::dataTableOutput("ied2"),
          br()
        )
    }
  })
  output$ied2 <- renderDataTable({
    if(input$iecontrol02){
     ccfFile <- system.file("extdata/", "HCC_LDC.ccf.tsv", package = "MesKit")
     ccf_data <- suppressWarnings(data.table::fread(
         ccfFile,
         quote = "",
         header = TRUE,
         fill = TRUE,
         sep = '\t',
         stringsAsFactors = FALSE
     ))
      d <- datatable(ccf_data, options = list(searching = TRUE, pageLength = 5, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
      return(d)
    }
  })
  
  output$datapreview <- renderUI({
      if(input$submit1){
        DT::dataTableOutput('maftable', width = '100%')
      }
  })
  output$maftable <- DT::renderDataTable({
    if(input$submit1){
        maf <- inputData()
        if(class(maf) == "Maf"){
            t <- datatable(inputData()@data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)  
        }else{
            t <- datatable(inputData()[[1]]@data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)  
        }
        return(t)
    }
  })

  ms <- eventReactive(input$submit_mathscore, {
    maf <- isolate(varsMaf$maf)
    validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
    )
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'MATH Score: Calculation in progress',
                 detail = 'This may take a while...')
    mathScore(maf, min.vaf = as.numeric(input$mathscore_minvaf),
              withinTumor = input$mathscore_withintumor)
  })
  output$mathScore <- DT::renderDataTable({
    ms()
  })
  
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
  output$DownloadMathScore <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          data <- ms()
          write.csv(data,file,row.names = F)
      },
      contentType = 'text/csv'
  )
  
  vafcluster <- eventReactive(input$submit3, {
      maf <- isolate(varsMaf$maf)
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      
      withProgress(min = 0, max = 1, value = 0,{
          setProgress(message = 'vafCluster: Calculation in progress',
                      detail = 'This may take a while...')
          vc <- vafCluster(maf, 
                           withinTumor = input$vafcluster_withintumor,
                           segCN.file = input$vafcluster_segfile,
                           min.vaf = as.numeric(input$vafcluster_minvaf) ,
                           max.vaf = as.numeric(input$vafcluster_maxvaf))        
          incProgress(amount = 1)
          setProgress(message = 'vafCluster: Calculation done!')
      })
          
        return(vc)
  })
  

  
  getpatient.vafcluster <- eventReactive(input$vafcluster.patientlist, {
      return(input$vafcluster.patientlist)
  })
  
  output$vafcluster.patientlist <- renderUI({
      if(!is.null(vafcluster())){
          if(!"cluster.plot" %in% names(vafcluster())){
              names <- names(vafcluster())
              selectInput("vafcluster.patientlist", "Patient",
                          choices = names, width = 600) 
          }
      }
  })
  
  output$vafcluster.samplelist <- renderUI({
      if(!is.null(vafcluster())){
          if("cluster.plot" %in% names(vafcluster())){
              sample.list <- vafcluster()$cluster.plot
              names <- names(sample.list) 
          }else{
              sample.list <- vafcluster()[[getpatient.vafcluster()]]$cluster.plot
              names <- names(sample.list)  
              print(names)
          }
          if(input$vafcluster_withintumor){
              tagList(
                  selectInput("vafcluster.sl", "Tumor ID",
                              choices = names, width = 600) 
              )
          }else{
              tagList(
                  selectInput("vafcluster.sl", "Tumor Sample Barcode",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getsample.vafcluster <- eventReactive(input$vafcluster.sl,{
      return(input$vafcluster.sl)
  })
  
  vafcluster_width <- reactive({
      return(input$vafcluster_width)
  })
  vafcluster_height <- reactive({
      return(input$vafcluster_height)
  })
  
  output$vaf <- renderPlot({
      if(!is.null(vafcluster())){
          if("cluster.plot" %in% names(vafcluster())){
              return(vafcluster()$cluster.plot[[getsample.vafcluster()]])
          }else{
              return(vafcluster()[[getpatient.vafcluster()]]$cluster.plot[[getsample.vafcluster()]])
          }
      }
  },width = vafcluster_width,
  height = vafcluster_height,
  res = 100)
  
  output$vcdb <- renderUI({
    if(!is.null(vafcluster())){
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
  
  output$DownloadVafPlot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$DownloadVafPlotCheck, sep='')
      },
      content = function(file) {
          if (input$DownloadVafPlotCheck == "png"){
              png(file,width = input$vafcluster_width , height = input$vafcluster_height,res = 100)
          }
          else if (input$DownloadVafPlotCheck == "pdf"){
              pdf(file,width = input$vafcluster_width/100 , height = input$vafcluster_height/100)
          }
          if("cluster.plot" %in% names(vafcluster())){
              print(vafcluster()$cluster.plot[[getsample.vafcluster()]])
          }else{
              print(vafcluster()[[getpatient.vafcluster()]]$cluster.plot[[getsample.vafcluster()]])
          }
          dev.off()
      },
      contentType = paste('image/',input$DownloadVafPlotCheck,sep="")
  )
  
  output$vafcluster_table <- DT::renderDataTable({
      if(!is.null(vafcluster())){
          if("cluster.data" %in% names(vafcluster())){
              t <- vafcluster()$cluster.data
          }else{
              t <- vafcluster()[[getpatient.vafcluster()]]$cluster.data
          }
          dt <- datatable(t, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5) 
          return(dt)
      }
  })
  
  
  output$vafcluster_table_ui <- renderUI({
      if(!is.null(vafcluster())){
          tagList(
              h4(strong('Cluster result')),
              br(),
              DT::dataTableOutput('vafcluster_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_vafcluster_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$Download_vafcluster_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if("cluster.data" %in% names(vafcluster())){
              t <- vafcluster()$cluster.data
          }else{
              t <- vafcluster()[[getpatient.vafcluster()]]$cluster.data
          }
          write.csv(t,file,row.names = F)
      },
      contentType = 'text/csv'
  )
  
  ccfauc <- eventReactive(input$submit_ccfauc,{
      maf <- varsMaf$maf
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      
      withProgress(min = 0, max = 1, value = 0, {
          setProgress(message = 'ccfAUC: Calculation in progress',
                      detail = 'This may take a while...')
          validate(
              need(!(is.null(maf)), "")
          )
          cc <- ccfAUC(maf,
                       min.ccf = as.numeric(input$ccfauc_minccf) ,
                       withinTumor = input$ccfauc_withintumor)
          incProgress(amount = 1)
          setProgress(message = 'ccfAUC: Calculation done!')
      })
      return(cc)
  })
  
  output$ccfauc.patientlist <- renderUI({
      if(!is.null(ccfauc())){
          if(!"CCF.density.plot" %in% names(ccfauc())){
              cc <- ccfauc()
              names <- names(cc)
              tagList(
                  selectInput("auc.pl", "Patient",
                              choices = names, width = 600) 
              ) 
          }
      }
  })
  
  getpatient.ccfauc <- eventReactive(input$auc.pl,{
      return(input$auc.pl)
  })
  
  ccfauc_width <- reactive({
      return(input$ccfauc_width)
  })
  ccfauc_height <- reactive({
      return(input$ccfauc_height)
  })
  
  output$ccfauc_plot <- renderPlot({
      if(!is.null(ccfauc())){
          if("CCF.density.plot" %in% names(ccfauc())){
              return(ccfauc()$CCF.density.plot)
          }else{
              return(ccfauc()[[getpatient.ccfauc()]]$CCF.density.plot) 
          }
      }
  },  
  width = ccfauc_width,
  height = ccfauc_height,
  res = 100)
  
  output$ccfauc_db_ui <- renderUI({
      if(!is.null(ccfauc())){
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
  
  output$Download_ccfauc_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_ccfauc_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_ccfauc_plot_check == "png"){
              png(file,width = input$ccfauc_width , height = input$ccfauc_height,res = 100)
          }
          else if (input$Download_ccfauc_plot_check == "pdf"){
              pdf(file,width = input$ccfauc_width/100 , height = input$ccfauc_height/100)
          }
          if("CCF.density.plot" %in% names(ccfauc())){
              print(ccfauc()$CCF.density.plot)
          }else{
              print(ccfauc()[[getpatient.ccfauc()]]$CCF.density.plot) 
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_ccfauc_plot_check,sep="")
  )
  
  output$ccfauc_table <- DT::renderDataTable({
      if(!is.null(ccfauc())){
          if("CCF.density.plot" %in% names(ccfauc())){
             t <- ccfauc()$AUC.value
          }else{
              t <- ccfauc()[[getpatient.ccfauc()]]$AUC.value
          }
          dt <- datatable(t,
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
  
  output$Download_ccfauc_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if("CCF.density.plot" %in% names(ccfauc())){
              t <- ccfauc()$AUC.value
          }else{
              t <- ccfauc()[[getpatient.ccfauc()]]$AUC.value
          }
          write.csv(t,file,row.names = F)
      },
      contentType = 'text/csv'
  )
  
  ## calfst sever
  calfst <- eventReactive(input$submit_calfst,{
      maf <- varsMaf$maf
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      withProgress(min = 0, max = 1, value = 0, {
          setProgress(message = 'calFst: Calculation in progress',
                      detail = 'This may take a while...')
          if(input$calfst_title == ""){
              title <- NULL
          }else{
              title <- input$calfst_title
          }
          fst <- calFst(maf,
                        min.vaf = as.numeric(input$calfst_minvaf) ,
                        min.total.depth = as.numeric(input$calfst_mintotaldepth),
                        withinTumor =input$calfst_withinTumor,
                        title = title,
                        use.circle = input$calfst_usecircle,
                        number.cex = as.numeric(input$calfst_numbercex),
                        number.col = input$calfst_numbercol)
          incProgress(amount = 1)
          setProgress(message = 'calFst: Calculation done!')
      })
      return(fst)
  })
  
  output$calfst.patientlist <- renderUI({
      if(!is.null(calfst())){
          if(!"Fst.plot" %in% names(calfst())){
              names <- names(calfst())
              tagList(
                  selectInput("calfst.pl", "Patient",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getpatient.calfst <- eventReactive(input$calfst.pl,{
      return(input$calfst.pl)
  })
  
  calfst_width <- reactive({
      return(input$calfst_width)
  })
  calfst_height <- reactive({
      return(input$calfst_height)
  })
  
  output$calfst_plot <- renderPlot({
      if(!is.null(calfst())){
          if(!"Fst.plot" %in% names(calfst())){
              return(calfst()[[getpatient.calfst()]]$Fst.plot)
          }else{
              return(calfst()$Fst.plot)
          }
      }
  },  
  width = calfst_width,
  height = calfst_height,
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
  
  
  output$Download_calfst_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_calfst_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_calfst_plot_check == "png"){
              png(file,width = input$calfst_width , height = input$calfst_height,res = 100)
          }
          else if (input$Download_calfst_plot_check == "pdf"){
              pdf(file,width = input$calfst_width/100 , height = input$calfst_height/100)
          }
          if(!"Fst.plot" %in% names(calfst())){
              print(calfst()[[getpatient.calfst()]]$Fst.plot)
          }else{
               print(calfst()$Fst.plot)
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_calfst_plot_check,sep="")
  )
  
  output$calfst_pair_table <- DT::renderDataTable({
      if(!is.null(calfst())){
          if(!"Fst.plot" %in% names(calfst())){
              m <- calfst()[[getpatient.calfst()]]$Fst.pair
          }else{
              m <- calfst()$Fst.pair
          }
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
  
  output$Download_calfst_pair_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if(!"Fst.plot" %in% names(calfst())){
              m <- calfst()[[getpatient.calfst()]]$Fst.pair
          }else{
              m <- calfst()$Fst.pair
          }
          rownames(m) <- colnames(m)
          m <- as.data.frame(m)
          colnames(m) <- rownames(m)
          write.csv(m,file,row.names = TRUE)
      },
      contentType = 'text/csv'
  )
  
  
  ## neidist sever
  calneidist <- eventReactive(input$submit_calneidist,{
      maf <- varsMaf$maf
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      withProgress(min = 0, max = 1, value = 0, {
          setProgress(message = 'calNeiDist: Calculation in progress',
                      detail = 'This may take a while...')
          if(input$calneidist_title==""){
              title <- NULL
          }else{
              title <- input$calneidist_title
          }
          cc <- calNeiDist(maf, 
                           min.ccf = as.numeric(input$calneidist_minccf) ,
                           title = title,
                           withinTumor = input$calneidist_withintumor,
                           use.circle = input$calneidist_usecircle,
                           number.cex = as.numeric(input$calneidist_numbercex),
                           number.col = input$calneidist_numbercol)
          incProgress(amount = 1)
          setProgress(message = 'calNeiDist: Calculation done!')
      })
      return(cc)
  })
  
  output$calneidist.patientlist <- renderUI({
      if(!is.null(calneidist())){
          if(!"Nei.plot" %in% names(calneidist())){
              names <- names(calneidist())
              tagList(
                  selectInput("calneidist.pl", "Patient",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getpatient.calneidist <- eventReactive(input$calneidist.pl,{
      return(input$calneidist.pl)
  })
  
  calneidist_width <- reactive({
      return(input$calneidist_width)
  })
  calneidist_height <- reactive({
      return(input$calneidist_height)
  })
  
  output$calneidist_plot <- renderPlot({
      if(!is.null(calneidist())){
          if(!"Nei.plot" %in% names(calneidist())){
              return(calneidist()[[getpatient.calneidist()]]$Nei.plot)
          }else{
              return(calneidist()$Nei.plot)
          }
      }
  },  
  width = calneidist_width,
  height = calneidist_height,
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
  
  output$Download_calneidist_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_calneidist_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_calneidist_plot_check == "png"){
              png(file,width = input$calneidist_width , height = input$calneidist_height,res = 100)
          }
          else if (input$Download_calneidist_plot_check == "pdf"){
              pdf(file,width = input$calneidist_width/100 , height = input$calneidist_height/100)
          }
          if(!"Nei.plot" %in% names(calneidist())){
              print(calneidist()[[getpatient.calneidist()]]$Nei.plot)
          }else{
              print(calneidist()$Nei.plot)
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_calneidist_plot_check,sep="")
  )
  
  output$calneidist_pair_table <- DT::renderDataTable({
      if(!is.null(calneidist())){
          if(!"Nei.plot" %in% names(calneidist())){
              m <- calneidist()[[getpatient.calneidist()]]$Nei.dist
          }else{
              m <- calneidist()$Nei.dist
          }
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
              h4(strong("Nei's distance")),
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
  
  output$Download_calneidist_pair_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if(!"Nei.plot" %in% names(calneidist())){
              m <- calneidist()[[getpatient.calneidist()]]$Nei.dist
          }else{
              m <- calneidist()$Nei.dist
          }
          rownames(m) <- colnames(m)
          m <- as.data.frame(m)
          colnames(m) <- rownames(m)
          write.csv(m,file,row.names = TRUE)
      },
      contentType = 'text/csv'
  )
  
  ## mutheatmap
  mutheatmap <- eventReactive(input$submit_mutheatmap,{
      maf <- varsMaf$maf
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      if(!is.null(input$mutheatmap_genelist$datapath)){
          genelist <- as.character(read.table(input$mutheatmap_genelist$datapath)$V1)
      }else{
          genelist <- NULL
      }
      withProgress(min = 0, max = 1, value = 0, {
          setProgress(message = 'mutHeatmap: Calculation in progress',
                      detail = 'This may take a while...')
          
          hm <- mutHeatmap(maf,
                           geneList = genelist,
                           min.vaf = as.numeric(input$mutheatmap_minvaf),
                           min.ccf = as.numeric(input$mutheatmap_minccf),
                           use.ccf = input$mutheatmap_useccf,
                           plot.geneList = input$mutheatmap_plotgenelist,
                           show.gene = input$mutheatmap_showgene,
                           show.geneList = input$mutheatmap_showgenelist,
                           mut.threshold = as.numeric(input$mutheatmap_mutthreshold),
                           sample.text.size = as.numeric(input$mutheatmap_sampletextsize),
                           legend.title.size = as.numeric(input$mutheatmap_legendtitlesize),
                           gene.text.size = as.numeric(input$mutheatmap_genetextsize))
          incProgress(amount = 1)
          setProgress(message = 'mutHeatmap: Calculation done!')
      })
      return(hm)
  })
  
  output$mutheatmap.patientlist <- renderUI({
      if(!is.null(mutheatmap())){
          if(!identical(c("gg","ggplot"),class(mutheatmap()))){
              names <- names(mutheatmap())
              tagList(
                  selectInput("mutheatmap.pl", "Patient",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getpatient.mutheatmap <- eventReactive(input$mutheatmap.pl,{
      return(input$mutheatmap.pl)
  })
  
  mutheatmap_width <- reactive({
      return(input$mutheatmap_width)
  })
  mutheatmap_height <- reactive({
      return(input$mutheatmap_height)
  })
  
  output$mutheatmap_plot <- renderPlot({
      if(!is.null(mutheatmap())){
          if(!identical(c("gg","ggplot"),class(mutheatmap()))){
              return(mutheatmap()[[getpatient.mutheatmap()]]) 
          }else{
              return(mutheatmap())
          }
          
      }
  },  
  width = mutheatmap_width,
  height = mutheatmap_height,
  res = 100)
  
  
  output$mutheatmap_db_ui <- renderUI({
      if(!is.null(mutheatmap())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_mutheatmap_plot_check', 
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
                  downloadBttn('Download_mutheatmap_plot', 'Download')
              )
          )
      }
  })
  
  output$Download_mutheatmap_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_mutheatmap_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_mutheatmap_plot_check == "png"){
              png(file,width = input$mutheatmap_width , height = input$mutheatmap_height,res = 100)
          }
          else if (input$Download_mutheatmap_plot_check == "pdf"){
              pdf(file,width = input$mutheatmap_width/100 , height = input$mutheatmap_height/100)
          }
          if(!identical(c("gg","ggplot"),class(mutheatmap()))){
              print(mutheatmap()[[getpatient.mutheatmap()]]) 
          }else{
              print(mutheatmap())
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_mutheatmap_plot_check,sep="")
  )
   
  ## comparejsi sever
  comparejsi <- eventReactive(input$submit_comparejsi,{
      maf <- varsMaf$maf
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      withProgress(min = 0, max = 1, value = 0, {
          setProgress(message = 'compareJSI : Calculation in progress',
                      detail = 'This may take a while...')

          if(input$comparejsi_title==""){
              title <- NULL
          }else{
              title <- input$comparejsi_title
          }
          cc <- compareJSI(maf, 
                           min.vaf = as.numeric(input$comparejsi_minccf) ,
                           pairByTumor = input$comparejsi_pairbytumor,
                           title = title,
                           use.circle = input$comparejsi_usecircle,
                           number.cex = as.numeric(input$comparejsi_numbercex),
                           number.col = input$comparejsi_numbercol)
          incProgress(amount = 1)
          setProgress(message = 'compareJSI: Calculation done!')
      })
      return(cc)
  })
  
  
  
  output$comparejsi.patientlist <- renderUI({
      if(!is.null(comparejsi())){
          if(!"JSI.plot" %in% names(comparejsi())){
              names <- names(comparejsi())
              tagList(
                  selectInput("comparejsi.pl", "Patient",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getpatient.comparejsi <- eventReactive(input$comparejsi.pl,{
      return(input$comparejsi.pl)
  })
  
  comparejsi_width <- reactive({
      return(input$comparejsi_width)
  })
  comparejsi_height <- reactive({
      return(input$comparejsi_height)
  })
  
  output$comparejsi_plot <- renderPlot({
      if(!is.null(comparejsi())){
          if(!"JSI.plot" %in% names(comparejsi())){
              return(comparejsi()[[getpatient.comparejsi()]]$JSI.plot)
          }else{
              return(comparejsi()$JSI.plot)
          }
      }
  },  
  width = comparejsi_width,
  height = comparejsi_height,
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
  
  output$Download_comparejsi_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_comparejsi_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_comparejsi_plot_check == "png"){
              png(file,width = input$comparejsi_width , height = input$comparejsi_height,res = 100)
          }
          else if (input$Download_comparejsi_plot_check == "pdf"){
              pdf(file,width = input$comparejsi_width/100 , height = input$comparejsi_height/100)
          }
          if(!"JSI.plot" %in% names(comparejsi())){
              print(comparejsi()[[getpatient.comparejsi()]]$JSI.plot)
          }else{
              print(comparejsi()$JSI.plot)
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_comparejsi_plot_check,sep="")
  )
  
  output$comparejsi_pair_table <- DT::renderDataTable({
      if(!is.null(comparejsi())){
          if(!"JSI.plot" %in% names(comparejsi())){
              m <- comparejsi()[[getpatient.comparejsi()]]$JSI.pair
          }else{
              m <- comparejsi()$JSI.pair
          }
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
  
  
  output$comparejsi_pair_table_ui <- renderUI({
      if(!is.null(comparejsi())){
          tagList(
              h4(strong("JSI pair")),
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
  
  output$Download_comparejsi_pair_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if(!"JSI.plot" %in% names(comparejsi())){
              m <- comparejsi()[[getpatient.comparejsi()]]$JSI.pair
          }else{
              m <- comparejsi()$JSI.pair
          }
          rownames(m) <- colnames(m)
          m <- as.data.frame(m)
          colnames(m) <- rownames(m)
          write.csv(m,file,row.names = TRUE)
      },
      contentType = 'text/csv'
  )
  
  ## classifymut sever
  classifymut <- eventReactive(input$submit_classifymut,{
      withProgress(min = 0, max = 2, value = 0, {
          setProgress(message = 'Mutational lanscape : Calculation in progress',
                      detail = 'This may take a while...')
          maf <- varsMaf$maf
          validate(
              need(!(is.null(maf)), "")
          )
          cc <- classifyMut(maf)
          incProgress(amount = 1)
          setProgress(message = 'Mutational lanscape: Calculation done!')
      })
      return(cc)
  })
  
  output$classifymut_plot <- renderPlot({
      if(!is.null(classifymut())){
          return(classifymut()$mut.profile.plot)
      }
  },  
  width = 560,
  height = 560,
  res = 100)
  
  output$classifymut_download_button_ui <- renderUI({
      if(!is.null(classifymut())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_classifymut_plot_check', 
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
                  downloadBttn('Download_classifymut_plot', 'Download')
              )
          )
      }
  })
  
  output$classifymut_table <- DT::renderDataTable({
      if(!is.null(classifymut())){
          t <- classifymut()$mut.class
          dt <- datatable(t,
                          options = list(searching = TRUE,
                                         pageLength = 10, 
                                         scrollX = TRUE,
                                         dom = "t",
                                         fixedHeader = TRUE),
                          rownames = F) 
          return(dt)
      }
  })
  
  
  output$classifymut_table_ui <- renderUI({
      if(!is.null(classifymut())){
          tagList(
              h4(strong('Mutation class')),
              br(),
              DT::dataTableOutput('classifymut_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_classifymut_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$plotmutprofile.patientlist <- renderUI({
      maf <- varsMaf$maf
      if(!is.null(maf)){
          if(class(maf) == "MafList"){
              patient.list <- names(maf) 
              tagList(
                  selectizeInput("plotmutprofile_patientid",
                                 label = div(style = "font-size:1.5em; font-weight:600;  ", "Select patients"),
                                 choices = patient.list,
                                 multiple = TRUE),
                  bsTooltip(id = "plotmutprofile_patientid",
                            title = 'Select the specific patients. Default: NULL, all patients are included',
                            placement = "top",
                            trigger = "hover"),
              )
          }
      }
  })
  
  plotmutprofile <- eventReactive(input$submit_plotmutprofile,{
      maf <- varsMaf$maf
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      withProgress(min = 0, max = 1, value = 0, {
          setProgress(message = 'plotMutProfile: Calculation in progress',
                      detail = 'This may take a while...')
          if(is.null(input$plotmutprofile_patientid)){
              patientid <- NULL
          }else{
              patientid <- input$plotmutprofile_patientid
          }
          cc <- plotMutProfile(maf,
                               patient.id = patientid,
                               classByTumor = input$plotmutprofile_classByTumor,
                               class = input$plotmutprofile_class,
                               topGenesCount = input$plotmutprofile_topGenesCount,
                               remove_empty_columns = input$plotmutprofile_remove_empty_columns,
                               remove_empty_rows = input$plotmutprofile_remove_empty_rows)
          incProgress(amount = 1)
          setProgress(message = 'plotMutProfile: Calculation done!')
      })
      return(cc)
  })
  
  
  getpatient.plotmutprofile <- eventReactive(input$plotmutprofile.pl,{
      return(input$plotmutprofile.pl)
  })
  
  plotmutprofile_width <- reactive({
      return(input$plotmutprofile_width)
  })
  plotmutprofile_height <- reactive({
      return(input$plotmutprofile_height)
  })
  
  output$plotmutprofile_plot <- renderPlot({
      if(!is.null(plotmutprofile())){
          return(plotmutprofile())
      }
  },  
  width = plotmutprofile_width,
  height = plotmutprofile_height,
  res = 100)
  
  
  output$plotmutprofile_download_button_ui <- renderUI({
      if(!is.null(plotmutprofile())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_plotmutprofile_plot_check', 
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
                  downloadBttn('Download_plotmutprofile_plot', 'Download')
              )
          )
      }
  })
  
  output$Download_plotmutprofile_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_plotmutprofile_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_plotmutprofile_plot_check == "png"){
              png(file,width = input$plotmutprofile_width , height = input$plotmutprofile_height,res = 100)
          }
          else if (input$Download_plotmutprofile_plot_check == "pdf"){
              pdf(file,width = input$plotmutprofile_width/100 , height = input$plotmutprofile_height/100)
          }
          print(plotmutprofile())
          dev.off()
      },
      contentType = paste('image/',input$Download_plotmutprofile_plot_check,sep="")
  )
  
  ## plotCNA sever 
  output$plotcna.patientlist <- renderUI({
      
      if(!is.null(input$plotcna_segfile$datapath)){
          seg <- readSegment(segCN.file = input$plotcna_segfile$datapath,
                             gisticAllLesionsFile = input$plotcna_gisticAllLesionsFile,
                             gisticAmpGenesFile = input$plotcna_gisticAmpGenesFile,
                             gisticDelGenesFile = input$plotcna_gisticDelGenesFile,
                             gistic.qval = input$plotcna_gisticqval)
          patient.list <- unique(seg$Patient_ID)
           tagList(
                  selectizeInput("plotcna_patientid",
                                 label = div(style = "font-size:1.5em; font-weight:600;  ", "Select patients"),
                                 choices = patient.list,
                                 multiple = TRUE),
                  bsTooltip(id = "plotcna_patientid",
                            title = 'Select the specific patients. Default: NULL, all patients are included',
                            placement = "top",
                            trigger = "hover"),
              )
      }
  })
  plotcna <- eventReactive(input$submit_plotcna,{
      
      validate(
          need(!is.null(input$plotcna_segfile$datapath),
               "PlotCNA need copy number information,upload seg file first"
                    )
          )
      
      withProgress(min = 0, max = 2, value = 0, {
          setProgress(message = 'plotCNA : Calculation in progress',
                      detail = 'This may take a while...')
          seg <- readSegment(segCN.file = input$plotcna_segfile$datapath,
                             gisticAllLesionsFile = input$plotcna_gisticAllLesionsFile,
                             gisticAmpGenesFile = input$plotcna_gisticAmpGenesFile,
                             gisticDelGenesFile = input$plotcna_gisticDelGenesFile,
                             gistic.qval = input$plotcna_gisticqval)
          
          incProgress(amount = 1)
          setProgress(message = 'Complete reading seg file!')
          if(is.null(input$plotcna_patientid)){
              patientid <- NULL
          }else{
              patientid <- input$plotcna_patientid
          }
          cna.plot <- plotCNA(seg,
                              patient.id = patientid,
                              refBuild = input$plotcna_refBuild,
                              show.GISTIC.gene = input$plotcna_showGISTICgene,
                              sample.text.size = as.numeric(input$plotcna_sampletextsize) ,
                              sample.bar.height = as.numeric(input$plotcna_samplebarheight) ,
                              legend.text.size = as.numeric(input$plotcna_legendtextsize) ,
                              legend.title.size = as.numeric(input$plotcna_legendtitlesize) ,
                              chrom.bar.height = as.numeric(input$plotcna_chrombarheight) )
          if(!is.null(patientid)){
              seg <- seg[Patient_ID %in% patientid]
          }
          incProgress(amount = 1)
          setProgress(message = 'PlotCNA done!')
      })
      return(list(seg = seg, cna.plot = cna.plot))
  })
  
  plotcna_width <- reactive({
      return(input$plotcna_width)
  })
  plotcna_height <- reactive({
      return(input$plotcna_height)
  })
  
  output$plotcna_plot <- renderPlot({
      if(!is.null(plotcna())){
          return(plotcna()$cna.plot)
      }
  },  
  width = plotcna_width,
  height = plotcna_height,
  res = 100)
  
  output$plotcna_download_button_ui <- renderUI({
      if(!is.null(plotcna())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_plotcna_plot_check', 
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
                  downloadBttn('Download_plotcna_plot', 'Download')
              )
          )
      }
  })
  
  output$Download_plotcna_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_plotcna_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_plotcna_plot_check == "png"){
              png(file,width = input$plotcna_width , height = input$plotcna_height,res = 100)
          }
          else if (input$Download_plotcna_plot_check == "pdf"){
              pdf(file,width = input$plotcna_width/100 , height = input$plotcna_width/100)
          }
          print(plotcna()$cna.plot)
          dev.off()
      },
      contentType = paste('image/',input$Download_plotcna_plot_check,sep="")
  )
  
  output$plotcna_table <- DT::renderDataTable({
      if(!is.null(plotcna())){
          t <- plotcna()$seg
          d <- datatable(t, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)  
          return(d)
      }
  })
  
  
  output$plotcna_table_ui <- renderUI({
      if(!is.null(plotcna())){
          tagList(
              h4(strong('Seg file ')),
              br(),
              DT::dataTableOutput('plotcna_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_plotcna_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$Download_plotcna_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          data <- plotcna()$seg
          write.csv(data,file,row.names = F)
      },
      contentType = 'text/csv'
  )
  
  ## testneutral sever
  testneutral <- eventReactive(input$submit_testneutral,{
      maf <- varsMaf$maf
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      withProgress(min = 0, max = 1, value = 0, {
          setProgress(message = 'testNeutral : Calculation in progress',
                      detail = 'This may take a while...')
          t <- testNeutral(maf,
                            min.vaf = as.numeric(input$testneutral_minvaf),
                            max.vaf = as.numeric(input$testneutral_maxvaf),
                           min.total.depth = as.numeric(input$testneutral_mintotaldepth),
                           withinTumor = input$testneutral_withintumor,
                           R2.threshold = as.numeric(input$testneutral_R2threshold),
                           min.mut.count = as.numeric(input$testneutral_minmutcount))
          incProgress(amount = 1)
          setProgress(message = 'testNeutral: Calculation done!')
      })
      validate(
          need(!(is.na(t)), "No result is generated")
      )
      return(t)
  })
  
  output$testneutral.patientlist <- renderUI({
      if(!is.null(testneutral())){
          if(!"neutrality.metrics" %in% names(testneutral())){
              patient.list <- testneutral()
              names <- names(patient.list)
              tagList(
                  selectInput("testneutral.pl", "Patient",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getpatient.testneutral <- eventReactive(input$testneutral.pl,{
      return(input$testneutral.pl)
  })
  
  output$testneutral.samplelist <- renderUI({
      if(!is.null(testneutral())){
          if(!"neutrality.metrics" %in% names(testneutral())){
              patient <- input$testneutral.pl
              sample.list <- testneutral()[[patient]]$model.fitting.plot
              names <- names(sample.list)
              tagList(
                  selectInput("testneutral.sl", "Tumor sample barcode",
                              choices = names, width = 600) 
              ) 
          }else{
              sample.list <- testneutral()$model.fitting.plot
              names <- names(sample.list)
              tagList(
                  selectInput("testneutral.sl", "Tumor sample barcode",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getsample.testneutral <- eventReactive(input$testneutral.sl,{
      return(input$testneutral.sl)
  })
  
  testneutral_width <- reactive({
      return(input$testneutral_width)
  })
  testneutral_height <- reactive({
      return(input$testneutral_height)
  })
  
  output$testneutral_plot <- renderPlot({
      if(!is.null(testneutral())){
          patient <- input$testneutral.pl
          sample <- input$testneutral.sl
          if(!"neutrality.metrics" %in% names(testneutral())){
              return(testneutral()[[patient]]$model.fitting.plot[[sample]])
          }else{
              return(testneutral()$model.fitting.plot[[sample]])
          }
      }
  },  
  width = testneutral_width,
  height = testneutral_height,
  res = 100)
  
  output$testneutral_db_ui <- renderUI({
      if(!is.null(testneutral())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_testneutral_plot_check', 
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
                  downloadBttn('Download_testneutral_plot', 'Download')
              )
          )
      }
  })
  
  output$Download_testneutral_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_testneutral_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_testneutral_plot_check == "png"){
              png(file,width = input$testneutral_width , height = input$testneutral_height,res = 100)
          }
          else if (input$Download_testneutral_plot_check == "pdf"){
              pdf(file,width = input$testneutral_width/100 , height = input$testneutral_height/100)
          }
          patient <- input$testneutral.pl
          sample <- input$testneutral.sl
          if(!"neutrality.metrics" %in% names(testneutral())){
              print(testneutral()[[patient]]$model.fitting.plot[[sample]])
          }else{
              print(testneutral()$model.fitting.plot[[sample]])
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_testneutral_plot_check,sep="")
  )
  
  output$testneutral_table <- DT::renderDataTable({
      if(!is.null(testneutral())){
          if(!"neutrality.metrics" %in% names(testneutral())){
              t <- testneutral()[[getpatient.testneutral()]]$neutrality.metrics
          }else{
              t <- testneutral()$neutrality.metrics 
          }
          dt <- datatable(t,
                          options = list(searching = TRUE,
                                         pageLength = 10, 
                                         scrollX = TRUE,
                                         dom = "t",
                                         fixedHeader = TRUE),
                          rownames = F) 
          return(dt)
      }
  })
  
  
  output$testneutral_table_ui <- renderUI({
      if(!is.null(testneutral())){
          tagList(
              h4(strong('Neutrality metrics')),
              br(),
              DT::dataTableOutput('testneutral_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_testneutral_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$Download_testneutral_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if(!"neutrality.metrics" %in% names(testneutral())){
              t <- testneutral()[[getpatient.testneutral()]]$neutrality.metrics
          }else{
              t <- testneutral()$neutrality.metrics 
          }
          write.csv(t,file,row.names = F)
      },
      contentType = 'text/csv'
  )

  ## compareccf sever
  compareccf <- eventReactive(input$submit_compareccf, {
      maf <- isolate(varsMaf$maf)
      validate(
          need(!(is.null(maf)), "Please upload data in 'Input Data'!")
      )
      progress <- Progress$new(session, min=0, max=1)
      on.exit(progress$close())
      progress$set(message = 'compareCCF: Calculation in progress',
                   detail = 'This may take a while...')
      cc <- compareCCF(maf,min.ccf = input$compareccf_minccf,
                       pairByTumor = input$compareccf_pairbytumor)
      progress$set(value = 1)
      return(cc)
  })
  
  output$compareccf.patientlist <- renderUI({
      if(!is.null(compareccf())){
          if(class(compareccf()[[1]]) != "data.frame"){
              patient.list <- compareccf()
              names <- names(patient.list)
              tagList(
                  selectInput("compareccf.pl", "Patient",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  
  getpatient.compareccf <- eventReactive(input$compareccf.pl,{
      return(input$compareccf.pl)
  })
  
  output$compareccf.samplelist <- renderUI({
      if(!is.null(compareccf())){
          if(class(compareccf()[[1]]) != "data.frame"){
              patient <- input$compareccf.pl
              sample.list <- compareccf()[[patient]]
              names <- names(sample.list)
              tagList(
                  selectInput("compareccf.sl", "pairwise",
                              choices = names, width = 600) 
              )
          }else{
              sample.list <- compareccf()
              names <- names(sample.list)
              tagList(
                  selectInput("compareccf.sl", "pairwise",
                              choices = names, width = 600) 
              )
          }
      }
  })
  
  getsample.compareccf <- eventReactive(input$compareccf.sl,{
      return(input$compareccf.sl)
  })
  
  
  output$compareccf_table_ui <- renderUI({
      if(!is.null(compareccf())){
          tagList(
              br(),
              DT::dataTableOutput('compareccf_table'),
              br(),
              fluidRow(
                  column(
                      width = 9
                  ),
                  column(
                      width = 3,
                      downloadBttn('Download_compareccf_table', 'Download')
                  )
              )
          )
      }
  })
  
  output$compareccf_table <- DT::renderDataTable({
      if(!is.null(compareccf())){
          if(class(compareccf()[[1]]) != "data.frame"){
              t <- compareccf()[[getpatient.compareccf()]][[getsample.compareccf()]]
          }else{
              t <- compareccf()[[getsample.compareccf()]]
          }
          dt <- datatable(t,rownames = F) 
          return(dt)
      }
  })
  
  output$Download_compareccf_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if(class(compareccf()[[1]]) != "data.frame"){
              t <- compareccf()[[getpatient.compareccf()]][[getsample.compareccf()]]
          }else{
              t <- compareccf()[[getsample.compareccf()]]
          }
          write.csv(t,file,row.names = F)
      },
      contentType = 'text/csv'
  )
  
  varsphyloTree <- reactiveValues()
  
  observeEvent(input$submit_getphylotree,{
      maf <- isolate(varsMaf$maf)
      validate(
          need(!is.null(maf), "Upload maf file in the section 'Input Data'")
      )
      withProgress(min = 0, max = 1, value = 0, {
          ## Rshiny: progress bar
          setProgress(message = 'Generating ', detail = paste("phyloTree/phyloTreeList Class", sep="")) 
          phyloTree <- getPhyloTree(maf,
                                    min.vaf = as.numeric(input$getphylotree_minvaf),
                                    min.ccf = as.numeric(input$getphylotree_minccf),
                                    bootstrap.rep.num = as.numeric(input$getphylotree_bootstraprepnum))
          varsphyloTree[["phyloTree"]] <- phyloTree
          incProgress(amount=1)
          
          ## Rshiny: progress bar
          # setProgress(message = 'Input data: Generating ', detail = paste("phyloTree from MAF ",sep="")) 
          # varsMaf[['phyloTree']] <-  phyloTree <- getPhyloTree(isolate(varsMaf$maf),method = input$method)
          # incProgress(amount=1)
          
          setProgress(message = paste("Input data: phyloTree/phyloTreeList Generation Done!", sep=""), detail = "") 
          Sys.sleep(1)
          
      })
  })
  
  phylotree <- eventReactive(input$submit_plotphylotree, {
      
      maf <- isolate(varsMaf$maf)
      validate(
          need(!is.null(maf), "Upload maf file in the section 'Input Data'")
      )
      withProgress(min = 0, max = 2, value = 0,{
          
          setProgress(message = 'Generating ', detail = paste("phyloTree/phyloTreeList Class", sep="")) 
          phyloTree <- getPhyloTree(maf,
                                    min.vaf = as.numeric(input$plotphylotree_getphylotree_minvaf),
                                    min.ccf = as.numeric(input$plotphylotree_getphylotree_minccf),
                                    method = input$plotphylotree_getphylotree_method,
                                    bootstrap.rep.num = as.numeric(input$plotphylotree_getphylotree_bootstraprepnum))
          
          incProgress(amount=1)
          setProgress(message = 'Phylogenetic tree: Calculation in progress',
                      detail = 'This may take a while...')
          
          if(input$plotphylotree_branchcol == "NULL"){
              branchCol <- NULL
          }else{
              branchCol <- input$plotphylotree_branchcol
          }
          plot.list <- plotPhyloTree(phyloTree,
                                     branchCol = branchCol,
                                     show.bootstrap = input$plotphylotree_showbootstrap,
                                     use.box = input$plotphylotree_usebox,
                                     min.ratio = as.numeric(input$plotphylotree_minratio) ,
                                     signaturesRef = input$plotphylotree_signatureref,
                                     min.mut.count = as.numeric(input$plotphylotree_minmutcount) )
          incProgress(amount=1)
          setProgress(message = paste("Plot phylotree done!", sep=""), detail = "") 
          
          Sys.sleep(1)
      })
      
      return(plot.list)
  })
  
  output$phylotree.patientlist <- renderUI({
      if(!is.null(phylotree())){
          if(!identical(c("gg","ggplot"),class(phylotree()))){
              names <- names(phylotree())
              selectInput("phylotree.pl", "Patient",
                          choices = names, width = 600)
          }
      }
  })
  
  getpatient.phylotree <- eventReactive(input$phylotree.pl,{
      return(input$phylotree.pl)
  })
  
  plotphylotree_width <- reactive({
      return(input$plotphylotree_width)
  })
  plotphylotree_height <- reactive({
      return(input$plotphylotree_height)
  })
  
  output$phylotree_plot <- renderPlot({
      if(!is.null(phylotree())){
          if(!identical(c("gg","ggplot"),class(phylotree()))){
              return(phylotree()[[getpatient.phylotree()]]) 
          }else{
              return(phylotree()) 
          }
      }  
  },
  width = plotphylotree_width,
  height = plotphylotree_height,
  res = 100
  )
  output$phylotree_downloadbutton_ui<- renderUI({
      if(!is.null(phylotree())){
          br()
          br()
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons('Download_phylotree_check', label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                               choiceNames = list(
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                               ),
                               choiceValues = c("png", "pdf"), inline = T)
              ),
              column(
                  width = 3,
                  downloadBttn('Download_phylotree', 'Download')
              )
          )
      }
  })
  
  output$Download_phylotree <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_phylotree_check, sep='')
      },
      content = function(file) {
          if (input$Download_phylotree_check == "png"){
              png(file,width = input$plotphylotree_width , height = input$plotphylotree_height,res = 100)
          }
          else if (input$Download_phylotree_check == "pdf"){
              pdf(file,width = input$plotphylotree_width/100 , height = input$plotphylotree_height/100)
          }
          if(!identical(c("gg","ggplot"),class(phylotree()))){
              print(phylotree()[[getpatient.phylotree()]]) 
          }else{
              print(phylotree()) 
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_phylotree_check,sep="")
  )

  comparetree <- eventReactive(input$submit_comparetree, {
      
      maf <- isolate(varsMaf$maf)
      validate(
          need(!is.null(maf), "Upload maf file in the section 'Input Data'")
      )
      
      withProgress(min = 0, max = 1, value = 0,{
          setProgress(message = 'Comparetree: Calculation in progress',
                      detail = 'This may take a while...')
          
          phylotree1 <- getPhyloTree(maf, patient.id = input$comparetree_patientid, 
                                     method = input$comparetree_getphylotree_method1,
                                     min.ccf = as.numeric(input$comparetree_getphylotree_minccf1) ,
                                     min.vaf = as.numeric(input$comparetree_getphylotree_minvaf1) ,
                                     bootstrap.rep.num = as.numeric(input$comparetree_getphylotree_bootstraprepnum1))
          phylotree2 <- getPhyloTree(maf, patient.id = input$comparetree_patientid, 
                                     method = input$comparetree_getphylotree_method2,
                                     min.ccf = as.numeric(input$comparetree_getphylotree_minccf2) ,
                                     min.vaf = as.numeric(input$comparetree_getphylotree_minvaf2) ,
                                     bootstrap.rep.num = as.numeric(input$comparetree_getphylotree_bootstraprepnum2) )
          
          ct <- compareTree(phylotree1, phylotree2,
                            common.col = input$comparetree_commoncol,
                            min.ratio = input$comparetree_minratio,
                            show.bootstrap = input$comparetree_showbootstrap,
                            use.box = input$comparetree_usebox,
                            plot = TRUE)
          
          incProgress(amount=1)
          setProgress(message = paste("Comparetree done!", sep=""), detail = "") 
          
          Sys.sleep(1)
      })
      
      return(ct)
  })
  
  output$comparetree.patientlist <- renderUI({
      maf <- varsMaf$maf
      if(!is.null(maf)){
          if(class(maf) == "MafList"){
              patient.list <- names(maf) 
          }else{
              patient.list <- unique(maf@data$Patient_ID)
          }
          tagList(
              selectInput("comparetree_patientid",
                             label = div(style = "font-size:1.6em; font-weight:600;  ", strong("Select patients")),
                             choices = patient.list),
              bsTooltip(id = "comparetree_patientid",
                        title = 'Select the specific patient',
                        placement = "top",
                        trigger = "hover"),
          )
      }
  })
  
  
  
  comparetree_width <- reactive({
      return(input$comparetree_width)
  })
  comparetree_height <- reactive({
      return(input$comparetree_height)
  })
  
  output$comparetree_plot <- renderPlot({
      if(!is.null(comparetree())){
          return(comparetree()$compare.plot)
      }  
  },
  width = comparetree_width,
  height = comparetree_height,
  res = 100
  )
  output$comparetree_dist <- renderPrint({
      if(!is.null(comparetree())){
          return(comparetree()$compare.dist)
      } 
  })
  output$comparetree_db_ui<- renderUI({
      if(!is.null(comparetree())){
          br()
          br()
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons('Download_comparetree_check', label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                               choiceNames = list(
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                               ),
                               choiceValues = c("png", "pdf"), inline = T)
              ),
              column(
                  width = 3,
                  downloadBttn('Download_comparetree', 'Download')
              )
          )
      }
  })
  
  output$Download_comparetree <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_comparetree_check, sep='')
      },
      content = function(file) {
          if (input$Download_comparetree_check == "png"){
              png(file,width = input$comparetree_width , height = input$comparetree_height,res = 100)
          }
          else if (input$Download_comparetree_check == "pdf"){
              pdf(file,width = input$comparetree_width/100 , height = input$comparetree_height/100)
          }
          print(comparetree()$compare.plot)
          dev.off()
      },
      contentType = paste('image/',input$Download_comparetree_check,sep="")
  )
  
  
  treemutsig <- eventReactive(input$submit_treemutsig, {
      
      maf <- isolate(varsMaf$maf)
      validate(
          need(!is.null(maf), "Upload maf file in the section 'Input Data'")
      )
      withProgress(min = 0, max = 4, value = 0,{
          setProgress(message = 'Generating ', detail = paste("phyloTree/phyloTreeList Class", sep="")) 
          
          phyloTree <- getPhyloTree(maf,
                                    method = input$treemutsig_getphylotree_method,
                                    min.vaf = as.numeric(input$treemutsig_getphylotree_minvaf),
                                    min.ccf = as.numeric(input$treemutsig_getphylotree_minccf),
                                    bootstrap.rep.num = as.numeric(input$treemutsig_getphylotree_bootstraprepnum))
          
          incProgress(amount=1)
          
          
          setProgress(message = 'triMatrix: calculation in progress',
                      detail = 'This may take a while...')
          
          tm <- triMatrix(phyloTree, withinTumor = input$treemutsig_withintumor)
          incProgress(amount = 1)
          
          setProgress(message = 'fitSignature: calculation in progress',
                      detail = 'This may take a while...')
          fs <- fitSignatures(tm, signaturesRef = input$treemutsig_signatureref,
                              min.mut.count = as.numeric(input$treemutsig_minmutcount) ,
                              signature.cutoff = as.numeric(input$treemutsig_signaturecutoff))
          incProgress(amount = 1)
          
          setProgress(message = 'plotMutSigProfile: drawing in progress',
                      detail = 'This may take a while...')
          
          if(input$treemutsig_mode == 'NULL'){
              mode <- NULL
          }else{
              mode <- input$treemutsig_mode
          }
          pms <- plotMutSigProfile(fs, mode = mode)
          incProgress(amount = 1)
          
          setProgress(message = 'plotMutSigProfile done!')
          Sys.sleep(1)
          
          return(pms)
          
      })
  })
  
  getpatient.treemutsig <- eventReactive(input$treemutsig.pl,{
      return(input$treemutsig.pl)
  })
  
  output$treemutsig.patientlist <- renderUI({
      if(!is.null(treemutsig())){
          if(class(treemutsig()) == "list"){
              plot.list <- treemutsig()
              names <- names(plot.list)
              selectInput("treemutsig.pl", "Patient_branches",
                          choices = names, width = 600) 
          }
      }
  })
  
  treemutsig_width <- reactive({
      return(input$treemutsig_width)
  })
  treemutsig_height <- reactive({
      return(input$treemutsig_height)
  })
  
  output$treemutsig_plot <- renderPlot({
      if(!is.null(treemutsig())){
          if(class(treemutsig()) == "list"){
              return(treemutsig()[[getpatient.treemutsig()]])
          }else{
              return(treemutsig())
          }
      }
  },
  width = treemutsig_width,
  height = treemutsig_height,
  res = 100
  )
  
  output$treemutsig_download_button_ui <- renderUI({
      if(!is.null(treemutsig())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons(inputId = 'Download_treemutsig_plot_check', 
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
                  downloadBttn('Download_treemutsig_plot', 'Download')
              )
          )
      }
  })
  
  output$Download_treemutsig_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_treemutsig_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_treemutsig_plot_check == "png"){
              png(file,width = input$treemutsig_width , height = input$treemutsig_height,res = 100)
          }
          else if (input$Download_treemutsig_plot_check == "pdf"){
              pdf(file,width = input$treemutsig_width/100 , height = input$treemutsig_height/100)
          }
          if(class(treemutsig()) == "list"){
              print(treemutsig()[[getpatient.treemutsig()]])
          }else{
              print(treemutsig())
          }
          dev.off()
      },
      contentType = paste('image/',input$Download_treemutsig_plot_check,sep="")
  )
  
  output$treemutsig_table_ui <- renderUI({
    if(!is.null(treemutsig())){
      tagList(
        h4(strong('Signature summary')),
        br(),
        DT::dataTableOutput('treemutsig_table'),
        br(),
        fluidRow(
          column(
            width = 9
          ),
          column(
            width = 3,
            downloadBttn('Download_treemutsig_table', 'Download')
          )
        )
      )
    }
  })

  
  muttrunkbranch <- eventReactive(input$submit_muttrunkbranch, {
      maf <- isolate(varsMaf$maf)
      validate(
          need(!is.null(maf), "Upload maf file in the section 'Input Data'")
      )
      withProgress(min = 0, max = 2, value = 0,{
          setProgress(message = 'Generating ', detail = paste("phyloTree/phyloTreeList Class", sep="")) 
          phyloTree <- getPhyloTree(maf,
                                    method = input$muttrunkbranch_getphylotree_method,
                                    min.vaf = as.numeric(input$muttrunkbranch_getphylotree_minvaf),
                                    min.ccf = as.numeric(input$muttrunkbranch_getphylotree_minccf),
                                    bootstrap.rep.num = as.numeric(input$muttrunkbranch_getphylotree_bootstraprepnum))
          incProgress(amount = 1)
          setProgress(message = 'mutTrunkBranch : Calculation in progress',
                      detail = 'This may take a while...')
          
          mtb <- mutTrunkBranch(phyloTree, CT = input$muttrunkbranch_ct,
                                pvalue = as.numeric(input$muttrunkbranch_pvalue))
          incProgress(amount = 1)
          setProgress(message = 'mutTrunkBranch done!')
          return(mtb)
      })
  })
  
  getpatient.muttrunkbranch <- eventReactive(input$muttrunkbranch.pl,{
      return(input$muttrunkbranch.pl)
  })
  
  
  
  output$muttrunkbranch.patientlist <- renderUI({
      if(!is.null(muttrunkbranch())){
          if(!"mutTrunkBranch.plot" %in% names(muttrunkbranch()) ){
              x <- muttrunkbranch()
              patient.list <- names(x)
              selectInput("muttrunkbranch.pl", "Patient",
                          choices = patient.list, width = 600)
          }
      }
  })
  
  muttrunkbranch_width <- reactive({
      return(input$muttrunkbranch_width)
  })
  muttrunkbranch_height <- reactive({
      return(input$muttrunkbranch_height)
  })
  
  output$muttrunkbranch_plot <- renderPlot({
      if(!is.null(muttrunkbranch())){
          if(!"mutTrunkBranch.plot" %in% names(muttrunkbranch())){
              p <- muttrunkbranch()[[getpatient.muttrunkbranch()]]$mutTrunkBranch.plot
          }else{
              p <- muttrunkbranch()$mutTrunkBranch.plot
          }
          return(p) 
      }
  },
  width = muttrunkbranch_width,
  height = muttrunkbranch_height,
  res = 100
  )
  
  output$muttrunkbranch_download_button_ui <- renderUI({
      if(!is.null(muttrunkbranch())){
          fluidRow(
              column(
                  width = 7
              ),
              column(
                  width = 2,
                  radioButtons('Download_muttrunkbranch_plot_check', label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                               choiceNames = list(
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                                   tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                               ),choiceValues = c("png", "pdf"), inline = T)
              ),
              column(
                  width = 3,
                  downloadBttn('Download_muttrunkbranch_plot', 'Download')
              )
          )
      }
  })
  
  output$Download_muttrunkbranch_plot <- downloadHandler(
      filename = function() {
          paste("Rplot.",input$Download_muttrunkbranch_plot_check, sep='')
      },
      content = function(file) {
          if (input$Download_muttrunkbranch_plot_check == "png"){
              png(file,width = input$muttrunkbranch_width , height = input$muttrunkbranch_height,res = 100)
          }
          else if (input$Download_muttrunkbranch_plot_check == "pdf"){
              pdf(file,width = input$muttrunkbranch_width/100 , height = input$muttrunkbranch_height/100)
          }
          if(!"mutTrunkBranch.plot" %in% names(muttrunkbranch())){
              p <- muttrunkbranch()[[getpatient.muttrunkbranch()]]$mutTrunkBranch.plot
          }else{
              p <- muttrunkbranch()$mutTrunkBranch.plot
          }
          print(p)
          dev.off()
      },
      contentType = paste('image/',input$Download_muttrunkbranch_plot_check,sep="")
  )
  
  output$muttrunkbranch_table_ui <- renderUI({
    if(!is.null(muttrunkbranch())){
        tagList(
            h4(strong('Summary')),
            br(),
            DT::dataTableOutput('muttrunkbranch_table'),
            fluidRow(
                column(
                    width = 9
                ),
                column(
                    width = 3,
                    downloadBttn('Download_muttrunkbranch_table', 'Download')
                )
            )
        ) 
    }
  })
  
  output$muttrunkbranch_table <- DT::renderDataTable({
    if(!is.null(muttrunkbranch())){
        if(!"mutTrunkBranch.plot" %in% names(muttrunkbranch())){
            t <- muttrunkbranch()[[getpatient.muttrunkbranch()]]$mutTrunkBranch.res
        }else{
            t <- muttrunkbranch()$mutTrunkBranch.res
        }
        dt <- datatable(t, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5) 
        return(dt)
    }
  })
  
  
  output$Download_muttrunkbranch_table <- downloadHandler(
      filename = "Rtable.csv",
      content = function(file){
          if(!"mutTrunkBranch.plot" %in% names(muttrunkbranch())){
              t <- muttrunkbranch()[[getpatient.muttrunkbranch()]]$mutTrunkBranch.res
          }else{
              t <- muttrunkbranch()$mutTrunkBranch.res
          }
          write.csv(t,file,row.names = F)
      },
      contentType = 'text/csv'
  )

  
  
  
  
   
})  
  
  ## Download control  