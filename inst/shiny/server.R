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
    name <- input$maf$name
    patientID <- strsplit(name,"\\.")[[1]][1]
    return(patientID)
  })
  

  
  inputSilent <- observe({
    if (!is.null(input$maf)) {
      mafFile <- input$maf$datapath
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
      mafFile <- system.file("extdata/maf", "HCC6046.maf", package = "MesKit")
      sampleInfoFile <- system.file("extdata", "HCC6046.sampleInfo.txt", package = "MesKit")
      ccfClusterTsvFile <- system.file("extdata/ccf", "HCC6046.cluster.tsv", package = "MesKit")
      ccfLociTsvFile <- system.file("extdata/ccf", "HCC6046.loci.tsv", package = "MesKit")
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
      
      if (!is.null(input$chrSilent)){
        ls.chrSilent <- strsplit(input$chrSilent, ",")
      } else {
        ls.chrSilent <- NULL
      }
      
      if(is.null(input$maf) | is.null(input$sampleInfo)){
        mafFile <- system.file("extdata/maf", "HCC6046.maf", package = "MesKit")
        sampleInfoFile <- system.file("extdata", "HCC6046.sampleInfo.txt", package = "MesKit")
        ccfClusterTsvFile <- system.file("extdata/ccf", "HCC6046.cluster.tsv", package = "MesKit")
        ccfLociTsvFile <- system.file("extdata/ccf", "HCC6046.loci.tsv", package = "MesKit")
        maf <- readMaf(mafFile = mafFile, 
                       sampleInfoFile = sampleInfoFile,
                       mutType=input$mutType, 
                       mutNonSilent=ls.mutNonSilent, 
                       chrSilent=ls.chrSilent, 
                       use.indel = input$useindel, 
                       ccfClusterTsvFile = ccfClusterTsvFile,
                       ccfLociTsvFile = ccfLociTsvFile, 
                       refBuild="hg19")
      } else {
        if(!is.null(input$ccf.cluster)&!is.null(input$ccf.loci)){
          name <- mafName()
          maf <- MesKit::readMaf(mafFile = input$maf$datapath,
                                 sampleInfoFile = input$sampleInfo$datapath,
                                 ccfClusterTsvFile =  input$ccf.cluster$datapath,
                                 ccfLociTsvFile = input$ccf.loci$datapath,name = name)
        }
        else{
          name <- mafName()
          maf <-  MesKit::readMaf(mafFile = input$maf$datapath,
                         sampleInfoFile = input$sampleInfo$datapath,name = name)
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
      setProgress(message = 'Input data: Generating ', detail = paste("NJtree from MAF ", isolate(varsLs$maf)@patientID, sep="")) 
      varsLs[['njtree']] <-  njtree <- MesKit::getNJtree(isolate(varsLs$maf))
      incProgress(amount=1)
      
      setProgress(message = paste("Input data: MAF and NJtree Generation for ", isolate(varsLs$maf)@patientID, " Done!", sep=""), detail = "") 
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
            p("Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2, VAF, Tumor_Sample_Barcode.",
              style = "font-size:16px; font-weight:500;line-height:30px;"),
            h4(strong("Example MAF file")),
            style = "width: 800px"
          ),
          DT::dataTableOutput("ied1",width = "100%"),
          br()
        )
      )
    }
  })
  output$ied1 <- renderDataTable({
    if(input$iecontrol01){
      maftable <- read.table('dom/maf1.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(maftable, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
    }
  })
  ## output Introduction of sampleinfo datatable
  output$ie2 <- renderUI({
    if(buttonValue$b == 1){
      box(
        width = NULL,
        tagList(
          div(
            h3(strong("Information of samples")),
            p("Below is an example of the first four rows of sample_info.txt. It should contain the sampleID, patientID, lesion and sampling time. The input files are located under the '/inst/extdata/' folder.",
              style = "font-size:16px; font-weight:500;line-height:30px;"),
            style = "width : 800px"
          ),
          tags$li(strong("tumors sampling across multiple spatially-distinct regions"),
                  style = "width : 800px;font-size:16px; font-weight:500;"),
          br(),
          DT::dataTableOutput("ied2_1",width = "70%"),
          br(),
          tags$li(strong("tumors sampling across multiple time points"),
                  style = "width : 800px;font-size:16px; font-weight:500;"),
          br(),
          DT::dataTableOutput("ied2_2",width = "70%"),
          br(),
          div(tags$span("Note:",style = "font-size:16px; font-weight:700;"),
              tags$span("'-' represents sampling at the same time or sampling from the same site.",
                        style = "font-size:16px; font-weight:500;")),
          br()
        )
      )
    }
  })
  output$ied2_1 <- renderDataTable({
    if(input$iecontrol02){
      spd1 <- read.table('dom/sampleinfo1.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(spd1, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
    }
  })
  output$ied2_2 <- renderDataTable({
    if(input$iecontrol02){
      spd1 <- read.table('dom/sampleinfo2.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(spd1, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
    }
  })
  ## output Introduction of ccf.cluster
  output$ie3 <- renderUI({
    if(buttonValue$c == 1){
      box(
        width = NULL,
        tagList(
          h3(strong("Example ccf.cluster file")),
          DT::dataTableOutput("ied3"),
          br()
        )
      )
    }
  })
  output$ied3 <- renderDataTable({
    if(input$iecontrol03){
      spd3 <- read.table('dom/ccf.cluster.CSV',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(spd3, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
    }
  })
  ## output Introduction of ccf.loci
  output$ie4 <- renderUI({
    if(buttonValue$d == 1){
      box(
        width = NULL,
        tagList(
          h3(strong("Example ccf.loci file")),
          DT::dataTableOutput("ied4"),
          br()
        )
      )
    }
  })
  output$ied4 <- renderDataTable({
    if(input$iecontrol04){
      spd4 <- read.table('dom/ccf.loci.CSV',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(spd4, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
    }
  })
  output$maftable <- DT::renderDataTable({
    datatable(inputData()@data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T, fixedColumns = TRUE, columnDefs=list(list(width="10em",targets="_all"))),rownames = FALSE, width=5)
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
      mathScore(maf,tsb = c("All"),
                      minvaf = input$minvaf, 
                      maxvaf = input$maxvaf)$sampleLevel
    }
  })
  output$mathScore <- DT::renderDataTable({
    ms()
  })
  ms2 <- eventReactive(input$submit0, {
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'TMB: Calculation in progress',
                 detail = 'This may take a while...')
    
    if(input$submit0){
      maf <- isolate(varsLs$maf)
      getTMB(maf,tsb = c("All"),
                        minvaf = input$minvaf, 
                        maxvaf = input$maxvaf)$sampleLevel
    }
  })
  output$mathScoreTMB <- DT::renderDataTable({
    ms2()
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
  
  output$msdbtmb <- renderUI({
    if(!is.null(ms2())){
      fluidRow(
        column(
          width = 9
        ),
        column(
          width = 3,
          downloadBttn('DownloadTMB', 'Download')
        )
      )
    }
  })
  
  stopButtonValue3 <- reactiveValues(a = 0)
  observeEvent(input$stop3,{
    stopButtonValue3$a <- 1
  })
  observeEvent(input$submit3,{
    stopButtonValue3$a <- 0
  })
  vc <- eventReactive(input$submit3, {
    if(input$submit3 & stopButtonValue3$a != 1){
      maf <- isolate(varsLs$maf)
      tsbmax <- length(unique(maf@data$Tumor_Sample_Barcode))
      picSep <- NULL
      
      withProgress(min = 0, max = tsbmax+1, value = 0, {
        setProgress(message = 'vafCluster: Calculation in progress',
                    detail = 'This may take a while...')
        if (input$plotOption == "separate") {
          picSep <- vafClusterRshiny(maf,
                                     plotOption = input$plotOption, 
                                     themeOption = input$themeOption,
                                     showMATH = input$showMATH)
          
          output$chooselistvaf <- renderUI({
            if(!is.null(picSep)){
              names <- names(picSep)
              selectInput("vsl", "Branch",
                          choices = names, width = 600)
            }
          })
          
          output$vaf <- renderPlot({
            print(picSep[[which(names(picSep) == getOption())]])
          }, 
          width = width1,
          height = 560,
          res = 100
          )
          
        } else {
          pic <- vafClusterRshiny(maf,
                                  plotOption = input$plotOption, 
                                  themeOption = input$themeOption,
                                  showMATH = input$showMATH)
          
          output$vaf <- renderPlot({
            print(pic)
          }, 
          width = width1,
          height = 560,
          res = 100
          )
        }

        })
    }
  })
  

  
  getOption <- eventReactive(input$vsl, {
    if (input$plotOption == "separate"){
      return(input$vsl)
    }
  })
  
  
  
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
  stopButtonValue4 <- reactiveValues(a = 0)
  observeEvent(input$stop4,{
    stopButtonValue4$a <- 1
  })
  observeEvent(input$submit4,{
    stopButtonValue4$a <- 0
  })
  msp <- eventReactive(input$submit4, {
    if(input$submit4 & stopButtonValue4$a != 1){
      progress <- Progress$new(session, min=1, max=30)
      on.exit(progress$close())
      progress$set(message = 'Mutsharedprivateplot: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:30) {
        progress$set(value = i)
        Sys.sleep(2)
      }
      maf <- isolate(varsLs$maf)
      return(MesKit::mutSharedPrivate(maf,show.num = input$show.num1))
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
      mutOncoTSG(maf, oncogeneListFile = oncogeneListFile,
                   tsgListFile = tsgListFile, 
                   show.percentage = input$show.percentage)
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
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        validate(
          need(!(is.null(input$ccf.cluster$datapath)), "click the button 'use ccf',Upload ccf.cluster in Session 'Input Data' ")
        )
        validate(
          need(!(is.null(input$ccf.loci$datapath)), "Upload ccf.loci in Session 'Input Data'")
        )
        maf <- isolate(varsLs$maf)
        tumorClonesPlot(maf)
      }
      else{
        maf <- isolate(varsLs$maf)
        tumorClonesPlot(maf)
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
      njtree <- isolate(varsLs$njtree)
      MesKit::GO.njtree(njtree, 
                        GO.type = input$GO.type, 
                        plotType = input$plotType, 
                        pAdjustMethod=input$pAdjustMethod, 
                        qval = as.numeric(input$qval1), 
                        pval = as.numeric(input$pval1), 
                        showCategory = input$showCategory)
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
    datatable(data,options = list(pageLength = 5, dom = 'tp', scrollX = TRUE, columnDefs=list(list(width="10em",targets="_all"))), rownames = FALSE, width = 5)
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
      njtree <- isolate(varsLs$njtree)
      list <- MesKit::Pathway.njtree(njtree, 
                                     pathway.type=input$pathway.type, 
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
    return(DT::datatable(data,options = list(pageLength = 5, dom = 'tp', scrollX = T), rownames = FALSE,width = 5)) 
  })
  output$Pathwayplot <- renderPlot({
    return(Path()[[2]][[which(names(Path()[[2]]) == input$pl)]]) 
  },
  width = width7,
  height = height7,
  res = 100
  )
  stopButtonValueSig <- reactiveValues(a = 0)
  observeEvent(input$stopSig,{
    stopButtonValueSig$a <- 1
  })
  observeEvent(input$submitSig,{
    stopButtonValueSig$a <- 0
  })
  sigOFA <- eventReactive(input$submitSig, {
    if (input$oncogeneMapping) {
      if(is.null(input$driverGenesFile$datapath)){
        driverGenesFile <- system.file("extdata", "putative_driver_genes.txt", package = "MesKit")
      } else{
        driverGenesFile <- input$driverGenesFile$datapath
      }
      
      if(input$submitSig & stopButtonValueSig$a != 1){
        progress <- Progress$new(session, min=1, max=15)
        on.exit(progress$close())
        progress$set(message = 'Signature data summary: Calculation in progress',
                     detail = 'This may take a while...')
        
        for (i in 1:15) {
          progress$set(value = i)
          Sys.sleep(0.01)
        }
        
        njtree <- isolate(varsLs$njtree)
        df.signature <- treeMutationalSig(njtree, 
                                          driverGenesFile=driverGenesFile, 
                                          mutThreshold=input$mutThreshold, 
                                          signaturesRef=input$signaturesRef,
                                          plot.signatures=FALSE, 
                                          plot.branchTrunk=FALSE, 
                                          signif.level=0.05)
        return(df.signature)
      }
    } else {
      if(input$submitSig & stopButtonValueSig$a != 1){
        progress <- Progress$new(session, min=1, max=15)
        on.exit(progress$close())
        progress$set(message = 'Signature data summary: Calculation in progress',
                     detail = 'This may take a while...')
        
        for (i in 1:15) {
          progress$set(value = i)
          Sys.sleep(0.01)
        }
        
        njtree <- isolate(varsLs$njtree)
        df.signature <- treeMutationalSig(njtree, 
                                          driverGenesFile=NULL, 
                                          mutThreshold=input$mutThreshold, 
                                          signaturesRef=input$signaturesRef,
                                          plot.signatures=FALSE, 
                                          plot.branchTrunk=FALSE, 
                                          signif.level=0.05)
        return(df.signature)
        
      }
    }
  })
  output$sigOFAt <- DT::renderDataTable({
    return(datatable(sigOFA(), 
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
    if(input$submitSig1 & stopButtonValueSig1$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'Signature Plot: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      njtree <- isolate(varsLs$njtree)
      df.signature <- MesKit::treeMutationalSig(njtree, driverGenesFile=input$driverGenesFile$datapath, mutThreshold=input$mutThreshold, 
                                                signaturesRef=input$signaturesRef,
                                                plot.signatures=FALSE, plot.branchTrunk=FALSE, 
                                                signif.level=0.05)
      df.signature.plot <- MesKit::treeMutationalSig(njtree,
                                                     driverGenesFile=input$driverGenesFile1$datapath,
                                                     mutThreshold=input$mutThreshold1, 
                                                     signaturesRef=input$signaturesRef1,
                                                     plot.signatures=TRUE, plot.branchTrunk=FALSE, 
                                                     signif.level=0.05)
      return(list(df.signature.plot, df.signature))
      # return(datatable(df.signature, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T)))
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
    return(sigOFA1()[[1]])
  },
  width = widthsig1,
  height = heightsig1,
  res = 100
  )
  output$sigOFATable1 <- DT::renderDataTable({
    data <- sigOFA1()[[2]][,c(1:2)]
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
    if(input$submitSig2 & stopButtonValueSig2$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'TrunkOrBranch summary: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      njtree <- isolate(varsLs$njtree)
      df.branchTrunk.plot <- MesKit::treeMutationalSig(njtree, driverGenesFile=input$driverGenesFile2$datapath,
                                                       mutThreshold=input$mutThreshold2, 
                                                       signaturesRef=input$signaturesRef2,
                                                       plot.signatures=FALSE, plot.branchTrunk=TRUE, 
                                                       signif.level=input$signiflevel)
      return(df.branchTrunk.plot)
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
    return(sigOFA2()) 
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
    if(input$submit10 &stopButtonValue10$a != 1){
      progress <- Progress$new(session, min=1, max=15)
      on.exit(progress$close())
      progress$set(message = 'PhyloTree: Calculation in progress',
                   detail = 'This may take a while...')
      
      for (i in 1:15) {
        progress$set(value = i)
        Sys.sleep(0.01)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        njtree <- isolate(varsLs$njtree)
        if(input$useccf == T){
          validate(
            need(input$heatmap.type == "CCF","switch heatmap type to CCF")
          )
        }
        if(njtree@patientID == "0"){
          id <- input$maf$name
          njtree@patientID <- strsplit(id,"\\.")[[1]][1]
        }
        p <- plotPhyloTree(njtree, heatmap.type = input$heatmap.type, sig.name = "default",
                                   show.mutSig = input$showmutSig, show.heatmap = input$showheatmap)
        return(p)
        # else{
        #   validate(
        #     need(!is.null(input$phylotree.dir),"Upload your phylotree file")
        #   )
        #   p <- MesKit::plotPhyloTree(phylotree.dat = input$phylotree.dir$datapath, 
        #                              phylotree.type = input$phyloTreeType)
        #   return(p)
        # }
      }
      else{
        njtree <- isolate(varsLs$njtree)
        if(njtree@patientID == "0"){
          id <- input$maf$name
          njtree@patientID <- strsplit(id,"\\.")[[1]][1]
        }
        p <- plotPhyloTree(njtree, heatmap.type = input$heatmap.type, sig.name = "default",
                                   show.mutSig = input$showmutSig, show.heatmap = input$showheatmap)
        return(p)
        # inputData()$phylotreeplot
      }
    }
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
    filename = function() {
      paste("MathScore_",Sys.Date(),".csv", sep = '')
    },
    content = function(file){
      data <- ms()
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadTMB <- downloadHandler(
    filename = function() {
      paste("TMB_",Sys.Date(),".csv", sep = '')
    },
    content = function(file){
      data <- ms2()
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadVafPlot <- downloadHandler(
    filename = function() {
      paste("VafPlot",'.',input$DownloadVafPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadVafPlotCheck == "png"){
        png(file,width = input$width1 , height = 560,res = 100)
      }
      else if (input$DownloadVafPlotCheck == "pdf"){
        pdf(file,width = input$width1/100 , height = 6)
      }
      vc()
      dev.off()
    },
    contentType = paste('image/',input$DownloadVafPlotCheck,sep="")
  )
  
  output$DownloadStackPlot <- downloadHandler(
    filename = function() {
      paste("mutOncoTSG",'.',input$DownloadStackPlotCheck, sep='')
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
      paste("JaccardIndex",'.',input$DownloadJaccardIndexCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadJaccardIndexCheck == "png"){
        png(file,width = input$width4 , height = input$width4,res = 100)
      }
      else if (input$DownloadJaccardIndexCheck == "pdf"){
        pdf(file,width = input$width4/100 , height = input$width4/100)
      }
      ji()
      dev.off()
    },
    contentType = paste('image/',input$DownloadJaccardIndexCheck,sep="")
  )
  output$DownloadSharedPlot <- downloadHandler(
    filename = function() {
      paste("mutPrivateSharedPlot",'.',input$DownloadSharedPlotCheck, sep='')
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
      paste("ClonePlot",'.',input$DownloadClonePlotCheck, sep='')
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
      paste("ClonePlot",'.',input$DownloadClonePlotCheck, sep='')
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
      paste("PhyloTree",'.',input$DownloadPhyloTreeCheck, sep='')
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
      paste("GOPlot","_",input$gl,".",input$DownloadGOPlotCheck, sep='')
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
    filename = function() {
      paste("GO_",input$gl,"_",Sys.Date(),'.csv', sep='')
    },
    content = function(file){
      data <- GO()[[1]][[which(names(GO()[[1]]) == input$gl)]]
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  output$DownloadPathPlot <- downloadHandler(
    filename = function() {
      paste("Pathwayplot","_",input$pl,".",input$DownloadPathPlotCheck, sep='')
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
    filename = function() {
      paste("Pathway_",input$pl,"_",Sys.Date(),'.csv', sep='')
    },
    content = function(file){
      data <- Path()[[1]][[which(names(Path()[[1]]) == input$pl)]]
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  output$DownloadSignatureSummary <- downloadHandler(
    filename = function() {
      paste("Signature_summary_",Sys.Date(),'.csv', sep='')
    },
    content = function(file){
      data <- sigOFA()
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  output$DownloadSigOFATable1 <- downloadHandler(
    filename = function() {
      paste("Signature_summary_",Sys.Date(),'.csv', sep='')
    },
    content = function(file){
      data <- sigOFA1()[[2]][,c(1:2)]
      write.csv(data,file,row.names = F)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadSignaturePlot1 <- downloadHandler(
    filename = function() {
      paste("SignaturePlot",'.',input$DownloadSignaturePlotCheck1, sep='')
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
      paste("Branch_trunck",'.',input$DownloadSignaturePlotCheck2, sep='')
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
  
})  