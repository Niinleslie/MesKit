suppressMessages(library(shiny))
suppressMessages(library(Meskit))
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
    return(input$widthsig2)
  })
  heightsig1 <- reactive({
    return(input$heightsig2)
  })
  widthsig2 <- reactive({
    return(input$widthsig2)
  })
  heightsig2 <- reactive({
    return(input$heightsig2)
  })

  inputData <- eventReactive(input$submit1,{
    if(is.null(input$maf) | is.null(input$sampleInfo)){
      mafFile <- './example/311252.maf'
      sampleInfoFile <- './example/sample_info.txt'
      ccfClusterTsvFile <- './example/311252.cluster.tsv'
      ccfLociTsvFile <- './example/311252.loci.tsv'
      maf <- Meskit::readMaf(mafFile = mafFile,
                             sampleInfoFile = sampleInfoFile,
                             ccfClusterTsvFile =  ccfClusterTsvFile,
                             ccfLociTsvFile = ccfLociTsvFile)
    }
    else{
      if(!is.null(input$ccf.cluster)&!is.null(input$ccf.loci)){
        maf <- Meskit::readMaf(mafFile = input$maf$datapath,
                               sampleInfoFile = input$sampleInfo$datapath,
                               ccfClusterTsvFile =  input$ccf.cluster$datapath,
                               ccfLociTsvFile = input$ccf.loci$datapath,
                               inputFileName = input$maf$name)
      }
      else{
        maf <- readMaf(mafFile = input$maf$datapath,
                       sampleInfoFile = input$sampleInfo$datapath,
                       inputFileName = input$maf$name)
      }
    }
  })
  # inputData <- reactive({
  #   if(input$submit1){
  #     if(is.null(input$maf) | is.null(input$sampleInfo)){
  #       mafFile <- './example/311252.maf'
  #       sampleInfoFile <- './example/sample_info.txt'
  #       ccfClusterTsvFile <- './example/311252.cluster.tsv'
  #       ccfLociTsvFile <- './example/311252.loci.tsv'
  #       maf <- Meskit::readMaf(mafFile = mafFile,
  #                              sampleInfoFile = sampleInfoFile, 
  #                              ccfClusterTsvFile =  ccfClusterTsvFile, 
  #                              ccfLociTsvFile = ccfLociTsvFile)
  #     }
  #     else{
  #       if(!is.null(input$ccf.cluster)&!is.null(input$ccf.loci)){
  #         maf <- Meskit::readMaf(mafFile = input$maf$datapath,
  #                                sampleInfoFile = input$sampleInfo$datapath, 
  #                                ccfClusterTsvFile =  input$ccf.cluster$datapath, 
  #                                ccfLociTsvFile = input$ccf.loci$datapath,
  #                                inputFileName = input$maf$name)
  #       }
  #       else{
  #         maf <- readMaf(mafFile = input$maf$datapath, 
  #                        sampleInfoFile = input$sampleInfo$datapath,
  #                        inputFileName = input$maf$name)
  #       }
  #     }
  #   }
  #   if(is.null(input$maf) | is.null(input$sampleInfo)){
  #     mafFile <- './example/311252.maf'
  #     sampleInfoFile <- './example/sample_info.txt'
  #     ccfClusterTsvFile <- './example/311252.cluster.tsv'
  #     ccfLociTsvFile <- './example/311252.loci.tsv'
  #     maf <- Meskit::readMaf(mafFile = mafFile,
  #                            sampleInfoFile = sampleInfoFile, 
  #                            ccfClusterTsvFile =  ccfClusterTsvFile, 
  #                            ccfLociTsvFile = ccfLociTsvFile)
  #   }
  #   else{
  #     if(!is.null(input$ccf.cluster)&!is.null(input$ccf.loci)){
  #       maf <- Meskit::readMaf(mafFile = input$maf$datapath,
  #                              sampleInfoFile = input$sampleInfo$datapath, 
  #                              ccfClusterTsvFile =  input$ccf.cluster$datapath, 
  #                              ccfLociTsvFile = input$ccf.loci$datapath,
  #                              inputFileName = input$maf$name)
  #     }
  #     else{
  #       maf <- readMaf(mafFile = input$maf$datapath, 
  #                      sampleInfoFile = input$sampleInfo$datapath,
  #                      inputFileName = input$maf$name)
  #     }
  #   }
  # })
  
  inputNJtree <- reactive({
    maf <- inputData()
    njtree <- Meskit::getNJtree(maf, use.indel = input$use.indel)
  })
  ## output Introduction of maf datatable
  output$ie1 <- renderUI({
    if((input$iecontrol01)%%2 != 0 ){
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
          DT::dataTableOutput("ied1",width = "70%"),
          br()
        )
      )
    }
  })
  output$ied1 <- renderDataTable({
    if(input$iecontrol01){
      md1 <- read.table('dom/maf1.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
      return(DT::datatable(md1, options = list(pageLength = 5, dom = 't', scrollX = T),rownames = FALSE)) 
    }
  })
  ## output Introduction of sampleinfo datatable
  output$ie2 <- renderUI({
    if((input$iecontrol02)%%2 != 0){
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
      datatable(spd1, options = list(pageLength = 5, dom = 't', scrollX = T), rownames = FALSE,width = 5)
    }
  })
  output$ied2_2 <- renderDataTable({
    if(input$iecontrol02){
      spd1 <- read.table('dom/sampleinfo2.csv',encoding = "UTF-8",sep = ",",header = T,fill = T)
      datatable(spd1, options = list(pageLength = 5, dom = 't', scrollX = T), rownames = FALSE,width = 5)
    }
  })
  ## output Introduction of ccf.cluster
  output$ie3 <- renderUI({
    if((input$iecontrol03)%%2 != 0){
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
      datatable(spd3, options = list(pageLength = 6, dom = 't', scrollX = T), rownames = FALSE,width = 5)
    }
  })
  ## output Introduction of ccf.loci
  output$ie4 <- renderUI({
    if((input$iecontrol04)%%2 != 0){
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
      datatable(spd4, options = list(pageLength = 6, dom = 't', scrollX = T), rownames = FALSE,width = 5)
    }
  })
  output$maftable <- DT::renderDataTable({
    datatable(inputData()@data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T))
  })

  ms <- eventReactive(input$submit2,{
    maf <- inputData()
    Meskit::mathScore(maf,tsb = c("All"),
                      minvaf = input$minvaf,maxvaf = input$maxvaf)$sampleLevel
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
  
  vc <- eventReactive(input$submit3,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    maf <- inputData()
    Meskit::vafCluster(maf,plotOption = input$plotOption,themeOption = input$themeOption)
  })
  output$chooselistvaf <- renderUI({
    names <- names(vc())
    selectInput("vsl","Branch",
                choices = names ,selected = names[1],width = 600)
  })
  output$vaf <- renderPlot({
    if(input$plotOption == "separate"){
      return(vc()[[which(names(vc()) == input$vsl)]])
      print(names(vc()))
    }
    else{
      vc()
    }
  }, 
  width = width1,
  height = 560,
  res = 100
  )
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
  msp <- eventReactive(input$submit4,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    maf <- inputData()
    return(Meskit::mutSharedPrivate(maf,show.num = input$show.num1))
  })
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
  })
  stk <- eventReactive(input$submit5,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    if(is.null(input$oncogeneListFile$datapath)){
      oncogeneListFile <- './example/oncogene.list.txt'
    }
    else{
      oncogeneListFile <- input$oncogeneListFile$datapath
    }
    if(is.null(input$tsgListFile$datapath)){
      tsgListFile <- './example/TSG.list.txt'
    }
    else{
      tsgListFile <- input$tsgListFile$datapath
    }
    maf <- inputData()
    Meskit::mutStackPlot(maf, oncogeneListFile = oncogeneListFile,
                         tsgListFile = tsgListFile, themeOption=input$themeOption2,
                         show.percentage = input$show.percentage)
  })
  output$stackplot <- renderPlot({
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
          width = 7
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
  ji <- eventReactive(input$submit6,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    maf <- inputData()
    return(Meskit::JaccardIndex(maf,type = input$JItype))
  })
  output$JaccardIndex <- renderPlot({
    ji()
  },
  width = width4,
  height = 560,
  res = 100
  )
  output$jidb <- renderUI({
    if(!is.null(ji())){
      fluidRow(
        column(
          width = 7
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
  clp <- eventReactive(input$submit7,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      validate(
        need(!(is.null(input$ccf.cluster$datapath)), "click the button 'use ccf',Upload ccf.cluster in Session 'Input Data' ")
      )
      validate(
        need(!(is.null(input$ccf.loci$datapath)), "Upload ccf.loci Session 'Input Data'")
      )
      maf <- inputData()
      Meskit::tumorClonesPlot(maf)
    }
    else{
      maf <- inputData()
      Meskit::tumorClonesPlot(maf)
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
          downloadBttn('DownloadClonePlot', 'Download')
        )
      )
    }
  })
  GO <- eventReactive(input$submit8,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    njtree <- inputNJtree()
    Meskit::GO.njtree(njtree, qval = as.numeric(input$qval1) ,pval = as.numeric(input$pval1))
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
            downloadBttn('DownloadGOTable', 'Download')
          )
        )
      )
    }
  })
  output$gotable <- renderDataTable({
    data <- GO()[[1]][[which(names(GO()[[1]]) == input$gl)]]
    datatable(data,options = list(pageLength = 5, dom = 'tp', scrollX = T), rownames = FALSE,width = 5)
  })
  output$chooselist1 <- renderUI({
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
          width = 7
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
  Path <- eventReactive(input$submit9,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    njtree <- inputNJtree()
    list <- Meskit::Pathway.njtree(njtree, qval = as.numeric(input$qval2) ,pval = as.numeric(input$pval2))
    return(list)
  })
  output$Pathdb <- renderUI({
    if(!is.null(Path())){
      fluidRow(
        column(
          width = 7
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
          downloadBttn('DownloadPathPlot', 'Download')
        )
      )
      
    }
  })
  output$chooselist2 <- renderUI({
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
  
  sigOFA1 <- eventReactive(input$submitSig1,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    njtree <- inputNJtree()
    df.signature <- Meskit::treeMutationalSig(njtree,
                                              driverGenesFile=input$driverGenesFile$datapath1,
                                              mutThreshold=input$mutThreshold1, 
                                              signaturesRef=input$signaturesRef1,
                                              plot.signatures=FALSE, plot.branchTrunk=FALSE, 
                                              signif.level=0.05)
    df.signature.plot <- Meskit::treeMutationalSig(njtree,
                                                   driverGenesFile=input$driverGenesFile1$datapath,
                                                   mutThreshold=input$mutThreshold1, 
                                                   signaturesRef=input$signaturesRef1,
                                                   plot.signatures=TRUE, plot.branchTrunk=FALSE, 
                                                   signif.level=0.05)
    return(list(df.signature.plot,df.signature))
    # return(datatable(df.signature, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T)))
  })
  output$sigOFATableUI1 <- renderUI({
    if(!is.null(sigOFA1()[[2]])){
      tagList(
        h4(strong('Signature summary')),
        br(),
        DT::dataTableOutput('sigOFATable1'),
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
    data <- sigOFA1()[[2]]
    datatable(data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T))
  })
  sigOFA2 <- eventReactive(input$submitSig2,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    njtree <- inputNJtree()
    df.branchTrunck <- Meskit::treeMutationalSig(njtree, driverGenesFile=input$driverGenesFile$datapath2,
                                                 mutThreshold=input$mutThreshold2, 
                                                 signaturesRef=input$signaturesRef2,
                                                 plot.signatures=FALSE, plot.branchTrunk=FALSE, 
                                                 signif.level=0.05)
    df.branchTrunk.plot <- Meskit::treeMutationalSig(njtree, driverGenesFile=input$driverGenesFile2$datapath,
                                                     mutThreshold=input$mutThreshold2, 
                                                     signaturesRef=input$signaturesRef2,
                                                     plot.signatures=FALSE, plot.branchTrunk=TRUE, 
                                                     signif.level=input$signiflevel)
    return(list(df.branchTrunk.plot,df.branchTrunck))
  })
  output$sigOFATableUI2 <- renderUI({
    if(!is.null(sigOFA2()[[2]]))
    tagList(
      h4(strong('Signature summary')),
      br(),
      DT::dataTableOutput('sigOFATable2'),
      fluidRow(
        column(
          width = 9
        ),
        column(
          width = 3,
          downloadBttn('DownloadSigOFATable2', 'Download')
        )
      )
    )
  })
  output$sigOFATable2 <- renderDataTable({
    if(!is.null(sigOFA2()[[2]])){
      return(datatable(sigOFA2()[[2]], options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T)))
    }
  })
  
  output$sigOFAPlot2 <- renderPlot({
    return(sigOFA2()[[1]]) 
  },
  width = widthsig2,
  height = heightsig2,
  res = 100
  )
  
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
  
  pht <- eventReactive(input$submit10,{
    progress <- Progress$new(session, min=1, max=15)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...')
    
    for (i in 1:15) {
      progress$set(value = i)
      Sys.sleep(0.1)
    }
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      if(input$phyloTreeType == 'njtree'){
        njtree <- inputNJtree()
        if(input$useccf == T){
          validate(
            need(input$heatmap.type == "CCF","switch heatmap type to CCF")
          )
        }
        p <- Meskit::plotPhyloTree(njtree, phylotree.type = input$phyloTreeType, 
                                   heatmap.type = input$heatmap.type, sig.name = "alias",
                                   show.mutSig = input$show.mutSig, show.heatmap = input$show.heatmap)
        return(p)
      }
      else{
        validate(
          need(!is.null(input$phylotree.dir),"Upload your phylotree file")
        )
        p <- Meskit::plotPhyloTree(phylotree.dat = input$phylotree.dir$datapath, 
                                   phylotree.type = input$phyloTreeType)
        return(p)
      }
    }
    else{
      njtree <- inputNJtree()
      p <- Meskit::plotPhyloTree(njtree, phylotree.type = input$phyloTreeType, 
                                 heatmap.type = input$heatmap.type, sig.name = "alias",
                                 show.mutSig = input$show.mutSig, show.heatmap = input$show.heatmap)
      return(p)
      # inputData()$phylotreeplot
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
      write.csv(data,file)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadVafPlot <- downloadHandler(
    filename = function() {
      paste("VafPlot",'.',input$DownloadVafPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadVafPlotCheck == "png"){
        png(file,width = input$width1 , height = 900,res = 144)
      }
      else if (input$DownloadVafPlotCheck == "pdf"){
        pdf(file,width = input$width1/100 , height = 9)
      }
      print(vc())
      dev.off()
    },
    contentType = paste('image/',input$DownloadVafPlotCheck,sep="")
  )
  output$DownloadStackPlot <- downloadHandler(
    filename = function() {
      paste("StackPlot",'.',input$DownloadStackPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadStackPlotCheck == "png"){
        png(file,width = input$width3 , height = 900,res = 144)
      }
      else if (input$DownloadStackPlotCheck == "pdf"){
        pdf(file,width = input$width3/100 , height = 9)
      }
      print(stk())
      dev.off()
    },
    contentType = paste('image/',input$DownloadVafPlotCheck,sep="")
  )
  output$DownloadJaccardIndex <- downloadHandler(
    filename = function() {
      paste("JaccardIndex",'.',input$DownloadJaccardIndexCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadJaccardIndexCheck == "png"){
        png(file,width = input$width4 , height = 900,res = 144)
      }
      else if (input$DownloadStackPlotCheck == "pdf"){
        pdf(file,width = input$width4/100 , height = 9)
      }
      print(ji())
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
        png(file,width = input$width2 , height = 900, res = 144)
        
      }
      else if (input$DownloadSharedPlotCheck == "pdf"){
        pdf(file,width = input$width2/100 , height = 9)
      }
      print(msp())
      dev.off()
    },
    contentType = paste('image/',input$DownloadSharedPlotCheck,sep="")
  )
  
  output$DownloadClonePlot <- downloadHandler(
    filename = function() {
      paste("ClonePlot",'.',input$DownloadClonePlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadClonePlotCheck == "png"){
        png(file,width = input$width5 , height = 800,res = 144)
      }
      else if (input$DownloadClonePlotCheck == "pdf"){
        pdf(file,width = input$width5/100 , height = 8)
      }
      clp()
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
        png(file,width = 1000, height = 650,res = 80)
      }
      else if (input$DownloadPhyloTreeCheck == "pdf"){
        pdf(file,width = 10, height = 6.5)
      }
      print(pht)
      dev.off()
    }
  )
  
  output$DownloadGOPlot <- downloadHandler(
    filename = function() {
      paste("GOPlot", '.',input$DownloadGOPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadGOPlotCheck == "png"){
        png(file,width = input$width6, height = input$height6,res = 144)
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
      write.csv(data,file)
    },
    contentType = 'text/csv'
  )
  output$DownloadPathPlot <- downloadHandler(
    filename = function() {
      paste("Pathwaytplot",'.',input$DownloadPathPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadPathPlotCheck == "png"){
        png(file,width = input$width7,height = input$height7,res = 144)
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
      write.csv(data,file)
    },
    contentType = 'text/csv'
  )
  
  output$DownloadSignaturePlot <- downloadHandler(
    filename = function() {
      paste("SignaturePlot",'.',input$DownloadSignaturePlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadSignaturePlotCheck == "png"){
        png(file,width = input$widthsig2, height = input$heightsig2,res = 144)
      }
      else if (input$DownloadSignaturePlotCheck == "pdf"){
        pdf(file,width = input$widthsig2/100, height = input$heightsig2/100)
      }
      print(sigOFA2())
      dev.off()
    },
    contentType = paste('image/',input$DownloadSignaturePlotCheck,sep="")
  )
  
})