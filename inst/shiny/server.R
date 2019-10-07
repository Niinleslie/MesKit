suppressMessages(library(shiny))



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
  
  inputData <- eventReactive(input$submit1,{
    if(is.null(input$maf) | is.null(input$sampleInfo)){
      maf <- readRDS("./example/MesKit_Example.rds")
    }
    else{
      if(!is.null(input$ccf.cluster)&!is.null(input$ccf.loci)){
        maf <- Meskit::readMaf(mafFile = input$maf$datapath,
                               sampleInfoFile = input$sampleInfo$datapath, 
                               ccfClusterTsvFile =  input$ccf.cluster$datapath, 
                               ccfLociTsvFile = input$ccf.loci$datapath)
      }
      else{
        maf <- Meskit::readMaf(mafFile = input$maf$datapath, 
                               sampleInfoFile = input$sampleInfo$datapath)
      }
    }
  })
  
  inputNJtree <- reactive({
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      maf <- inputData()
      njtree <- Meskit::getNJtree(maf, use.indel = input$use.indel)
    }
    else{
      maf <- inputData()$maf
      njtree <- Meskit::getNJtree(maf, use.indel = input$use.indel)
    }
  })
  
  output$maftable <- DT::renderDataTable({
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      datatable(inputData()@data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T))
    } else {
      datatable(inputData()$maf@data, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T))
    }
  })
  
  ms <- eventReactive(input$submit2,{
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      maf <- inputData()
      Meskit::mathScore(maf,tsb = c("All"),
                        minvaf = input$minvaf,maxvaf = input$maxvaf)$sampleLevel
    }
    else{
      maf <- inputData()$maf
      Meskit::mathScore(maf,tsb = c("All"),
                        minvaf = input$minvaf,maxvaf = input$maxvaf)$sampleLevel
    }
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
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      maf <- inputData()
      Meskit::vafCluster(maf,plotOption = input$plotOption,themeOption = input$themeOption)
    }
    else{
      maf <- inputData()$maf
      Meskit::vafCluster(maf,plotOption = input$plotOption,themeOption = input$themeOption)
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
  output$vaf <- renderPlot({
    vc()
  }, 
  width = width1,
  height = 560,
  res = 100
  )
  msp <- eventReactive(input$submit4,{
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      maf <- inputData()
      Meskit::mutSharedPrivate(maf,show.num = input$show.num)
    }
    else{
      maf <- inputData()$maf
      Meskit::mutSharedPrivate(maf,show.num = input$show.num)
    }
  })
  output$mut.share_private <- renderPlot({
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
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      maf <- inputData()
      Meskit::mutStackPlot(maf, oncogeneListFile = oncogeneListFile,
                           tsgListFile = tsgListFile, themeOption=input$themeOption2, show.percentage = input$show.percentage)
    }
    else{
      maf <- inputData()$maf
      Meskit::mutStackPlot(maf, oncogeneListFile = oncogeneListFile,
                           tsgListFile = tsgListFile, themeOption=input$themeOption2, show.percentage = input$show.percentage)
    }
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
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      maf <- inputData()
      Meskit::JaccardIndex(maf,type = input$JItype)
    }
    else{
      maf <- inputData()$maf
      Meskit::JaccardIndex(maf,type = input$JItype)
    }
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
      maf <- inputData()$maf
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
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      njtree <- inputNJtree()
      Meskit::GO.njtree(njtree, qval = as.numeric(input$qval1) ,pval = as.numeric(input$pval1))
    }
    else{
      njtree <- inputNJtree()
      Meskit::GO.njtree(njtree, qval = as.numeric(input$qval1) ,pval = as.numeric(input$pval1))
    }
  })
  output$chooselist1 <- renderUI({
    selectInput("gl","Branch",
                choices = names(GO()[[2]]) ,selected = names(GO()[[2]])[1],width = 600)
  })
  output$GOplot <- renderPlot({
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      return(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
    }
    else{
      return(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
    }
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
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      njtree <- inputNJtree()
      list <- Meskit::Pathway.njtree(njtree, qval = as.numeric(input$qval2) ,pval = as.numeric(input$pval2))
      return(list)
    }
    else{
      
    }
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
    selectInput("pl","Branch",
                choices = names(Path()[[2]]) ,selected = names(Path()[[2]])[1],width = 600)
  })
  output$Pathwayplot <- renderPlot({
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      return(Path()[[2]][[which(names(Path()[[2]]) == input$pl)]])
    }
    else{
      return(Path()[[2]][[which(names(Path()[[2]]) == input$pl)]])
    }
  },
  width = width7,
  height = height7,
  res = 100
  )
  
  sigOFA <- eventReactive(input$submitSig,{
    njtree <- inputNJtree()
    df.signature <- Meskit::treeMutationalSig(njtree, driverGenesFile=input$driverGenesFile$datapath, mutThreshold=input$mutThreshold, 
                                              signaturesRef=input$signaturesRef,
                                              plot.signatures=FALSE, plot.branchTrunk=FALSE, 
                                              signif.level=0.05)
    return(datatable(df.signature, options = list(searching = TRUE, pageLength = 10, lengthMenu = c(5, 10, 15, 18), scrollX = T)))
  })
  
  output$sigOFA <- DT::renderDataTable({
    sigOFA()
  })
  
  sigOFA2 <- eventReactive(input$submitSig2,{
    njtree <- inputNJtree()
    if (input$sigplot == "signaturesProb") {
      Meskit::treeMutationalSig(njtree, driverGenesFile=input$driverGenesFile2$datapath, mutThreshold=input$mutThreshold2, 
                                signaturesRef=input$signaturesRef2,
                                plot.signatures=TRUE, plot.branchTrunk=FALSE, 
                                signif.level=0.05)
    } else if (input$sigplot == "branchTrunk"){
      Meskit::treeMutationalSig(njtree, driverGenesFile=input$driverGenesFile2$datapath, mutThreshold=input$mutThreshold2, 
                                signaturesRef=input$signaturesRef2,
                                plot.signatures=FALSE, plot.branchTrunk=TRUE, 
                                signif.level=input$signifLevel)
    }
  })
  
  ## TODO: how to extract the print pictures in R shiny.
  
  pht <- eventReactive(input$submit10,{
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      if(input$phyloTreeType == 'njtree'){
        njtree <- inputNJtree()
        if(input$useccf == T){
          validate(
            need(input$heatmap.type == "CCF","switch heatmap type to CCF")
          )
        }
        p <- Meskit::plotPhyloTree(njtree, phylotree.type = input$phyloTreeType, 
                                   heatmap.type = input$heatmap.type, sig.name = input$sig.name,
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
                                 heatmap.type = input$heatmap.type, sig.name = input$sig.name,
                                 show.mutSig = input$show.mutSig, show.heatmap = input$show.heatmap)
      return(p)
      # inputData()$phylotreeplot
    }
  })
  
  output$phylotree <- renderPlot({
    pht()
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
                       ),inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadPhyloTree', 'Download')
        )
      )
    }
  })
  sigp <- eventReactive(input$submit11,{
    if(!is.null(input$maf) & !is.null(input$sampleInfo)){
      if(input$phyloTreeType == 'njtree'){
        njtree <- inputNJtree()
        Meskit::treeMutationalSig(njtree,plot.Signatures = T)
      }
      else{
        message("Drawing signature plot when using njtree")
      }
    }
    else{
      njtree <- inputNJtree()
      Meskit::treeMutationalSig(njtree,plot.Signatures = T)
    }
  })
  output$signature <- renderPlot({
    sigp()
  },
  res = 100
  )
  output$sigpdb <- renderUI({
    if(!is.null(sigp())){
      fluidRow(
        column(
          width = 7
        ),
        column(
          width = 2,
          radioButtons('DownloadSignaturePlotCheck', label = div(style = "font-size:18px; font-weight: bold; ", 'Save type as:'),
                       choiceNames = list(
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "png"), 
                         tags$span(style = "font-size:14.5px; font-weight:400; ", "pdf")
                       ),inline = T)
        ),
        column(
          width = 3,
          downloadBttn('DownloadSignaturePlot', 'Download')
        )
      )
    }
  })
  output$DownloadMathScore <- downloadHandler(
    filename = function() {
      paste("MathScore",'.',"csv", sep='')
    },
    content = function(file) {
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        maf <- inputData()
        data <- Meskit::mathScore(maf = maf)$sampleLevel
      }
      else{
        maf <- inputData()$maf
        data <- Meskit::mathScore(maf = maf)$sampleLevel
      }
      
      write.csv(data,file)
    }
  )
  
  output$DownloadVafPlot <- downloadHandler(
    filename = function() {
      paste("VafPlot",'.',input$DownloadVafPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadVafPlotCheck == "png"){
        png(file,width = 1200 , height = 900,res = 144)
      }
      else if (input$DownloadVafPlotCheck == "pdf"){
        pdf(file,width = 12 , height = 9)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        maf <- inputData()
        Meskit::vafCluster(maf)
      }
      else{
        maf <- inputData()$maf
        Meskit::vafCluster(maf)
      }
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
        png(file,width = 1200 , height = 900,res = 144)
      }
      else if (input$DownloadStackPlotCheck == "pdf"){
        pdf(file,width = 12 , height = 9)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        maf <- inputData()
        validate(
          need(!((is.null(input$oncogeneListFile$datapath) & is.null(input$tsgListFile$datapath))), 
               "Upload oncogeneListFile and tsgListFile in 'Setting&Upload' ")
        )
        Meskit::mutStackPlot(maf, oncogeneListFile = input$oncogeneListFile$datapath,
                             tsgListFile = input$tsgListFile$datapath, themeOption="npg", show.percentage = TRUE)
      }
      else{
        maf <- inputData()$maf
        Meskit::mutStackPlot(maf, oncogeneListFile = oncogeneListFile,
                             tsgListFile = tsgListFile, themeOption=input$themeOption2, show.percentage = input$show.percentage)
      }
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
        png(file,width = 1200 , height = 900,res = 144)
      }
      else if (input$DownloadStackPlotCheck == "pdf"){
        pdf(file,width = 12 , height = 9)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        maf <- inputData()
        Meskit::JaccardIndex(maf)
      }
      else{
        maf <- inputData()$maf
        Meskit::JaccardIndex(maf)
      }
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
        png(file,width = 1300 , height = 1000, res = 144)
        
      }
      else if (input$DownloadSharedPlotCheck == "pdf"){
        pdf(file,width = 13 , height = 10)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        maf <- inputData()
        Meskit::mutSharedPrivate(maf)
      }
      else{
        maf <- inputData()$maf
        Meskit::mutSharedPrivate(maf) 
      }
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
        png(file,width = 1500 , height = 800,res = 144)
      }
      else if (input$DownloadClonePlotCheck == "pdf"){
        pdf(file,width = 15 , height = 8)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        maf <- inputData()
        Meskit::tumorClonesPlot(maf)
      }
      else{
        maf <- inputData()$maf
        Meskit::tumorClonesPlot(maf)
      }
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
        png(file,width = 1400, height = 800,res = 80)
      }
      else if (input$DownloadPhyloTreeCheck == "pdf"){
        ggsave(file,inputData()$phylotreeplot,width = 14, height = 8)
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
                                     heatmap.type = input$heatmap.type, sig.name = input$sig.name,
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
                                   heatmap.type = input$heatmap.type, sig.name = input$sig.name,
                                   show.mutSig = input$show.mutSig, show.heatmap = input$show.heatmap)
        return(p)
        # inputData()$phylotreeplot
      }
      dev.off()
    }
  )
  
  output$DownloadGOPlot <- downloadHandler(
    filename = function() {
      paste("GOPlot", '.',input$DownloadGOPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadGOPlotCheck == "png"){
        png(file,width = 2000, height = 1600,res = 144)
      }
      else if (input$DownloadGOPlotCheck == "pdf"){
        pdf(file,width = 20, height = 16)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        return(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
      }
      else{
        return(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
      }
      dev.off()
    },
    contentType = paste('image/',input$DownloadGOPlotCheck,sep="")
  )
  output$DownloadPathPlot <- downloadHandler(
    filename = function() {
      paste("pathwatplot",'.',input$DownloadGOPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadPathPlotCheck == "png"){
        png(file,width = 2000, height = 1600,res = 144)
      }
      else if (input$DownloadPathPlotCheck == "pdf"){
        pdf(file,width = 20, height = 16)
      }
      if(!is.null(input$maf) & !is.null(input$sampleInfo)){
        return(Path()[[2]][[which(names(Path()[[2]]) == input$pl)]])
      }
      else{
        return(Path()[[2]][[which(names(Path()[[2]]) == input$pl)]])
      }
      dev.off()
    },
    contentType = paste('image/',input$DownloadPathPlotCheck,sep="")
  )
  output$DownloadSignaturePlot <- downloadHandler(
    filename = function() {
      paste("SignaturePlot",'.',input$DownloadSignaturePlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadSignaturePlotCheck == "png"){
        png(file,width = 1400, height = 800,res = 144)
      }
      else if (input$DownloadSignaturePlotCheckk == "pdf"){
        pdf(file,width = 1400, height = 800)
      }
      njtree <- inputNJtree()
      Meskit::treeMutationalSig(njtree,plot.Signatures = T)
      dev.off()
    },
    contentType = paste('image/',input$DownloadSignaturePlotCheck,sep="")
  )
  
  
})