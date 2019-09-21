suppressMessages(library(shiny))



# Define server logic required to plot various variables against mpg
shinyServer(function(input, output){
  
  phylotree.type <- reactive({
    return(input$phylotTreeType)
  })
  width1 <- reactive({
    return(input$width1)
  })
  height1 <- reactive({
    return(input$height1)
  })
  width2 <- reactive({
    return(input$width2)
  })
  height2 <- reactive({
    return(input$height2)
  })
  width3 <- reactive({
    return(input$width3)
  })
  height3 <- reactive({
    return(input$height3)
  })

  inputData <- eventReactive(input$submit1,{
    if(input$dataset == "default"){
      example.object <- readRDS("./example/MesKit_Example.rds")
    }
    else{
      
      if(!is.null(input$ccf.cluster)&!is.null(input$ccf.loci)){
        maf <- readMaf(mafFile = input$maf$datapath,
                       sampleInfoFile = input$sampleInfo$datapath, 
                       ccfClusterTsvFile =  input$ccf.cluster$datapath, 
                       ccfLociTsvFile = input$ccf.loci$datapath)
      }
      else{
        maf <- readMaf(mafFile = input$maf$datapath, 
                       sampleInfoFile = input$sampleInfo$datapath)
      }
    }
    switch(
      input$dataset,
      "default" = example.object,
      "upload"  = maf,
    )
    })
  
  inputNJtree <- reactive({
    if(input$dataset == "upload"){
      maf <- inputData()
      njtree <- getNJtree(maf, use.indel = input$use.indel)
    }
  })
  
  output$mafSummary <- renderPlot({
    if(input$dataset == "upload"){
      maf <- inputData()
      plotmafSummary(maf = maf)
    }
    else{
      maf <- inputData()$maf
      plotmafSummary(maf = maf)
    }
  },
  width = 1000,
  height = 800,
  res = 120
  )
  output$mathScore <- renderDataTable({
    if(input$dataset == "upload"){
      maf <- inputData()
      mathScore(maf)$sampleLevel
    }
    else{
      inputData()$mathscore$sampleLevel
    }
  })
  output$vaf.cluster <- renderPlot({
    if(input$dataset == "upload"){
      maf <- inputData()
      vafCluster(maf)
    }
    else{
      inputData()$vafplot
    }
  }, 
    width = width1,
    height = 560,
    res = 100
  )
  output$mut.share_private <- renderPlot({
    if(input$dataset == "upload"){
      maf <- inputData()
      mutSharedPrivate(maf,show.num = input$show.num)
    }
    else{
      inputData()$privateplot
    }
  },
  width = width1,
  height = 560,
  res = 100
  )
  output$stackplot <- renderPlot({
    if(input$dataset == "upload"){
      maf <- inputData()
      validate(
        need(!((is.null(input$oncogeneListFile$datapath) & is.null(input$tsgListFile$datapath))), 
             "Upload oncogeneListFile and tsgListFile in 'Setting&Upload' ")
      )
      mutStackPlot(maf, oncogeneListFile = input$oncogeneListFile$datapath,
                   tsgListFile = input$tsgListFile$datapath, themeOption="npg", show.percentage = TRUE)
    }
    else{
      inputData()$stackplot
    }
  },
  width = width1,
  height = 560,
  res = 100
  )
  output$JaccardIndex <- renderPlot({
    if(input$dataset == "upload"){
      maf <- inputData()
      JaccardIndex(maf)
    }
    else{
      maf <- inputData()$maf
      JaccardIndex(maf)
    }
  },
  width = width1,
  height = 560,
  res = 100
  )
  output$cloneplot <- renderPlot({
    if(input$dataset == "upload"){
      validate(
        need(!(is.null(input$ccf.cluster$datapath)), "click the button 'use ccf',Upload ccf.cluster in Session 'Input Data' ")
      )
      validate(
        need(!(is.null(input$ccf.loci$datapath)), "Upload ccf.loci Session 'Input Data'")
      )
      maf <- inputData()
      tumorClonesPlot(maf)
    }
    else{
      maf <- inputData()$maf
      tumorClonesPlot(maf)
    }
  },
  width = width2,
  height = 560,
  res = 100
  )
  
  output$GOplot <- renderPlot({
    if(input$dataset == "upload"){
      njtree <- inputNJtree()
      GO.njtree(njtree, savePlot = input$saveplot1, writeTable = input$writetable1, qval = input$qval ,pval = input$pval)
    }
    else{
      plot(inputData()$GOplot)
    }
  },
  width = width3,
  height = height3
  )
  output$Pathwayplot <- renderPlot({
    if(input$dataset == "upload"){
      njtree <- inputNJtree()
      Pathway.njtree(njtree, savePlot = input$saveplot1, writeTable = input$writetable1, qval = input$qval ,pval = input$pval)
    }
    else{
      plot(inputData()$Pathwayplot)
    }
  },
  width = width3,
  height = height3,
  res = 100
  )
  output$phylotree <- renderPlot({
    if(input$dataset == "upload"){
      if(input$phyloTreeType == 'njtree'){
        njtree <- inputNJtree()
        if(!is.null(input$useccf)){
          validate(
            need(input$heatmap.type == "CCF","switch heatmap type to CCF")
          )
        }
        p <- plotPhyloTree(njtree, phylotree.type = input$phyloTreeType, 
                           heatmap.type = input$heatmap.type, sig.name = input$sig.name,
                           show.mutSig = input$show.mutSig, show.heatmap = input$show.heatmap)
        return(p)
      }
      else{
        validate(
          need(!is.null(input$phylotree.dir),"Upload your phylotree file")
        )
        p <- plotPhyloTree(phylotree.dat = input$phylotree.dir$datapath, 
                           phylotree.type = input$phyloTreeType)
        return(p)
      }
    }
    else{
      inputData()$phylotreeplot
    }
  },
  res = 100
)
  output$signature <- renderPlot({
    if(input$dataset == "upload"){
      if(input$phyloTreeType == 'njtree'){
        njtree <- inputNJtree()
        treeMutationalSig(njtree,plot.Signatures = T)
      }
      else{
        message("Drawing signature plot when using njtree")
      }
    }
    else{
      njtree <- inputData()$njtree
      treeMutationalSig(njtree,plot.Signatures = T)
    }
  },
  res = 100
  )
  output$DownloadMafSummary <- downloadHandler(
    
    filename = function() {
      paste("MafSummary",".",input$DownloadMafSummaryCheck,sep='')
    },
    content = function(file) {
      if (input$DownloadMafSummaryCheck == "png"){
        png(file,width = 1400 , height =1400,res = 144)
      }
      else if (input$DownloadMafSummaryCheck == "pdf"){
        pdf(file,width = 14 , height = 14)
      }
      if(input$dataset == "upload"){
        maf <- inputData()
        plotmafSummary(maf = maf)
      }
      else{
        maf <- inputData()$maf
        plotmafSummary(maf = maf)
      }
      dev.off()
    }
  )
  output$DownloadMathScore <- downloadHandler(
    filename = function() {
      paste("MathScore",'.',"csv", sep='')
    },
    content = function(file) {
      maf <- inputData()
      data <- mathScore(maf = maf)
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
      if(input$dataset == "upload"){
        maf <- inputData()
        vafCluster(maf)
      }
      else{
        print(inputData()$vafplot)
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
      if(input$dataset == "upload"){
        mutStackPlot(maf, oncogeneListFile = input$oncogeneListFile$datapath,
                     tsgListFile = input$tsgListFile$datapath, themeOption="npg", show.percentage = TRUE)
      }
      else{
        print(inputData()$stackplot)
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
      if(input$dataset == "upload"){
        maf <- inputData()
        JaccardIndex(maf)
      }
      else{
        maf <- inputData()$maf
        JaccardIndex(maf)
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
      if(input$dataset == "upload"){
        maf <- inputData()
        mutSharedPrivate(maf)
      }
      else{
        print(inputData()$privateplot) 
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
      if(input$dataset == "upload"){
        maf <- inputData()
        tumorClonesPlot(maf)
      }
      else{
        maf <- inputData()$maf
        tumorClonesPlot(maf)
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
    if(input$dataset == "upload"){
      maf <- inputData()
      plotPhyloTree(maf)
    }
    else{
      print(inputData()$phylotreeplot)
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
    if(input$dataset == "upload"){
      njtree <- inputNJtree()
      GO.njtree(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
    }
    else{
      plot(inputData()$GOplot)
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
    if(input$dataset == "upload"){
      njtree <- upload_njtree()
      Pathway.njtree.shiny(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
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
    njtree <- upload_njtree()
    treeMutationalSig(njtree,plot.Signatures = T)
    dev.off()
  },
  contentType = paste('image/',input$DownloadSignaturePlotCheck,sep="")
)


})