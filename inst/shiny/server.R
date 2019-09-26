suppressMessages(library(shiny))



# Define server logic required to plot various variables against mpg
shinyServer(function(input, output){
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
    if(input$dataset == "default"){
      maf_example <- readRDS("./example/MesKit_Example.rds")
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
      "default" = maf_example,
      "upload"  = maf,
    )
    })
  
  inputNJtree <- reactive({
    if(input$dataset == "upload"){
      maf <- inputData()
      njtree <- getNJtree(maf, use.indel = input$use.indel,
                          ccf.mutation.id = input$ccf.mutation.id,ccf.mutation.sep = input$ccf.mutation.sep)
    }
  })

  ms <- eventReactive(input$submit2,{
    if(input$dataset == "upload"){
      maf <- inputData()
      mathScore(maf,tsb = input$tsb,
                minvaf = input$vafrange[1],maxvaf = input$vafrange[2])$sampleLevel
    }
    else{
      maf <- inputData()
      mathScore(maf,tsb = input$tsb,
                minvaf = input$vafrange[1],maxvaf = input$vafrange[2])$sampleLevel
    }
  })
  output$mathScore <- renderDataTable({
     ms()
  })
  vc <- eventReactive(input$submit3,{
    if(input$dataset == "upload"){
      maf <- inputData()
      vafCluster(maf,plotOption = input$plotOption)
    }
    else{
      maf <- inputData()
      vafCluster(maf,plotOption = input$plotOption)
    }
  })
  output$vaf.cluster <- renderPlot({
    vc()
  }, 
    width = width1,
    height = 560,
    res = 100
  )
  msp <- eventReactive(input$submit4,{
    if(input$dataset == "upload"){
      maf <- inputData()
      mutSharedPrivate(maf,show.num = input$show.num)
    }
    else{
      maf <- inputData()
      mutSharedPrivate(maf,show.num = input$show.num)
    }
  })
  output$mut.share_private <- renderPlot({
    msp()
  },
  width = width2,
  height = 560,
  res = 100
  )
  stk <- eventReactive(input$submit5,{
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
      maf <- inputData()
      validate(
        need(!((is.null(input$oncogeneListFile$datapath) & is.null(input$tsgListFile$datapath))), 
             "Upload oncogeneListFile and tsgListFile in 'Setting&Upload' ")
      )
      mutStackPlot(maf, oncogeneListFile = input$oncogeneListFile$datapath,
                   tsgListFile = input$tsgListFile$datapath, themeOption="npg", show.percentage = TRUE)
    }
  })
  output$stackplot <- renderPlot({
     stk()
  },
  width = width3,
  height = 560,
  res = 100
  )
  ji <- eventReactive(input$submit6,{
    if(input$dataset == "upload"){
      maf <- inputData()
      JaccardIndex(maf)
    }
    else{
      maf <- inputData()
      JaccardIndex(maf)
    }
  })
  output$JaccardIndex <- renderPlot({
     ji()
  },
  width = width4,
  height = 560,
  res = 100
  )
  clp <- eventReactive(input$submit7,{
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
      maf <- inputData()
      tumorClonesPlot(maf)
    }
  })
  output$cloneplot <- renderPlot({
      clp()
  },
  width = width5,
  height = 560,
  res = 100
  )
  GO <- eventReactive(input$submit8,{
    if(input$dataset == "upload"){
      njtree <- inputNJtree()
      GO.njtree(njtree, qval = input$qval1 ,pval = input$pval1)
    }
    else{
      njtree <- inputNJtree()
      GO.njtree(njtree, qval = input$qval1 ,pval = input$pval1)
    }
  })
  output$chooselist1 <- renderUI({
    selectInput("gl","chooose branch to show",
                choices = names(GO()[[2]]) ,selected = names(GO()[[2]])[1],width = 600)
  })
  output$GOplot <- renderPlot({
    if(input$dataset == "upload"){
      return(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
    }
    else{
      return(GO()[[2]][[which(names(GO()[[2]]) == input$gl)]])
    }
  },
  width = width6,
  height = height6
  )
  Path <- eventReactive(input$submit9,{
    if(input$dataset == "upload"){
      njtree <- inputNJtree()
      list <- Pathway.njtree(njtree, qval = input$qval2 ,pval = input$pval2)
      return(list)
    }
    else{
      njtree <- inputNJtree()
      list <- Pathway.njtree(njtree, qval = input$qval2 ,pval = input$pval2)
      return(list)
    }
  })
  output$chooselist2 <- renderUI({
    selectInput("pl","chooose branch to show",
                choices = names(Path()[[2]]) ,selected = names(Path()[[2]])[1],width = 600)
  })
  output$Pathwayplot <- renderPlot({
    if(input$dataset == "upload"){
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
  pht <- eventReactive(input$submit10,{
    if(input$dataset == "upload"){
      if(input$phyloTreeType == 'njtree'){
        njtree <- inputNJtree()
        if(input$useccf == T){
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
      njtree <- inputData()$njtree
      p <- plotPhyloTree(njtree, phylotree.type = input$phyloTreeType, 
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
  sigp <- eventReactive(input$submit11,{
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
  })
  output$signature <- renderPlot({
     sigp()
  },
  res = 100
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
        maf <- inputData()
        vafCluster(maf)
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
        mutStackPlot(maf, oncogeneListFile = input$oncogeneListFile$datapath,
                     tsgListFile = input$tsgListFile$datapath, themeOption="npg", show.percentage = TRUE)
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
        maf <- inputData()
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
        maf <- inputData()
        mutSharedPrivate(maf) 
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
        maf <- inputData()
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
      if(input$phyloTreeType == 'njtree'){
        njtree <- inputNJtree()
        if(input$useccf == T){
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
      njtree <- inputData()$njtree
      p <- plotPhyloTree(njtree, phylotree.type = input$phyloTreeType, 
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
    if(input$dataset == "upload"){
      njtree <- inputNJtree()
      GO.njtree(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
    }
    else{
      njtree <- inputNJtree()
      GO.njtree(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
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

    njtree <- upload_njtree()
    Pathway.njtree(njtree, savePlot = F, qval = input$qval ,pval = input$pval)


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