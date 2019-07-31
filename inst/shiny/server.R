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

  upload_maf <- eventReactive(input$submit1,{
    if(input$useccf){
      read.Maf(patientID = input$patientid,maf.dir =input$maf$datapath, sample_info.dir = input$sampleInfo$datapath,
               ccf.cluster.dir = input$ccf.cluster$datapath, ccf.loci.dir = input$ccf.loci$datapath)
    }
    else{
      read.Maf(patientID = input$patientid,maf.dir = input$maf$datapath, sample_info.dir = input$sampleInfo$datapath)
    }
    })
  
  upload_njtree <- eventReactive(input$submit3,{
    maf <- upload_maf()
    njtree <- NJtree(maf,use.indel = input$use.indel,use.ccf = input$useccf )
  })
  
  output$mafSummary <- renderPlot({
    maf <- upload_maf()
    plotmafSummary(maf = maf)
  },
  width = width1,
  height = height1
  )
  output$mathScore <- renderDataTable({
    maf <- upload_maf()
    MATH_score(maf)
  })
  output$vaf.cluster <- renderPlot({
    maf <- upload_maf()
    VAF_plot(maf)
    }, 
    width = width1,
    height = height1
  )
  output$mut.share_private <- renderPlot({
    maf <- upload_maf()
    mutSharedPrivate(maf)
  },
  width = width1,
  height = height1
  )
  output$cloneplot <- renderPlot({
    validate(
      need(!(is.null(input$ccf.cluster1$datapath)), "Upload ccf.cluster")
    )
    validate(
      need(!(is.null(input$ccf.loci1$datapath)), "Upload ccf.loci ")
    )
    validate(
      need(input$okk != 0, "Press"))
    TumorClones_plot.shiny(patientID = input$patientid,
                           ccf.cluster.dir = input$ccf.cluster1$datapath,
                           ccf.loci.dir = input$ccf.loci1$datapath,
                           out.dir = '.')
  },
  width = width2,
  height = height2
  )
  
  output$GOplot <- renderPlot({
    njtree <- upload_njtree()
    GO.njtree.shiny(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
  },
  width = width3,
  height = height3
  )
  output$Pathwayplot <- renderPlot({
    njtree <- upload_njtree()
    Pathway.njtree.shiny(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
  },
  width = width3,
  height = height3
  )
  
  output$phylotree <- renderPlot({
    if(input$phyloTreeType == 'njtree'){
      maf <- upload_maf()
      validate(
        need(!(input$submit.tree),"Press tht button if you choose parameter")
      )
      plot.PhyloTree(maf,use.indel = input$use.indel,heatmap.type = input$heatmap.type ,show.mutSig = input$show.mutSig,
                show.heatmap = input$show.heatmap, output.dir = '',phylotree.type = input$phyloTreeType)
    }
    else{
      validate(
        need(!is.null(input$phylotree.dir),"Upload your phylotree file")
      )
      validate(
        need((input$submit.tree),"press button") 
      )
      plot.PhyloTree(phylotree.dir = input$phylotree.dir$datapath,phylotree.type = input$phyloTreeType,output.dir = '')
    }
    
  })
  
  output$DownloadMafSummary <- downloadHandler(
    
    filename = function() {
      paste("MafSummary",'.',input$DownloadMafSummaryCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadMafSummaryCheck == "png"){
        png(file,width = width1() , height = height1())
      }
      else if (input$DownloadMafSummaryCheck == "pdf"){
        pdf(file,width = 14 , height = 14)
      }
      maf <- upload_maf()
      plotmafSummary(maf = maf)
      dev.off()
    },
    contentType = paste('image/',input$DownloadMafSummaryCheck,sep="")
  )
  
  output$DownloadVafPlot <- downloadHandler(
    filename = function() {
      paste("VafPlot", Sys.time(), '.',input$DownloadVafPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadVafPlotCheck == "png"){
        png(file,width = width1() , height = width1())
      }
      else if (input$DownloadVafPlotCheck == "pdf"){
        pdf(file,width = 12 , height = 9)
      }
      maf <- upload_maf()
      VAF_plot(maf = maf)
      dev.off()
    },
    contentType = paste('image/',input$DownloadVafPlotCheck,sep="")
  )
  
  output$DownloadSharedPlot <- downloadHandler(
    filename = function() {
      paste("SharedPlot", Sys.time(), '.',input$DownloadSharedPlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadSharedPlotCheck == "png"){
        png(file,width = width1() , height = height1())
        
      }
      else if (input$DownloadSharedPlotCheck == "pdf"){
        pdf(file,width = 13 , height = 10)
      }
      dev.off()
    }
  )
  output$DownloadClonePlot <- downloadHandler(
    filename = function() {
      paste("ClonePlot", Sys.time(), '.',input$DownloadClonePlotCheck, sep='')
    },
    content = function(file) {
      if (input$DownloadClonePlotCheck == "png"){
        png(file,width = width2() , height = height2())
      }
      else if (input$DownloadClonePlotCheck == "pdf"){
        pdf(file,width = 9 , height = 6.5)
      }
      maf <- upload_maf()
      TumorClones_plot(maf = maf)
      dev.off()
    },
    contentType = paste('image/',input$DownloadClonePlotCheck,sep="")
  )
output$DownloadPhyloTree <- downloadHandler(
  filename = function() {
    paste("PhyloTree", Sys.time(), '.',input$DownloadPhyloTreeCheck, sep='')
  },
  content = function(file) {
    if (input$DownloadPhyloTreeCheck == "png"){
      png(file,width = "1100px", height = "850px")
    }
    else if (input$DownloadPhyloTreeCheck == "pdf"){
      pdf(file,width = 14, height = 8)
    }
    maf <- upload_maf()
    plot.PhyloTree(maf = maf)
    dev.off()
  },
  contentType = paste('image/',input$DownloadPhyloTreeCheck,sep="")
)

output$DownloadGOPlot <- downloadHandler(
  filename = function() {
    paste("GOPlot", Sys.time(), '.',input$DownloadGOPlotCheck, sep='')
  },
  content = function(file) {
    if (input$DownloadGOPlotCheck == "png"){
      png(file,width = width3(), height = width3())
    }
    else if (input$DownloadGOPlotCheck == "pdf"){
      pdf(file,width = width3()/100, height = height3()/100)
    }
    njtree <- upload_njtree()
    GO.njtree.shiny(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
    dev.off()
  },
  contentType = paste('image/',input$DownloadGOPlotCheck,sep="")
)
output$DownloadPathPlot <- downloadHandler(
  filename = function() {
    paste("EnrichPlot",'.',input$DownloadGOPlotCheck, sep='')
  },
  content = function(file) {
    if (input$DownloadPathPlotCheck == "png"){
      png(file,width = width3(), height = width3())
    }
    else if (input$DownloadPathPlotCheck == "pdf"){
      pdf(file,width = width3()/100, height = height3()/100)
    }
    njtree <- upload_njtree()
    Pathway.njtree.shiny(njtree, savePlot = F, qval = input$qval ,pval = input$pval)
    dev.off()
  },
  contentType = paste('image/',input$DownloadPathPlotCheck,sep="")
)


})