options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))

#sider bar----

sidebar <- dashboardSidebar(
    width = 300,
  sidebarMenu(id="sidername",selected='home',
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("ITH evaluation", tabName = "ITH", icon = icon("location-arrow")),
    menuItem("Clonal analysis", tabName = "clone", icon = icon("th",lib = "glyphicon")),
    menuItem("Functional analysis", tabName = "function", icon = icon("bar-chart")),
    menuItem("PhyloTree", tabName = "Survival", icon = icon("line-chart"))
  )
)


#bodyHome ----
#tabItem:创建标签页；子页面；与sidebar里的tabName相对应
#fluidRow代表一行，
#box是基本容器，可以存放图形或其他输出内容(非必须；但无法设定宽度，maxWidth = 12)
#shinydashboar
bodyHome <- tabItem("home",
              fluidRow(
                box(
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  title = strong("Wellcome to the MesKit reporter"),
                  includeMarkdown("dom/Introduction.md")
                )
              ),
                    
              fluidRow(
                box(
                  width = 6,
                  status = "info",
                  solidHeader = TRUE,
                  title =  strong("Overview of MesKit package"),
                  p("The typical workflow begins with MAF object creation by reading an MAF file combind with sample information. Based on Maf object, both ITH assessment and clonal analysis can be conducted. Furthermore, MesKit can perform function analysis and mutation analsysi on njtree object, which is converted from Maf object.",
                    style = "font-si16pt"),
                  br(),
                  div(img(src = "images/pipeline.png", width=600), style="text-align: center;")
                  ),
                box(
                  width = 6,
                  status = "info",
                  solidHeader = TRUE,
                  title =  strong("MesKit result viewer"),
                  includeMarkdown("dom/Results_viewer.md"),
                  br(),
                  div(img(src = "images/result_view.png", width=730, height = 422), style="text-align: center;")
                  #img(src = "images/result_view.png", align = "center", width="100%")
                  )

                )
)


bodyITH <- tabItem("ITH",
                   h2(strong('ITH evaluation')),
                   fluidRow(
                     column(
                       width = 3,
                      box(
                       status = 'primary',
                       width = NULL,
                       id = 'upload',
                       h3(strong('Upload file /Setting parameters')),
                       fileInput('maf','Upload your maf file'),
                       fileInput('sampleInfo','Upload sampleInfo' ),
                       textInput('patientid','Input patient ID'),
                       selectInput('ref','select reference genome',
                                   choices = c(
                                     'hg19','hg38'
                                   )),
                       checkboxInput('useccf','choose if  use ccf',value = F),
                       conditionalPanel(
                         condition = "input.useccf == true",
                         fileInput('ccf.cluster','Upload ccf.cluster'),
                         fileInput('ccf.loci','Upload ccf.loci' )
                       ),
                       actionBttn('submit1',div(
                         strong("Click ME to start analysing"),align = 'center',
                         icon("hand-right", lib = "glyphicon")))
                     ),
                     box(
                      width = NULL,
                      status = 'primary',
                      h3(strong('Adjust image height/width')),
                      sliderInput('width1','adjust image width',min = 700,max = 1100, value = 850),
                       sliderInput('height1','adjust image height',min = 400,max = 570, value = 560)
                     )
                     ),
                     column(
                       width = 9,                     
                      box(
                       h3(strong('Visuliaztion')),
                       width = NULL,
                       height = 900,
                       tabBox(                           
                         height = "100%", width = "100%",
                         selected = "caInput01",
                         side = "left",
                         tabPanel(
                           title = "MafSummary",
                           value = "caInput01",
                           withSpinner(plotOutput('mafSummary',height = '100%')),
                           box(
                             status = 'success',
                             width = NULL,
                             radioButtons('DownloadMafSummaryCheck','Choose file type to download:',
                                        c('png' ='png','pdf' = 'pdf'),inline = T),
                             downloadBttn('DownloadMafSummary', 'Download')
                           )
                         ),
                         tabPanel(
                           title = "Math.Score",
                           value = "caInput02",
                           withSpinner(dataTableOutput('mathScore'))                 
                         ),
                         tabPanel(
                           title = "VAF plot",
                           value = "caInput03",
                           plotOutput("vaf.cluster",height = "100%",width = "100%"),
                           box(
                             width = NULL,
                             status = 'success',
                             radioButtons('DownloadVafPlotCheck','Choose file type to download:',
                                        c('png' ='png','pdf' = 'pdf'),inline = T),
                             downloadBttn('DownloadVafPlot', 'Download')
                           )
                         ),
                         tabPanel(
                           title = "SharePrivate",
                           value = "caInput04",
                           plotOutput("mut.share_private",height = "100%"),
                           box(
                             width = NULL,
                             status = 'success',
                             radioButtons('DownloadSharePlotCheck','Choose file type to download:',
                                        c('png' ='png','pdf' = 'pdf'),inline = T),
                             downloadBttn('DownloadSharePlot', 'Download')
                           )

                         )
                         )
                     )
                     )

                   )
)

bodyclone <- tabItem('clone',
                     h2(strong('Clonal analysis')),
                    fluidRow(
                      column(
                        width = 3,
                      box(
                        status = 'primary',
                        width = NULL,
                        h3(strong("Upload ccf")),
                        fileInput('ccf.cluster1','Upload ccf.cluster'),
                        fileInput('ccf.loci1','Upload ccf.loci' ),
                        actionBttn('okk',div(
                          strong("Click ME to visualize result"),align = 'center',
                          icon("hand-right", lib = "glyphicon")))
                      ),
                      box(
                        status = 'primary',
                        width = NULL,
                        h3(strong("Adjust image")),
                        sliderInput('width2','adjust image width',min = 700,max = 1100, value = 850),
                        sliderInput('height2','adjust image height',min = 400,max = 570, value = 560)
                      )
                      ),
                     column(
                       width = 9,
                      box(
                        h3(strong('Visuliaztion')),
                        width = NULL,
                        tabBox(
                          selected = 'c01',
                          side = 'left',
                          height = "100%",
                          width = "100%",
                          tabPanel(
                            value = 'c01',
                            title = "TumorClonePlot",
                            withSpinner(plotOutput('cloneplot',height = "100%")),
                            box(
                              width = NULL,
                              status = 'success',
                              radioButtons('DownloadClonePlotCheck','Choose file type to download:',
                                     c('png' ='png','pdf' = 'pdf'),inline = T),
                              downloadBttn('DownloadClonePlot', 'Download')
                            )

                          )
                        )
                      )
                     )    

                    )
)

bodyfunction <- tabItem('function',
                        h2(strong('Functional analysis')),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                             width = NULL,
                             status = 'primary',
                             h3(strong('Adjust pval/qval')),
                             sliderInput('pval','adjust pval',min = 0,max = 1, value = 0.05),
                             sliderInput('qval','adjust image qval',min = 0,max = 1, value = 0.2),
                             actionBttn('submit3',div(strong("Click me Start/Refresh"),align = 'center'))
                           ),
                          box(
                             width = NULL,
                             status = 'primary',
                             h3(strong('Adjust image height/width')),
                             sliderInput('width3','adjust image width',min = 700,max = 1100, value = 850),
                             sliderInput('height3','adjust image height',min = 800,max = 2000, value = 1300)
                          )
                          ),
                          column(
                           width = 9,
                           box(
                             h3(strong('Visuliaztion')),
                            width = NULL,
                            height = 2100,
                            tabBox(
                              side = 'left',
                              selected = 'F02',
                              width = "100%",
                              height = "100%",
                              tabPanel(
                                title = 'GO analysis',
                                value = 'F01',
                                withSpinner(plotOutput('GOplot',height = "100%",width = "100%")),
                                box(
                                  width = NULL,
                                  status = 'success',
                                  radioButtons('DownloadGOPlotCheck','Choose file type to download:',
                                             c('png' ='png','pdf' = 'pdf'),inline = T),
                                  downloadBttn('DownloadGOPlot', 'Download')
                                )

                              ),
                              tabPanel(
                                title = 'pathway analysis',
                                value = 'F02',
                                withSpinner(plotOutput('Pathwayplot',height = "100%",width = "100%")),
                                box(
                                  width = NULL,
                                  status = 'success',
                                  radioButtons('DownloadPathPlotCheck','Choose file type to download:',
                                             c('png' ='png','pdf' = 'pdf'),inline = T),
                                  downloadBttn('DownloadPathPlot', 'Download')
                                )

                              )
                            )
                          )
                          )
                         

                        ))

bodySurvival <- tabItem('Survival',
                        h2(strong('Phylotree Visualiaztion')),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                            width = NULL,
                            status = 'primary',
                            radioButtons('phyloTreeType',h3(strong('select tye of your phylotree ')),
                                         c( 'NJtree(default)' = 'njtree',
                                            'Newick' = 'newick',
                                            'beast' = 'beast',
                                            'PAML' = 'PAML'),selected = 'njtree',inline = T),
                            h3(strong("parameters")),
                            checkboxInput('use.indel','choose if  use indel',value = FALSE),
                            checkboxInput('show.mutSig','choose if show mutationignature',value = TRUE),
                            checkboxInput('show.heatmap','choose if show heatmap',value = TRUE),
                            textInput('heatmap.type','heatmap type','binary'),
                            conditionalPanel(
                              condition = "input.phyloTreeType != 'njtree'",
                              fileInput('phylotree.dir','upload your phylotree file')
                            ),
                              actionBttn('submit.tree',div(
                              strong("Click ME to visualize result"),align = 'center',
                              icon("hand-right", lib = "glyphicon")))
                          )
                          ),
                          column(
                            width = 9,
                            box(
                            width = NULL,
                            height = 1200 ,
                            withSpinner(plotOutput("phylotree",height = 1000)), 
                            box(
                              width = NULL,
                              status = 'success',
                              radioButtons('DownloadPhyloTreeCheck','Choose file type to download:',
                                         c('png' ='png','pdf' = 'pdf'),inline = T),
                              downloadBttn('DownloadPhyloTree', 'Download')
                            )

                          )
                          )

                        ))
                    

#Main function----
shinyUI(
  dashboardPage(skin = "blue",
    dashboardHeader(title = "Meskit: Analysis and visualize multi-sample whole-exome sequencing data",
                  titleWidth = 650),
    sidebar,
    dashboardBody(
      tags$head(
        tags$style(HTML(".shiny-output-error-validation {color: brown;}")),
        tags$link(rel = "stylesheet", type = "text/css", href = "css/main.css")
      ),
      tabItems(
        bodyHome,
        bodyITH,
        bodyclone,
        bodyfunction,
        bodySurvival
      )
    )
  )
)
