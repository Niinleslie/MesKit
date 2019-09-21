options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))

#sider bar----

sidebar <- dashboardSidebar(
    width = 400,
  sidebarMenu(id="sidername",selected='home',
    menuItem(strong("Home"), tabName = "home", icon = icon("home")),
    menuItem(strong("Input Data"), tabName = "input", icon = icon("th",lib = "glyphicon")),
    menuItem(strong("ITH evaluation"), tabName = "ITH", icon = icon("location-arrow")),
    menuItem(strong("Clonal analysis"), tabName = "clone", icon = icon("th",lib = "glyphicon")),
    menuItem(strong("Functional analysis"), tabName = "function", icon = icon("bar-chart")),
    menuItem(strong("PhyloTree"), tabName = "Survival", icon = icon("line-chart"))
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

bodyIP <- tabItem("input",
                  h2('Input section'),
                  fluidRow(
                    column(
                      width = 3,
                      box(
                        title = div(shiny::icon("gear"), "Adjust image height & width", inline =TRUE),
                        width = NULL,
                        tabBox(
                          width = "100%",
                          height = "100%",
                          id = "InputSession",
                          title = "",
                          selected = "IS1",
                          tabPanel(
                            value = "IS1",
                            title = div(icon("book"), "Upload"),
                            radioGroupButtons(
                              "dataset",
                              label = "Data Set",
                              choices = c(default = "default", upload = "upload"),
                              selected = "default",
                              status = "primary"
                            ),
                            conditionalPanel(
                              condition = "input.dataset == 'upload'",
                              fileInput('maf','Upload maf file'),
                              fileInput('sampleInfo','Upload sampleInfo' ),
                              checkboxInput('useccf','choose if  use ccf',value = F),
                              conditionalPanel(
                                condition = "input.useccf == true",
                                fileInput('ccf.cluster','Upload ccf.cluster'),
                                fileInput('ccf.loci','Upload ccf.loci' )
                              )
                            ),
                            # fileInput('maf','Upload maf file'),
                            # fileInput('sampleInfo','Upload sampleInfo' ),
                            # checkboxInput('useccf','choose if  use ccf',value = F),
                            # conditionalPanel(
                            #   condition = "input.useccf == true",
                            #   fileInput('ccf.cluster','Upload ccf.cluster'),
                            #   fileInput('ccf.loci','Upload ccf.loci' )
                            # ),
                            actionBttn('submit1',div(
                              strong("Click ME to start analysing"),align = 'center',
                              icon("hand-right", lib = "glyphicon")))
                          ),
                          tabPanel(
                            value = "IS2",
                            title = div(icon("book"), "Setting"),
                            selectInput('ref','select reference genome(hg19/hg38)',
                                        choices = c(
                                          'hg19','hg38'
                                        )),
                            checkboxInput('use.indel','choose if  use indel',value = FALSE)
                          )
                        )
                      )
                    ),
                    column(
                      width = 9,
                      box(
                        title = "Maf Summary",
                        width = 1000,
                        height = 850,
                        # tabBox(
                        #   title = "",
                        #   selected = "MS2",
                        #   tabPanel(
                        #     id = "MS1",
                        #     title = div(icon("book"), "Readme")
                        #   ),
                        #   tabPanel(
                        #     id = "MS2",
                        #     title = div(icon("th", lib = "glyphicon"), "Maf Summary"),
                        #     withSpinner(plotOutput('mafSummary',height = "100%",width = "100%"))
                        #   )
                        # )
                        withSpinner(plotOutput('mafSummary'))
                      ),
                      box(
                        status = 'success',
                        width = NULL,
                        radioButtons('DownloadMafSummaryCheck','Choose file type to download:',
                                     c('png' ='png','pdf' = 'pdf'),inline = T),
                        downloadBttn('DownloadMafSummary', 'Download')
                      )
                      # box(
                      #   status = 'success',
                      #   width = NULL,
                      #   radioButtons('DownloadMafSummaryCheck','Choose file type to download:',
                      #                c('png' ='png','pdf' = 'pdf'),inline = T),
                      #   downloadBttn('DownloadMafSummary', 'Download')
                      # )
                    )
                  )
                  )


bodyITH <- tabItem("ITH",
                   h2('ITH evaluation'),
                   fluidRow(
                    column(
                     width = 3,
                     box(
                       width = NULL,
                       tabBox(
                         width = "100%",
                         height = "100%",
                         selected = "IE2",
                         tabPanel(
                           value = "IE1",
                           title = div(shiny::icon("gear"), "height&width", inline =TRUE),
                           sliderInput('width1','adjust image width',min = 700,max = 1100, value = 850)
                         ),
                         tabPanel(
                           value = "IE2",
                           title = div(shiny::icon("eye"), "setting&upload", inline =TRUE),
                           h4(strong("Setting for MutSharedPrivatedPlot")),
                           checkboxInput('show.num','show mutational number',value = FALSE),
                           h4(strong("Setting for StackPlot")),
                           fileInput('oncogeneListFile','Upload oncogeneListFile'),
                           fileInput('tsgListFile','Upload tsgListFile')
                         )
                       )
                       )
                       # title = div(shiny::icon("gear"), "Adjust image height & width", inline =TRUE),
                       # sliderInput('width1','adjust image width',min = 700,max = 1100, value = 850)                     )
                     ),
                     column(
                       width = 9,                     
                      box(
                       title = 'Visuliaztion',
                       width = NULL,
                       tabBox(                           
                         height = "100%", width = "100%",
                         selected = "caInput02",
                         side = "left",
                         tabPanel(
                           title = div(icon("chart-bar"), "mathScore"),
                           value = "caInput02",
                           withSpinner(dataTableOutput('mathScore')),
                           downloadBttn('DownloadMathScore', 'Download')
                         ),
                         tabPanel(
                           title = div(icon("image"), "vafPlot"),
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
                           title = div(icon("map"), "mutSharedPrivatePlot"),
                           value = "caInput04",
                           plotOutput("mut.share_private",height = "100%"),
                           box(
                             width = NULL,
                             status = 'success',
                             radioButtons('DownloadSharedPlotCheck','Choose file type to download:',
                                        c('png' ='png','pdf' = 'pdf'),inline = T),
                             downloadBttn('DownloadSharedPlot', 'Download')
                           )
                         ),
                         tabPanel(
                           title = div(icon("camara"), "stackPlot"),
                           value = "caInput05",
                           plotOutput("stackplot",height = "100%"),
                           box(
                             width = NULL,
                             status = 'success',
                             radioButtons('DownloadStackPlotCheck','Choose file type to download:',
                                          c('png' ='png','pdf' = 'pdf'),inline = T),
                             downloadBttn('DownloadStackPlot', 'Download')
                           )
                         ),
                         tabPanel(
                           title = div(icon("box"), "JaccardIndex"),
                           value = "caInput06",
                           plotOutput("JaccardIndex",height = "100%"),
                           box(
                             width = NULL,
                             status = 'success',
                             radioButtons('DownloadJaccardIndexCheck','Choose file type to download:',
                                          c('png' ='png','pdf' = 'pdf'),inline = T),
                             downloadBttn('DownloadJaccardIndex', 'Download')
                           )
                         )
                         )
                     )
                     )

                   )
)

bodyclone <- tabItem('clone',
                     h2('Clonal analysis'),
                    fluidRow(
                      column(
                        width = 3,
                      box(
                        status = 'primary',
                        width = NULL,
                        title = div(shiny::icon("gear"), "Adjust image height & width", inline =TRUE),
                        sliderInput('width2','adjust image width',min = 700,max = 1100, value = 850)
                      )
                      ),
                     column(
                       width = 9,
                      box(
                        title = 'Visuliaztion',
                        width = NULL,
                        tabBox(
                          selected = 'c01',
                          side = 'left',
                          height = "100%",
                          width = "100%",
                          tabPanel(
                            value = 'c01',
                            title = div(icon("newspaper"), "TumorClone Plot"),
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
                        h2('Functional analysis'),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                             width = NULL,
                             status = 'primary',
                             title = div(shiny::icon("gear"), "Adjust parameters", inline =TRUE),
                             sliderInput('pval','adjust pval',min = 0,max = 1, value = 0.05),
                             sliderInput('qval','adjust image qval',min = 0,max = 1, value = 0.2),
                             # checkboxInput('saveplot1',strong('click if save plot of result'),value = F),
                             # checkboxInput('writetable1',strong('click if write table of result'),value = F),
                             actionBttn('submit3',div(strong("Click me Start/Refresh"),align = 'center'))
                           ),
                          box(
                             width = NULL,
                             status = 'primary',
                             title = div(shiny::icon("gear"), "Data upload & Configuration", inline =TRUE),
                             sliderInput('width3','adjust image width',min = 700,max = 1100, value = 850),
                             sliderInput('height3','adjust image height',min = 600,max = 2000, value = 1300)
                          )
                          ),
                          column(
                           width = 9,
                           box(
                            title = 'Visuliaztion',
                            width = NULL,
                            height = 2100,
                            tabBox(
                              side = 'left',
                              selected = 'F02',
                              width = "100%",
                              height = "100%",
                              tabPanel(
                                title = div(icon("lightbulb"), "GO Analyse"),
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
                                title = div(icon("microsoft"), "Pathway Analyse"),
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
                        h2('Phylotree Visualiaztion'),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                            width = NULL,
                            status = 'primary',
                            title = div(shiny::icon("gear"), "Adjust parameters", inline =TRUE),
                            radioGroupButtons('phyloTreeType','select tye of your phylotree',
                                         c( 'NJtree' = 'njtree',
                                            'Newick' = 'newick',
                                            'beast' = 'beast',
                                            'PAML' = 'PAML'),
                                         status = "primary",
                                         selected = 'njtree'),
                            conditionalPanel(
                              condition = "input.phyloTreeType != 'njtree'",
                              fileInput('phylotree.dir','upload your phylotree file')
                            ),
                            checkboxInput('show.mutSig','choose if show mutationignature',value = TRUE),
                            checkboxInput('show.heatmap','choose if show heatmap',value = TRUE),
                            radioGroupButtons(
                              inputId = "sig.name",
                              label = "signature name",
                              choices = c(default = "Default", alias = "Alias"),
                              selected = "Default",
                              status = "primary"
                            ),
                            radioGroupButtons(
                              inputId = "heatmap.type",
                              label = "heatmap type",
                              choices = c(binary = "binary", CCF = "CCF"),
                              selected = "binary",
                              status = "primary"
                            )
                              # actionBttn('submit.tree',div(
                              # strong("Click ME to visualize result"),align = 'center',
                              # icon("hand-right", lib = "glyphicon")))
                          )
                          ),
                          column(
                            width = 9,
                            box(
                              title = "visulization",
                              width = NULL,
                              height = 1000,
                              tabBox(
                                side = 'left',
                                selected = 'S01',
                                width = "100%",
                                height = "100%",
                                tabPanel(
                                  title = div(icon("lightbulb"), "Phylotree"),
                                  value = "S01",
                                  withSpinner(plotOutput("phylotree",height = 700)),
                                  box(
                                    width = NULL,
                                    status = 'success',
                                    radioButtons('DownloadPhyloTreeCheck','Choose file type to download:',
                                                 c('png' ='png','pdf' = 'pdf'),inline = T),
                                    downloadBttn('DownloadPhyloTree', 'Download')
                                  )
                                ),
                                tabPanel(
                                  title = div(icon("align-left"), "signature"),
                                  value = "S02",
                                  withSpinner(plotOutput("signature",height = 700)),
                                  box(
                                    width = NULL,
                                    status = 'success',
                                    radioButtons('DownloadSignaturePlotCheck','Choose file type to download:',
                                                 c('png' ='png','pdf' = 'pdf'),inline = T),
                                    downloadBttn('DownloadSignaturePlot', 'Download')
                                  )
                                )
                            )
                            )

                            # box(
                            # title = "Visulization",
                            # width = NULL,
                            # height = 800,
                            # withSpinner(plotOutput("phylotree"))), 
                            # box(
                            #   width = NULL,
                            #   status = 'success',
                            #   radioButtons('DownloadPhyloTreeCheck','Choose file type to download:',
                            #              c('png' ='png','pdf' = 'pdf'),inline = T),
                            #   downloadBttn('DownloadPhyloTree', 'Download')
                            # )
                          )
                        )
                        )
                    

#Main function----
shinyUI(
  dashboardPage(skin = "blue",
    dashboardHeader(title = "Meskit: Analysis and visualize multi-sample whole-exome sequencing data",
                  titleWidth = 650),
    sidebar,
    dashboardBody(
      theme = 
      tags$head(
        tags$style(HTML(".shiny-output-error-validation {color: brown;}")),
        tags$link(rel = "stylesheet", type = "text/css", href = "css/main.css")
      ),
      tabItems(
        bodyHome,
        bodyIP,
        bodyITH,
        bodyclone,
        bodyfunction,
        bodySurvival
      )
    )
  )
)
