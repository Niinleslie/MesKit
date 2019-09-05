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
                            textInput('patientid','Input patient ID'),
                            selectInput('ref','select reference genome(hg19/hg38)',
                                        choices = c(
                                          'hg19','hg38'
                                        ))
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
                       status = 'primary',
                       title = div(shiny::icon("gear"), "Adjust image height & width", inline =TRUE),
                       sliderInput('width1','adjust image width',min = 700,max = 1100, value = 850)                     )
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
                           title = div(icon("chart-bar"), "Math Score"),
                           value = "caInput02",
                           withSpinner(dataTableOutput('mathScore')),
                           downloadBttn('DownloadMathScore', 'Download')
                         ),
                         tabPanel(
                           title = div(icon("image"), "VAF Plot"),
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
                           title = div(icon("map"), "mut Shared private plot"),
                           value = "caInput04",
                           plotOutput("mut.share_private",height = "100%"),
                           box(
                             width = NULL,
                             status = 'success',
                             radioButtons('DownloadSharedPlotCheck','Choose file type to download:',
                                        c('png' ='png','pdf' = 'pdf'),inline = T),
                             downloadBttn('DownloadSharedPlot', 'Download')
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
                             checkboxInput('saveplot1',strong('click if save plot of result'),value = F),
                             checkboxInput('writetable1',strong('click if write table of result'),value = F),
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
                            radioButtons('phyloTreeType',h4('select tye of your phylotree'),
                                         c( 'NJtree(default)' = 'njtree',
                                            'Newick' = 'newick',
                                            'beast' = 'beast',
                                            'PAML' = 'PAML'),selected = 'njtree',inline = T),
                            conditionalPanel(
                              condition = "input.phyloTreeType != 'njtree'",
                              fileInput('phylotree.dir','upload your phylotree file')
                            ),
                            checkboxInput('use.indel','choose if  use indel',value = FALSE),
                            checkboxInput('show.mutSig','choose if show mutationignature',value = TRUE),
                            checkboxInput('show.heatmap','choose if show heatmap',value = TRUE),
                            textInput('heatmap.type','heatmap type','binary'),
                            checkboxInput('saveplot3','choose if saveplot',value = F),
                              actionBttn('submit.tree',div(
                              strong("Click ME to visualize result"),align = 'center',
                              icon("hand-right", lib = "glyphicon")))
                          )
                          ),
                          column(
                            width = 9,
                            box(
                            title = "Visulization",
                            width = NULL,
                            height = 800,
                            withSpinner(plotOutput("phylotree"))), 
                            box(
                              width = NULL,
                              status = 'success',
                              radioButtons('DownloadPhyloTreeCheck','Choose file type to download:',
                                         c('png' ='png','pdf' = 'pdf'),inline = T),
                              downloadBttn('DownloadPhyloTree', 'Download')
                            )
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
