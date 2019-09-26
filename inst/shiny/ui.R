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
    menuItem(strong("PhyloTree"), tabName = "Survival", icon = icon("line-chart")),
    menuItem(strong("Signature analysis"), tabName = "signature", icon = icon("bar-chart"))
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
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  title =  strong("Overview of MesKit package"),
                  # p("The typical workflow begins with MAF object creation by reading an MAF file combind with sample information. Based on Maf object, both ITH assessment and clonal analysis can be conducted. Furthermore, MesKit can perform function analysis and mutation analsysi on njtree object, which is converted from Maf object.",
                  #   style = "font-si16pt"),
                  # br(),
                  fluidRow(
                    column(
                      width = 7,
                      div(img(src = "images/pipeline.png", width=950,height = 700),style="text-align: left;")
                    ),
                    column(
                      width = 5,
                      br(),
                      p("The typical workflow begins with MAF object creation by reading an MAF file combind with sample information. Based on Maf object, both ITH assessment and clonal analysis can be conducted. Furthermore, MesKit can perform function analysis and mutation analsysi on njtree object, which is converted from Maf object.",
                        style = "font-si16pt"),
                      br(),
                      includeMarkdown("dom/Results_viewer.md")
                    )
                  )
                  )
                )
)

bodyIP <- tabItem("input",
                  h2('Input section'),
                  fluidRow(
                    column(
                      width = 4,
                      box(
                        title = div(shiny::icon("gear"), "Upload", inline =TRUE),
                        width = NULL,
                        radioGroupButtons(
                          "dataset",
                          label = "Load data",
                          choices = c(Example = "default", Upload = "upload"),
                          selected = "default",
                          status = "primary"
                        ),
                        conditionalPanel(
                          condition = "input.dataset == 'upload'",
                          fileInput('maf','Upload maf file'),
                          fileInput('sampleInfo','Upload sampleInfo' ),
                          conditionalPanel(
                            condition = "input.useccf == true",
                            fileInput('ccf.cluster','Upload ccf.cluster'),
                            fileInput('ccf.loci','Upload ccf.loci' )
                          )
                        ),
                        actionBttn('submit1',div(
                          strong("Click ME to start analysing"),align = 'center',
                          icon("thumbs-up", lib = "glyphicon")))
                      )
                    ),
                    column(
                      width = 4,
                      box(
                        width = NULL,
                        title = div(icon("book"), "Parameters"),
                        # selectInput('ref','select reference genome(hg19/hg38)',
                        #             choices = c('hg19','hg38'),selected = "hg19"),
                        switchInput('useccf','use ccf',value = F,
                                    onLabel = "TRUE",offLabel = "FALSE",labelWidth = 200),
                        switchInput('use.indel','use indel',value = FALSE,
                                    onLabel = "TRUE",offLabel = "FALSE",labelWidth = 200),
                        selectInput('ref','Select reference genome(hg19/hg38)',
                                    choices = c('hg19','hg38'),selected = "hg19")
                        # conditionalPanel(
                        #   condition = "input.useccf == true",
                        #   textInput("ccf.mutation.id","ccf.mutation.id",value = "Hugo_Symbol','Chromosome','Start_Position'"),
                        #   textInput("ccf.mutation.sep","ccf.mutation.sep",":")
                        #   )
                      )
                    )
                  )
                  )


bodyITH <- tabItem("ITH",
                   h2('ITH evaluation'),
                   fluidRow(
                     column(
                       width = 11,                     
                      box(
                       width = NULL,
                       tabBox(                           
                         height = "100%", 
                         width = "100%",
                         selected = "caInput02",
                         side = "left",
                         tabPanel(
                           title = div(icon("chart-bar"), "MathScore"),
                           value = "caInput02",
                           fluidRow(
                             column(
                               width = 3,
                               textInput("tsb",h4(strong("tsb")),value = "All",width = 300)
                             ),
                             column(
                               width = 5,
                               sliderInput("vafrange",h4(strong("VAF range")),min = 0,max = 1,value = c(0,1),width = 600)
                             )
                           ),
                           actionBttn('submit2',div(
                             strong("Click ME to start analysing"),align = 'center',
                             icon("hand-point-down", lib = "glyphicon"))),
                           br(),
                           br(),
                           withSpinner(DT::dataTableOutput('mathScore')),
                           br(),
                           downloadBttn('DownloadMathScore', 'Download')
                         ),
                         tabPanel(
                           title = div(icon("image"), "Vafplot"),
                           value = "caInput03",
                           fluidRow(
                             column(
                               width = 3,
                               selectInput("plotOption",h4(strong("Plot Option")),
                                           choices = c(
                                             compare = "compare",
                                             combine = "combine",
                                             separate = "separate"
                                           ), selected = "compare",width = 300)
                             ),
                             column(
                               width = 5,
                               sliderInput('width1','Adjust Image Width',min = 700,max = 1100, value = 850,width = 500)
                             )
                           ),
                           actionBttn('submit3',div(
                             strong("Click ME to start analysing"),align = 'center',
                             icon("hand-point-down", lib = "glyphicon"))),
                           br(),
                           br(),
                           plotOutput("vaf.cluster",height = "100%",width = "100%"),
                           br(),
                           fluidRow(
                             column(
                               width = 3,
                               radioButtons('DownloadVafPlotCheck','Choose file type to download:',
                                            c('png' ='png','pdf' = 'pdf'),inline = T)
                             ),
                             column(
                               width = 3,
                               downloadBttn('DownloadVafPlot', 'Download')
                             )
                           )
                         ),
                         tabPanel(
                           title = div(icon("map"), "Mutsharedprivateplot"),
                           value = "caInput04",
                           fluidRow(
                             column(
                               width = 5,
                               sliderInput('width2','Adjust Image Width',min = 700,max = 1100, value = 850,width = 500)
                             ),
                             column(
                               width = 5,
                               br(),
                               br(),
                               switchInput('show.num','Show Mutational Number',value = FALSE,
                                           onLabel = "TRUE",offLabel = "FALSE",labelWidth = 200)
                             )
                           ),
                           br(),
                           actionBttn('submit4',div(
                             strong("Click ME to start analysing"),align = 'center',
                             icon("hand-point-down", lib = "glyphicon"))),
                           br(),
                           br(),
                           plotOutput("mut.share_private",height = "100%"),
                           br(),
                           fluidRow(
                             column(
                               width = 3,
                               radioButtons('DownloadSharedPlotCheck','Choose file type to download:',
                                            c('png' ='png','pdf' = 'pdf'),inline = T)
                             ),
                             column(
                               width = 3,
                               downloadBttn('DownloadSharedPlot', 'Download')
                             )
                           )
                         ),
                         tabPanel(
                           title = div(icon("camara"), "Stackplot"),
                           value = "caInput05",
                           fluidRow(
                             column(
                               width = 3,
                               fileInput('oncogeneListFile','Upload oncogeneListFile',width = 300)
                             ),
                             column(
                               width = 3,
                               fileInput('tsgListFile','Upload tsgListFile',width = 300)
                             )
                           ),
                           fluidRow(
                             column(
                               width = 5,
                               sliderInput('width3','Adjust Image Width',min = 700,max = 1100, value = 850,width = 500)
                             ),
                             column(
                               width = 5,
                               br(),
                               br(),
                               actionBttn('submit5',div(
                                 strong("Click ME to start analysing"),align = 'center',
                                 icon("hand-point-down", lib = "glyphicon")))
                             )
                           ),
                           br(),
                           br(),
                           plotOutput("stackplot",height = "100%"),
                           br(),
                           fluidRow(
                             column(
                               width = 3,
                               radioButtons('DownloadStackPlotCheck','Choose file type to download:',
                                            c('png' ='png','pdf' = 'pdf'),inline = T)
                             ),
                             column(
                               width = 3,
                               downloadBttn('DownloadStackPlot', 'Download')
                             )
                           )
                         ),
                         tabPanel(
                           title = div(icon("box"), "Jaccardindex"),
                           value = "caInput06",
                           fluidRow(
                             column(
                               width = 5,
                               sliderInput('width4','Adjust Image Width',min = 700,max = 1100, value = 850,width = 500)
                             ),
                             column(
                               width = 5,
                               br(),
                               br(),
                               actionBttn('submit6',div(
                                 strong("Click ME to start analysing"),align = 'center',
                                 icon("hand-point-down", lib = "glyphicon")))
                             )
                           ),
                           br(),
                           br(),
                           plotOutput("JaccardIndex",height = "100%"),
                           br(),
                           fluidRow(
                             column(
                               width = 3,
                               radioButtons('DownloadJaccardIndexCheck','Choose file type to download:',
                                            c('png' ='png','pdf' = 'pdf'),inline = T)
                             ),
                             column(
                               width = 3,
                               downloadBttn('DownloadJaccardIndex', 'Download')
                             )
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
                       width = 10,
                      box(
                        width = NULL,
                        tabBox(
                          selected = 'c01',
                          side = 'left',
                          height = "100%",
                          width = "100%",
                          tabPanel(
                            value = 'c01',
                            title = div(icon("newspaper"), "TumorCloneplot"),
                            fluidRow(
                              column(
                                width = 5,
                                sliderInput('width5','Adjust image width',min = 700,max = 1100, value = 850,width = 500)
                              ),
                              column(
                                width = 5,
                                br(),
                                br(),
                                actionBttn('submit7',div(
                                  strong("Click ME to start analysing"),align = 'center',
                                  icon("hand-point-down", lib = "glyphicon")))
                              )
                            ),
                            # sliderInput('width5','Adjust image width',min = 700,max = 1100, value = 850,width = 500),
                            # actionBttn('submit7',div(
                            #   strong("Click ME to start analysing"),align = 'center',
                            #   icon("hand-point-down", lib = "glyphicon"))),
                            br(),
                            br(),
                            withSpinner(plotOutput('cloneplot',height = "100%")),
                            br(),
                            fluidRow(
                              column(
                                width = 3,
                                radioButtons('DownloadClonePlotCheck','Choose file type to download:',
                                             c('png' ='png','pdf' = 'pdf'),inline = T)
                              ),
                              column(
                                width = 3,
                                downloadBttn('DownloadClonePlot', 'Download')
                              )
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
                           width = 12,
                           box(
                            width = NULL,
                            height = 800,
                            tabBox(
                              side = 'left',
                              selected = 'F02',
                              width = "100%",
                              height = "100%",
                              tabPanel(
                                title = div(icon("lightbulb"), "GO analysis"),
                                value = 'F01',
                                fluidRow(
                                  column(
                                    width = 5,
                                    sliderInput('pval1','Adjust pval',min = 0,max = 1, value = 0.05),
                                    sliderInput('qval1','Adjust image qval',min = 0,max = 1, value = 0.2)
                                  ),
                                  column(
                                    width = 5,
                                    sliderInput('width6','Adjust image width',min = 400,max = 800, value = 800),
                                    sliderInput('height6','Adjust image height',min = 400,max = 600, value = 500)
                                  )
                                ),
                                fluidRow(
                                  column(
                                    width = 5,
                                    br(),
                                    br(),
                                    uiOutput("chooselist1"),
                                    actionBttn('submit8',div(
                                      strong("Click ME to start analysing"),align = 'center',
                                      icon("hand-right", lib = "glyphicon")))
                                  ),
                                  column(
                                    width = 7,
                                    withSpinner(plotOutput('GOplot',height = "100%",width = "100%")),
                                    br(),
                                    radioButtons('DownloadGOPlotCheck','Choose file type to download:',
                                                 c('png' ='png','pdf' = 'pdf'),inline = T),
                                    downloadBttn('DownloadGOPlot', 'Download')
                                  )
                                )
                              ),
                              tabPanel(
                                title = div(icon("microsoft"), "Pathway analysis"),
                                value = 'F02',
                                fluidRow(
                                  column(
                                    width = 5,
                                    sliderInput('pval2','Adjust pval',min = 0,max = 1, value = 0.05),
                                    sliderInput('qval2','Adjust image qval',min = 0,max = 1, value = 0.2)
                                  ),
                                  column(
                                    width = 5,
                                    sliderInput('width7','Adjust image width',min = 400,max = 800, value = 800),
                                    sliderInput('height7','Adjust image height',min = 400,max = 600, value = 500)
                                  )
                                ),
                                fluidRow(
                                  column(
                                    width = 5,
                                    br(),
                                    br(),
                                    uiOutput("chooselist2"),
                                    actionBttn('submit9',div(
                                      strong("Click ME to start analysing"),align = 'center',
                                      icon("hand-right", lib = "glyphicon")))
                                  ),
                                  column(
                                    width = 7,
                                    withSpinner(plotOutput('Pathwayplot',height = "100%",width = "100%")),
                                    br(),
                                    radioButtons('DownloadPathPlotCheck','Choose file type to download:',
                                                 c('png' ='png','pdf' = 'pdf'),inline = T),
                                    downloadBttn('DownloadPathPlotPlot', 'Download')
                                  )
                                )
                              )
                            )
                          )
                          )
                         

                        ))

bodySurvival <- tabItem('Survival',
                        h2('Phylotree visualiaztion'),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                            width = NULL,
                            status = 'primary',
                            title = div(shiny::icon("gear"), "Parameters", inline =TRUE),
                            selectInput('phyloTreeType','Select tye of your phylotree',
                                         c( 'njtree','newick','beast','PAML'),
                                         selected = 'njtree'),
                            conditionalPanel(
                              condition = "input.phyloTreeType != 'njtree'",
                              fileInput('phylotree.dir','upload your phylotree file')
                            ),
                            switchInput('show.mutSig','Choose if show mutationignature',value = TRUE,
                                        width = 500,onLabel = "TRUE",offLabel = "FALSE"),
                            switchInput('show.heatmap','Choose if show heatmap',value = TRUE,width = 500,onLabel = "TRUE",offLabel = "FALSE"),
                            radioGroupButtons(
                              inputId = "sig.name",
                              label = "Signature name",
                              choices = c(default = "default", alias = "alias"),
                              selected = "default",
                              status = "primary"
                            ),
                            radioGroupButtons(
                              inputId = "heatmap.type",
                              label = "Heatmap type",
                              choices = c(binary = "binary", CCF = "CCF"),
                              selected = "binary",
                              status = "primary"
                            ),
                            actionBttn('submit10',div(
                              strong("Click ME to start analysing"),align = 'center',
                              icon("hand-right", lib = "glyphicon")))
                          )
                          ),
                          column(
                            width = 9,
                            box(
                              title = div(icon("eye"), "Visulization"),
                              width = NULL,
                              height = 1000,
                              withSpinner(plotOutput("phylotree",height = 700)),
                              br(),
                              fluidRow(
                                column(
                                  width = 3,
                                  radioButtons('DownloadPhyloTreeCheck','Choose file type to download:',
                                               c('png' ='png','pdf' = 'pdf'),inline = T)
                                ),
                                column(
                                  width = 3,
                                  downloadBttn('DownloadPhyloTree', 'Download')
                                )
                              )
                                )
                            )
                            
                          
                        )
                        )
bodysignature <- tabItem('signature',
                         h2('Signature analysis'),
                         fluidRow(
                           column(
                             width = 10,
                             box(
                               width = NULL,
                               tabBox(
                                 side = 'left',
                                 selected = 'sa1',
                                 width = "100%",
                                 height = "100%",
                                 tabPanel(
                                   title = div(icon("bar-chart"), "Signature"),
                                   value = "sa1",
                                   actionBttn('submit11',div(
                                     strong("Click ME to start analysing"),align = 'center',
                                     icon("hand-point-down", lib = "glyphicon"))),
                                   br(),
                                   br(),
                                   withSpinner(plotOutput("signature",height = 700)),
                                   br(),
                                   fluidRow(
                                     column(
                                       width = 3,
                                       radioButtons('DownloadSignaturePlotCheck','Choose file type to download:',
                                                    c('png' ='png','pdf' = 'pdf'),inline = T)
                                     ),
                                     column(
                                       width = 3,
                                       downloadBttn('DownloadSignaturePlot', 'Download')
                                     )
                                   )
                                   # box(
                                   #   width = NULL,
                                   #   status = 'success',
                                   #   radioButtons('DownloadSignaturePlotCheck','Choose file type to download:',
                                   #                c('png' ='png','pdf' = 'pdf'),inline = T),
                                   #   downloadBttn('DownloadSignaturePlot', 'Download')
                                   # )
                                 )
                               )
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
      tags$head(
        tags$style(HTML("
                        .shiny-output-error-validation {
                            color: brown;
                        }
                        # .pipediv{
                        #     width:900px;
                        #     height:500px;
                        # }
                        # .pipediv .pipe{
                        #     float:left;
                        # }
                        # .pipediv .pipe img{
                        #     width:500px;
                        #     height:500px;
                        # }
                        # .pipetext{
                        #     
                        # }
                        ")),
        tags$link(rel = "stylesheet", type = "text/css", href = "css/main.css")
      ),
      tabItems(
        bodyHome,
        bodyIP,
        bodyITH,
        bodyclone,
        bodyfunction,
        bodySurvival,
        bodysignature
      )
    )
  )
)
