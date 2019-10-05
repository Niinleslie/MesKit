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
      menuItem(strong("Signature analysis"), tabName = "signature", icon = icon("bar-chart")), 
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
                            label = "",
                            choices = c(Example = "default", Upload = "upload"),
                            selected = "default",
                            status = "primary"
                          ),
                          br(),
                          conditionalPanel(
                            condition = "input.dataset == 'upload'",
                            fileInput('maf','maf file'),
                            fileInput('sampleInfo','sampleInfo' ),
                            checkboxInput('useccf',strong('use ccf'),value = F,width = 200),
                            conditionalPanel(
                              condition = "input.useccf == true",
                              fileInput('ccf.cluster','ccf.cluster'),
                              fileInput('ccf.loci','ccf.loci' )
                            ),
                            checkboxInput('use.indel',strong('use indel'),value = FALSE,width = 400),
                            selectInput('ref','Select reference genome(hg19/hg38)',
                                        choices = c('hg19','hg38'),selected = "hg19")
                          ),
                          actionBttn('submit1',div(
                            strong("Click ME to start analysing"),align = 'center',
                            icon("thumbs-up", lib = "glyphicon")))
                        )
                      )
                    )
                    )
  
  
  bodyITH <- tabItem("ITH",
                     h2('ITH evaluation'),
                     fluidRow(
                       column(
                         width = 12,                     
                        box(
                         width = NULL,
                         tabBox(                           
                           height = "100%", 
                           width = "100%",
                           selected = "caInput02",
                           side = "left",
                           tabPanel(
                             title = div(icon("chart-bar"), "MathScore"),
                             h3("Parameter: "), 
                             value = "caInput02",
                             fluidRow(
                               column(
                                 width = 3,
                                 textInput("tsb", 
                                           h4(strong("Tumor sample barcode")), 
                                           value = "All",width = 300),
                                 div(
                                   style="display:inline-block", 
                                   textInput(inputId="minvaf", 
                                             label=h4(strong("minvaf")), 
                                             value = 0.0, width = 150)),
                                 div(style="display:inline-block", 
                                     textInput(inputId="maxvaf", 
                                               label=h4(strong("maxvaf")), 
                                               value = 1.0, width = 150)),
                                 br(),
                                 actionBttn('submit2',div(
                                   strong("Click ME to start analysing"),align = 'center',
                                   icon("hand-point-down", lib = "glyphicon")))
                                 ),
                               column(
                                 width = 8,
                                 withSpinner(DT::dataTableOutput('mathScore')),
                                 br(),
                                 br(),
                                 uiOutput("msdb")
                               )
                             )
                           ),
                           tabPanel(
                             title = div(icon("image"), "Vafplot"),
                             h3("Parameter: "), 
                             value = "caInput03",
                             fluidRow(
                               column(
                                 width = 3,
                                 selectInput("plotOption",h4(strong("Plot option")),
                                             choices = c(
                                               Compare = "compare",
                                               Combine = "combine",
                                               Separate = "separate"
                                             ), selected = "compare",width = 300),
                                 selectInput("themeOption",h4(strong("Theme option")),
                                             choices = c(NPG = "npg",
                                                         AAAS = "aaas",
                                                         NEJM = "nejm",
                                                         Lancet = "lancet",
                                                         JAMA = "jama",
                                                         JCO = "jco",
                                                         UCSCGB = "ucscgb",
                                                         D3 = "d3",
                                                         LocusZoom = "locuszoom",
                                                         IGV = "igv",
                                                         UChicago = "uchicago",
                                                         'Star Trek' = "startrek",
                                                         'Tron Legacy' = 'tron',
                                                         Futurama = "futurama",
                                                         'Rick and Morty' = 'rickandmorty',
                                                         'The Simpsons' = 'simpsons',
                                                         GSEA = 'gsea'),
                                             selected = "aaas",width = 300),
                                 sliderInput('width1','Adjust Image Width',min = 700,max = 1100, value = 850,width = 500),
                                 br(),
                                 actionBttn('submit3', div(
                                   strong("Click ME to start analysing"),align = 'center',
                                   icon("hand-point-down", lib = "glyphicon")),
                                   )
                               ),
                               column(
                                 width = 8,
                                 withSpinner(plotOutput("vaf",height = "100%")), 
                                 uiOutput("vcdb")
                               )
                             )
                           ),
                           tabPanel(
                             title = div(icon("map"), "Mutsharedprivateplot"),
                             h3("Parameter: "), 
                             value = "caInput04",
                             fluidRow(
                               column(
                                 width = 3,
                                 br(),
                                 checkboxInput('show.num',strong('Show mutation number'),width = 200),
                                 br(),
                                 sliderInput('width2','Adjust image width',min = 700,max = 1100, value = 850,width = 500),
                                 br(),
                                 actionBttn('submit4',div(
                                   strong("Click ME to start analysing"),align = 'center',
                                   icon("hand-point-down", lib = "glyphicon")))
                               ),
                               column(
                                 width = 8,
                                 withSpinner(plotOutput("mut.share_private",height = "100%")),   
                                 uiOutput("mspdb")
                               )
                             )
                           ),
                           tabPanel(
                             title = div(icon("camara"), "Stackplot"),
                             h3("Parameter: "), 
                             value = "caInput05",
                             fluidRow(
                               column(
                                 width = 3,
                                 fileInput('oncogeneListFile','Upload oncogeneListFile',width = 300),
                                 fileInput('tsgListFile','Upload tsgListFile',width = 300),
                                 selectInput("themeOption2","Theme option",
                                             choices = c(NPG = "npg",
                                                         AAAS = "aaas",
                                                         NEJM = "nejm",
                                                         Lancet = "lancet",
                                                         JAMA = "jama",
                                                         JCO = "jco",
                                                         UCSCGB = "ucscgb",
                                                         D3 = "d3",
                                                         LocusZoom = "locuszoom",
                                                         IGV = "igv",
                                                         UChicago = "uchicago",
                                                         'Star Trek' = "startrek",
                                                         'Tron Legacy' = 'tron',
                                                         Futurama = "futurama",
                                                         'Rick and Morty' = 'rickandmorty',
                                                         'The Simpsons' = 'simpsons',
                                                         GSEA = 'gsea'),
                                             selected = "aaas",width = 300),
                                 br(),
                                 checkboxInput('show.percentage',label = strong('Show percentage'),value = T),
                                 br(),
                                 sliderInput('width3','Adjust image width',min = 600,max = 1100, value = 700,width = 500),
                                 br(),
                                 actionBttn('submit5',div(
                                   strong("Click ME to start analysing"),align = 'center',
                                   icon("hand-point-down", lib = "glyphicon")))
                               ),
                               column(
                                 width = 8,
                                 withSpinner(plotOutput("stackplot",height = "100%")),
                                 br(),
                                 uiOutput("stkdb")
                               )
                             )
                           ),
                           tabPanel(
                             title = div(icon("box"), "Jaccardindex"),
                             h3("Parameter: "), 
                             value = "caInput06",
                             fluidRow(
                               column(
                                 width = 3,
                                 selectInput("JItype",h4(strong("Type")),
                                             choices = c(
                                               Lower = "lower",
                                               Upper = "upper",
                                               Full = "full"
                                             ), selected = "lower",width = 300),
                                 sliderInput('width4','Adjust image width',min = 700,max = 1100, value = 850,width = 500),
                                 br(),
                                 br(),
                                 actionBttn('submit6',div(
                                   strong("Click ME to start analysing"),align = 'center',
                                   icon("hand-point-down", lib = "glyphicon")))
                               ),
                               column(
                                 width = 8,
                                 withSpinner(plotOutput("JaccardIndex",height = "100%")) ,
                                 uiOutput("jidb")
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
                         width = 11,
                        box(
                          width = NULL,
                          tabBox(
                            selected = 'c01',
                            side = 'left',
                            height = "100%",
                            width = "100%",
                            tabPanel(
                              value = 'c01',
                              title = div(icon("newspaper"), "Tumorcloneplot"),
                              h3("Parameter: "), 
                              fluidRow(
                                column(
                                  width = 3,
                                  sliderInput('width5','Adjust image width',min = 700,max = 1100, value = 850,width = 500),
                                  br(),
                                  br(),
                                  actionBttn('submit7',div(
                                    strong("Click ME to start analysing"),align = 'center',
                                    icon("hand-point-down", lib = "glyphicon")))
                                ),
                                column(
                                  width = 8,
                                  withSpinner(plotOutput('cloneplot',height = "100%")),
                                  uiOutput("clpdb")
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
                                  h3("Parameter: "), 
                                  value = 'F01',
                                  fluidRow(
                                    column(
                                      width = 3,
                                      div(
                                        style="display:inline-block", 
                                        textInput(inputId="pval1", 
                                                  label="Pval", 
                                                  value = 0.05, width = 150)),
                                      div(style="display:inline-block", 
                                          textInput(inputId="qval1", 
                                                    label="Qval", 
                                                    value = 0.2, width = 150)),
                                      sliderInput('width6','Adjust image width',min = 400,max = 1000, value = 800),
                                      sliderInput('height6','Adjust image height',min = 400,max = 600, value = 500),
                                      br(),
                                      br(),
                                      actionBttn('submit8',div(
                                        strong("Click ME to start analysing"),align = 'center',
                                        icon("hand-right", lib = "glyphicon")))
                                    ),
                                    column(
                                      width = 9,
                                      uiOutput("chooselist1"),
                                      withSpinner(plotOutput('GOplot',height = "100%",width = "100%")),
                                      uiOutput("GOdb")
                                    )
                                  )
                                ),
                                tabPanel(
                                  title = div(icon("microsoft"), "Pathway analysis"),
                                  h3("Parameter: "), 
                                  value = 'F02',
                                  fluidRow(
                                    column(
                                      width = 3,
                                      div(
                                        style="display:inline-block", 
                                        textInput(inputId="pval2", 
                                                  label="Pval", 
                                                  value = 0.05, width = 150)),
                                      div(style="display:inline-block", 
                                          textInput(inputId="qval2", 
                                                    label="Qval", 
                                                    value = 0.2, width = 150)),
                                      sliderInput('width7','Adjust image width',min = 400,max = 1000, value = 800),
                                      sliderInput('height7','Adjust image height',min = 400,max = 600, value = 500),
                                      br(),
                                      actionBttn('submit9',div(
                                        strong("Click ME to start analysing"),align = 'center',
                                        icon("hand-right", lib = "glyphicon")))
                                    ),
                                    column(
                                      width = 9,
                                      uiOutput("chooselist2"),
                                      withSpinner(plotOutput('Pathwayplot',height = "100%",width = "100%")),
                                      uiOutput("Pathdb")
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
                              width = 12,
                              box(
                                width = NULL,
                                fluidRow(
                                  column(
                                    h3("Parameter: "), 
                                    width = 3,
                                    selectInput('phyloTreeType',h4(strong('Type')),
                                                c( 'njtree','newick','beast','PAML'),
                                                selected = 'njtree'),
                                    conditionalPanel(
                                      condition = "input.phyloTreeType != 'njtree'",
                                      fileInput('phylotree.dir','upload your phylotree file')
                                    ),
                                    checkboxInput('show.mutSig',strong('Show mutation signature'),value = TRUE),
                                    checkboxInput('show.heatmap',strong('Show heatmap'),value = TRUE),
                                    radioGroupButtons(
                                      inputId = "sig.name",
                                      label = "Signature name",
                                      choices = c(Default = "default", Alias = "alias"),
                                      selected = "default",
                                      status = "primary"
                                    ),
                                    radioGroupButtons(
                                      inputId = "heatmap.type",
                                      label = "Heatmap type",
                                      choices = c(Binary = "binary", CCF = "CCF"),
                                      selected = "binary",
                                      status = "primary"
                                    ),
                                    actionBttn('submit10',div(
                                      strong("Click ME to start analysing"),align = 'center',
                                      icon("hand-right", lib = "glyphicon")))
                                  ),
                                  column(
                                    width = 9,
                                    withSpinner(plotOutput("phylotree",height = 700)),
                                    br(),
                                    uiOutput("phtdb")
                                  )
                                )
                              )
                            )
                          )
                          )
  
  bodySignature <- tabItem('signature',
                           h2('Mutational signature analysis'),
                           fluidRow(
                             column(
                               width = 12,
                               box(
                                 width = NULL,
                                 height = 800,
                                 tabBox(
                                   side = 'left',
                                   selected = 'S01',
                                   width = "100%",
                                   height = "100%",
                                   tabPanel(title = div(icon("lightbulb"), "Mutational signature data"), 
                                            h3("Parameter: "), 
                                            value = 'S01',
                                            fluidRow(
                                              column(
                                                width = 3,
                                                textInput('pval1','Pval',value = 0.05),
                                                textInput('qval1','Qval',value = 0.2),
                                                sliderInput('width6','Adjust image width',min = 400,max = 1000, value = 800),
                                                sliderInput('height6','Adjust image height',min = 400,max = 600, value = 500),
                                                br(),
                                                br(),
                                                actionBttn('submit8',div(
                                                  strong("Click ME to start analysing"),align = 'center',
                                                  icon("hand-right", lib = "glyphicon")))
                                              ))
                                            ), 
                                   tabPanel("Mutation probability plot",
                                            h3("Parameter: ")
                                            ), 
                                   tabPanel("Trunk-branch plot", 
                                            h3("Parameter: ")
                                            )
                                   
                                   
                                   
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
          bodySignature,
          bodySurvival
        )
      )
    )
  )
