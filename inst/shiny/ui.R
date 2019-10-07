options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(deconstructSigs))
suppressMessages(library(DT))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))
suppressMessages(library(shinyjs))

#sider bar----

sidebar <- dashboardSidebar(
  width = 300,
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
                      width = 12, 
                      column(
                        width = 3,
                        box(
                          title = div(shiny::icon("gear"), "Upload", inline =TRUE),
                          width = NULL,
                          conditionalPanel(
                            condition = "true",
                            fileInput(inputId = 'maf', 
                                      label = div(style = "font-size:18px; font-weight:600; ", 'MAF file'), 
                                      placeholder = "example data: 311252.maf", 
                                      width = 400),
                            fileInput(inputId = 'sampleInfo', 
                                      label = div(style = "font-size:18px; font-weight:600; ", 'Sample information document'), 
                                      placeholder = "example data: sample_info.txt", 
                                      width = 400),
                            checkboxInput(inputId = 'useccf', label = div(style = "font-size:15px; ", 'use ccf'),value = FALSE, width = 200),
                            conditionalPanel(
                              condition = "input.useccf == true",
                              fileInput('ccf.cluster',label = div(style = "font-size:18px; font-weight:600; ", 'ccf.cluster')),
                              fileInput('ccf.loci',label = div(style = "font-size:18px; font-weight:600; ", 'ccf.loci'))
                            ),
                            checkboxInput('use.indel', label = div(style = "font-size:15px; ", 'use indel'),value = FALSE,width = 400),
                            selectInput('ref', label = div(style = "font-size:18px; font-weight:600; ", 'Select reference genome(hg19/hg38)'),
                                        choices = c('hg19','hg38'),selected = "hg19", width = 400)
                          ),
                          actionBttn('submit1',div(
                            strong("Click ME to start analysing"),align = 'center',
                            icon("thumbs-up", lib = "glyphicon")))
                        )
                      ), 
                      column(
                        width = 9, 
                        withSpinner(DT::dataTableOutput('maftable', width = '100%')),
                        uiOutput("mafdb")
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
                             h3(strong("Parameter: ")), 
                             value = "caInput02",
                             fluidRow(
                               column(
                                 width = 3,
                                 tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:18px; font-weight:600; ", "Min VAF:")),
                                           tags$td(width = "70%", textInput(inputId = "minvaf", value = 0.00, label = NULL)))
                                 ), 
                                 br(),
                                 tags$table(
                                   tags$tr(id = "inline",
                                           width = "100%",
                                           tags$td(width = "30%", tags$div(style = "font-size:18px; font-weight:600; ", "Max VAF:")),
                                           tags$td(width = "70%", textInput(inputId = "maxvaf", value = 1.00, label = NULL)))
                                 ), 
                                 br(),
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
                             h3(strong("Parameter: ")), 
                             value = "caInput03",
                             fluidRow(
                               column(
                                 width = 3,
                                 selectInput("plotOption", label = div(style = "font-size:18px; font-weight:600;  ", "Plot option"),
                                             choices = c(
                                               Compare = "compare",
                                               Combine = "combine",
                                               Separate = "separate"
                                             ), selected = "compare",width = 300),
                                 selectInput("themeOption", label = div(style = "font-size:18px; font-weight:600;  ", "Theme option"),
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
                                 sliderInput('width1', label = div(style = "font-size:18px; font-weight:600; ", 'Image width'), min = 700,max = 1100, value = 850,width = 500),
                                 br(),
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
                             h3(strong("Parameter: ")), 
                             value = "caInput04",
                             fluidRow(
                               column(
                                 width = 3,
                                 br(),
                                 checkboxInput('show.num',label = div(style = "font-size:15px; font-weight:400; ", 'Show mutation number'),width = 200),
                                 br(),
                                 sliderInput('width2', label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
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
                             h3(strong("Parameter: ")), 
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
                                 sliderInput('width3',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 600,max = 1100, value = 650,width = 500),
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
                             h3(strong("Parameter: ")), 
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
                                 sliderInput('width4',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
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
                         width = 12,
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
                               h3(strong("Parameter: ")), 
                               fluidRow(
                                 column(
                                   width = 3,
                                   sliderInput('width5',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
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
                                  h3(strong("Parameter: ")), 
                                  value = 'F01',
                                  fluidRow(
                                    column(
                                      width = 3,
                                      tags$table(
                                        tags$tr(id = "inline", 
                                                width = "100%",
                                                tags$td(width = "20%", div(style = "font-size:18px; font-weight:600;  ", "Pval:")),
                                                tags$td(width = "70%", textInput(inputId = "pval1", value = 0.05, label = NULL)))
                                      ), 
                                      br(),
                                      tags$table(
                                        tags$tr(id = "inline",
                                                width = "100%",
                                                tags$td(width = "20%", tags$div(style = "font-size:18px; font-weight:600; ", "Qval:")),
                                                tags$td(width = "70%", textInput(inputId = "qval1", value =  0.20, label = NULL)))
                                      ), 
                                      br(),
                                      sliderInput('width6',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800),
                                      sliderInput('height6',label = div(style = "font-size:18px; font-weight:600; ", 'Image height'),min = 400,max = 600, value = 500),
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
                                  h3(strong("Parameter: ")), 
                                  value = 'F02',
                                  fluidRow(
                                    column(
                                      width = 3,
                                      tags$table(
                                        tags$tr(id = "inline", 
                                                width = "100%",
                                                tags$td(width = "20%", div(style = "font-size:18px; font-weight:600; ", "Pval:")),
                                                tags$td(width = "70%", textInput(inputId = "pval2", value = 0.05, label = NULL)))
                                      ), 
                                      br(),
                                      tags$table(
                                        tags$tr(id = "inline",
                                                width = "100%",
                                                tags$td(width = "20%", tags$div(style = "font-size:18px; font-weight:600; ", "Qval:")),
                                                tags$td(width = "70%", textInput(inputId = "qval2", value =  0.20, label = NULL)))
                                      ), 
                                      br(),
                                      sliderInput('width7',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800),
                                      sliderInput('height7',label = div(style = "font-size:18px; font-weight:600; ", 'Image height'),min = 400,max = 600, value = 500),
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
                                 tabPanel(title = div(icon("lightbulb"), "Mutational signature summary"), 
                                          h3(strong("Parameter: ")), 
                                          value = 'S01',
                                          fluidRow(
                                            column(
                                              width = 3,
                                              fileInput('driverGenesFile',label = div(style = "font-size:18px; font-weight:600; ", 'Upload driverGenesFile')), 
                                              numericInput('mutThreshold', div(style = "font-size:18px; font-weight:600;  ", 'Mutation quantity threshold'), value = 50),
                                              selectInput("signaturesRef", label = div(style = "font-size:18px; font-weight:600;  ", "Signautre reference"),
                                                          choices = c(signatures.cosmic = "signatures.cosmic",
                                                                      signatures.nature2013 = "signatures.nature2013"),
                                                          selected = "signatures.cosmic"),
                                              br(),
                                              br(),
                                              actionBttn('submitSig',div(
                                                strong("Click ME to start analysing"),align = 'center',
                                                icon("hand-right", lib = "glyphicon")))
                                            ), 
                                            column(
                                              width = 9, 
                                              withSpinner(DT::dataTableOutput('sigOFA'))
                                            )
                                            )
                                 ), 
                                 tabPanel(title = div(icon("lightbulb"), "Mutational signature analysis"), 
                                          h3(strong("Parameter: ")), 
                                          value = 'S02',
                                          fluidRow(
                                            column(
                                              width = 3,
                                              fileInput('driverGenesFile2',label = div(style = "font-size:18px; font-weight:600; ", 'Upload driverGenesFile')), 
                                              numericInput('mutThreshold2', div(style = "font-size:18px; font-weight:600;  ", 'Mutation quantity threshold'), value = 50, step=10),
                                              selectInput("signaturesRef2", label = div(style = "font-size:18px; font-weight:600;  ", "Signautre reference"),
                                                          choices = c(signatures.cosmic = "signatures.cosmic",
                                                                      signatures.nature2013 = "signatures.nature2013"),
                                                          selected = "signatures.cosmic"),
                                              radioButtons(
                                                inputId = "sigplot", 
                                                label = div(style = "font-size:18px; font-weight:600; ", "Plot option"), 
                                                choiceNames = list(
                                                  tags$span(style = "font-size:14.5px; font-weight:400; ", "Signature probability"), 
                                                  tags$span(style = "font-size:14.5px; font-weight:400; ", "Branch-trunk")
                                                ),
                                                choiceValues = c("signaturesprob", "branchtrunk"),
                                                selected = "signaturesprob", 
                                                inline = TRUE), 
                                                conditionalPanel(
                                                  condition = "input.sigplot == 'branchtrunk'",
                                                  numericInput('signiflevel', div(style = "font-size:18px; font-weight:600;  ", 'Significant level'), value = 0.05, min=0, max=1, step=0.1)
                                                ), 
                                                sliderInput(inputId='widthsig2',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800, width = 500),
                                                sliderInput(inputId='heightsig2',label = div(style = "font-size:18px; font-weight:600; ", 'Image height'),min = 400,max = 600, value = 500, width = 500), 
                                              br(),
                                              br(),
                                              actionBttn('submitSig2',div(
                                                strong("Click ME to start analysing"),align = 'center',
                                                icon("hand-right", lib = "glyphicon")))
                                            ), column(
                                              width = 9,
                                              withSpinner(plotOutput('sigOFA2', height = "100%", width = "100%")),
                                              br(),
                                              uiOutput("sigpdb")
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
                                  h3(strong("Parameter: ")), 
                                  width = 3,
                                  selectInput('phyloTreeType',div(style = "font-size:18px; font-weight:600; ", 'Type'),
                                              c( 'njtree','newick','beast','PAML'),
                                              selected = 'njtree'),
                                  conditionalPanel(
                                    condition = "input.phyloTreeType != 'njtree'",
                                    fileInput('phylotree.dir','upload your phylotree file')
                                  ),
                                  checkboxInput('show.mutSig',div(style = "font-size:15px; font-weight:400; ", 'Show mutation signature'),value = TRUE),
                                  checkboxInput('show.heatmap',div(style = "font-size:15px; font-weight:400; ", 'Show heatmap'),value = TRUE),
                                  radioButtons(
                                    inputId = "sig.name", 
                                    label = div(style = "font-size:18px; font-weight:600; ", "Signature name"), 
                                    choiceNames = list(
                                      tags$span(style = "font-size:14.5px; font-weight:400; ", "Default"), 
                                      tags$span(style = "font-size:14.5px; font-weight:400; ", "Alias")
                                    ),
                                    choiceValues = c("default", "alias"),
                                    selected = "default", 
                                    inline = TRUE),
                                  radioButtons(
                                    inputId = "heatmap.type",
                                    label = div(style = "font-size:18px; font-weight:600; ", "Heatmap type"),
                                    choiceNames = list(
                                      tags$span(style = "font-size:14.5px; font-weight:400; ", "Binary"), 
                                      tags$span(style = "font-size:14.5px; font-weight:400; ", "CCF")
                                    ),
                                    choiceValues = c("binary", "CCF"),
                                    selected = "binary", 
                                    inline = TRUE
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


#Main function----
dbHeader <- dashboardHeader(title = "Meskit", titleWidth = 300, 
                tags$li(class = "dropdown", actionLink(inputId = "help", label = div(style = "font-size:15px; font-weight:400; ", "Help"))), 
                tags$li(class = "dropdown", actionLink(inputId = "contact", label = div(style = "font-size:15px; font-weight:400; ", "Contact"))))
dbHeader$children[[2]]$children <-  tags$a(href='https://github.com/Niinleslie/MesKit',
                                           tags$img(src='image/logo.png',height='65',width='250'))
shinyUI(
  dashboardPage(
    skin = "blue",
    header=dbHeader ,
    sidebar=sidebar,
    body=dashboardBody(
      ## add text behind the sidebar (design error)
      tags$head(tags$style(HTML(
        '/* logo */
        .skin-blue .main-header .logo {
                              background-color: #3c8dbc;
        }
        /* logo when hovered */
        .skin-blue .main-header .logo:hover {
                              background-color: #3c8dbc;
                              }
        .textnvbar { 
        font-size: 20px;
        line-height: 50px;
        text-align: left;
        font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
        padding: 0 15px;
        overflow: hidden;
        color: white;
        }'))),
      
      ## change the style of the progress bar.
      tags$head(tags$style(
        type="text/css", 
        ".progress.shiny-file-input-progress {
            height:3.5px;
        }
          
        .progress-bar {
          background-image: linear-gradient(to right, #77C7FF, #3c8dbc ) !important;
          background-size: auto !important;
          font-size:0px;
          height:3.5px;
        }"
      )),
      tags$script(HTML('
      $(document).ready(function() {
        $("header").find("nav").append(\'<span class="textnvbar"> Meskit: Analysis and visualize multi-sample whole-exome sequencing data</span>\');
      })
     ')), 
      ## used for inline setting
      tags$head(
        tags$style(type="text/css", "#inline label{ display: table-cell; text-align: centers; vertical-align: middle; width=400; } 
                #inline .form-group { display: table-row; width=400; }")
      ),
      tags$head(
        tags$style(HTML(".shiny-output-error-validation {
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
