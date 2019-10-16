options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(deconstructSigs))
suppressMessages(library(DT))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))
suppressMessages(library(shinyjs))
suppressMessages(library(Meskit))


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
                        h3(strong("Introduction")),
                        p("Malignant tumor, is considered one of the most serious threats to human health. Among numerous causes attributable to tumor, the intra-tumoral heterogeneity (ITH), once ignored, is now thought to have been a key factor contributing to the therapeutic failure. Today, with the rapid development of high-throughput sequencing technologies, sequencing of spatially or temporally distinct tumor regions has begun to uncover the bewildering extent of diversity within tumors. To facilitate the rapid analysis, we present an R package, MesKit, providing comprehensive analysis that are commonly used in cancer genomic ITH studies, including ITH evaluation, enrichment, signature, clone evolution analysis, also allowing visualizing phylogenetic trees.",
                          style = "font-size:18px; font-weight:500;line-height:40px;"),
                        br()
                      )
                    ),
                    
                    fluidRow(
                      box(
                        width = 12,
                        status = "info",
                        solidHeader = TRUE,
                        title =  strong("Overview of MesKit package"),
                        fluidRow(
                          column(
                            width = 7,
                            div(img(src = "images/pipeline.png", width=950,height = 500),style="text-align: left;")
                          ),
                          column(
                            width = 5,
                            h3(strong("With this MesKit's Shiny APP:")),
                            p("- Evaluate the intra-tumoral heterogeneity (ITH) .",br(),
                              " - Perform clonal analysis.",br(),
                              "- Perform enrichment analysis of GO ontology and pathway.",br(),
                              "- Perform mutational signature analysis.",br(),
                              "- Construct phylogenetic trees.",
                              style = "font-size:18px; font-weight:500;line-height:50px;")
                            # br(),
                            # p("The typical workflow begins with MAF object creation by reading an MAF file combind with sample information. Based on Maf object, both ITH assessment and clonal analysis can be conducted. Furthermore, MesKit can perform function analysis and mutation analsysi on njtree object, which is converted from Maf object.",
                            #   style = "font-si16pt"),
                            # br(),
                            # includeMarkdown("dom/Results_viewer.md")
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
                                      label = div(style = "font-size:18px; font-weight:600;",'MAF file',
                                                  tags$button(
                                                    Id = "iecontrol01",
                                                    type = "button",
                                                    class = "bttn-material-circle",
                                                    class = "btn action-button",
                                                    list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                    style = " background-position: center;padding:0;margin-bottom:7px;"
                                                  )
                                                 ), 
                                      placeholder = "example data: 311252.maf", 
                                      width = 400),
                            fileInput(inputId = 'sampleInfo', 
                                      label = div(style = "font-size:18px; font-weight:600; ", 'Sample information document',
                                                  tags$button(
                                                    Id = "iecontrol02",
                                                    type = "button",
                                                    class = "bttn-material-circle",
                                                    class = "btn action-button",
                                                    list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                    style = " background-position: center;padding:0;margin-bottom:7px;"
                                                  )
                                                  ), 
                                      placeholder = "example data: sample_info.txt", 
                                      width = 400),
                            checkboxInput(inputId = 'useccf', label = div(style = "font-size:15px; ", 'use ccf'),value = FALSE, width = 200),
                            conditionalPanel(
                              condition = "input.useccf == true",
                              fileInput('ccf.cluster',label = div(style = "font-size:18px; font-weight:600; ", 'ccf.cluster file',
                                                                  tags$button(
                                                                    Id = "iecontrol03",
                                                                    type = "button",
                                                                    class = "bttn-material-circle",
                                                                    class = "btn action-button",
                                                                    list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                                    style = " background-position: center;padding:0;margin-bottom:7px;"
                                                                  )
                                                                  ),
                                        placeholder = "example data: 311252.cluster.tsv", width = 400),
                              fileInput('ccf.loci',label = div(style = "font-size:18px; font-weight:600; ", 'ccf.loci file',
                                                               tags$button(
                                                                 Id = "iecontrol04",
                                                                 type = "button",
                                                                 class = "bttn-material-circle",
                                                                 class = "btn action-button",
                                                                 list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                                 style = " background-position: center;padding:0;margin-bottom:7px;"
                                                               )
                                                               ),
                                        placeholder = "example data: 311252.loci.tsv", width = 400)
                            ),
                            checkboxInput('use.indel', label = div(style = "font-size:15px; ", 'use indel'),value = FALSE,width = 400),
                            selectInput('ref', label = div(style = "font-size:18px; font-weight:600; ", 'Select reference genome(hg19/hg38)'),
                                        choices = c('hg19','hg38'),selected = "hg19", width = 400)
                          ),
                          actionBttn('submit1',div(
                            strong("Start analysing"),align = 'center',
                            icon("hand-right", lib = "glyphicon")))
                          # uiOutput("pb1")
                          # progressBar(
                          #   id = "pb1",
                          #   value = 0,
                          #   total = 100,
                          #   title = "Wait",
                          #   display_pct = FALSE, 
                          #   status = "custom"
                          # )
                        )
                      ), 
                      column(
                        width = 9, 
                        # box(
                        #   width = NULL,
                        #   withSpinner(DT::dataTableOutput('maftable', width = '100%')),
                        #   uiOutput("ie1"),
                        #   br(),
                        #   uiOutput("ie2")
                        # )
                        uiOutput("ie1"),
                        uiOutput("ie2"),
                        uiOutput("ie3"),
                        uiOutput("ie4")
                      )
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
                         conditionalPanel(
                           condition = "input.tith == 'caInput02'",
                           h3(strong("Parameter ")), 
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
                             strong("Start analysing"),align = 'center',
                             icon("hand-right", lib = "glyphicon")))
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput03'",
                           h3(strong("Parameter ")), 
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
                             strong("Start analysing"),align = 'center',
                             icon("hand-right", lib = "glyphicon"))
                           )
                           # progressBar(
                           #   id = "pb3",
                           #   value = 0,
                           #   total = 100,
                           #   title = "",
                           #   display_pct = TRUE, 
                           #   status = "custom"
                           # )
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput04'",
                           h3(strong("Parameter ")), 
                           br(),
                           checkboxInput('show.num1',label = div(style = "font-size:15px; font-weight:400; ", 'Show mutation number'),width = 200),
                           br(),
                           sliderInput('width2', label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
                           br(),
                           actionBttn('submit4',div(
                             strong("Start analysing"),align = 'center',
                             icon("hand-right", lib = "glyphicon")))
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput05'",
                           h3(strong("Parameter ")), 
                           fileInput(inputId = 'oncogeneListFile', 
                                     label = div(style = "font-size:18px; font-weight:600; ", 'Oncogene list file'), 
                                     placeholder = "Defalut file: oncogene.list.txt", 
                                     width = 400),
                           fileInput(inputId = 'tsgListFile', 
                                     label = div(style = "font-size:18px; font-weight:600; ", 'TSG list file'), 
                                     placeholder = "Defalut file: TSG.list.txt", 
                                     width = 400),
                           selectInput("themeOption2",label = div(style = "font-size:18px; font-weight:600;  ", "Theme option"),
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
                           checkboxInput('show.percentage',label = div(style = "font-size:15px; font-weight:400; ", 'Show Percentage'),value = T),
                           sliderInput('width3',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 600,max = 1100, value = 650,width = 500),
                           br(),
                           actionBttn('submit5',div(
                             strong("Start analysing"),align = 'center',
                             icon("hand-right", lib = "glyphicon")))
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput06'",
                           h3(strong("Parameter ")), 
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
                             strong("Start analysing"),align = 'center',
                             icon("hand-right", lib = "glyphicon")))
                         )
                       )
                     ),
                     column(
                       width = 9,
                       box(
                         width = NULL,
                         tabBox(
                           id = 'tith',
                           height = "100%", 
                           width = "100%",
                           selected = "caInput02",
                           side = "left",
                           tabPanel(
                             title = div(icon("chart-bar"), "MathScore"),
                             value = "caInput02",
                             DT::dataTableOutput('mathScore'),
                             br(),
                             br(),
                             uiOutput("msdb")
                           ),
                           tabPanel(
                             title = div(icon("image"), "Vafplot"),
                             value = "caInput03",
                             conditionalPanel(
                               condition = "input.plotOption == 'separate' ",
                               uiOutput("chooselistvaf")
                             ),
                             plotOutput("vaf",height = "100%"), 
                             uiOutput("vcdb")
                           ),
                           tabPanel(
                             title = div(icon("map"), "Mutsharedprivateplot"),
                             value = "caInput04",
                             plotOutput("mutSharedPrivatePlot",height = "100%"),   
                             uiOutput("mspdb")
                           ),
                           tabPanel(
                             title = div(icon("chart-bar"), "Stackplot"),
                             value = "caInput05",
                             plotOutput("stackplot",height = "100%"),
                             br(),
                             uiOutput("stkdb")
                           ),
                           tabPanel(
                             title = div(icon("box"), "Jaccardindex"),
                             value = "caInput06",
                             plotOutput("JaccardIndex",height = "100%",width = "100%") ,
                             uiOutput("jidb")
                           )
                         )
                         
                       )
                       # tabBox(
                       #   id = 'tith',
                       #   height = "100%", 
                       #   width = "100%",
                       #   selected = "caInput02",
                       #   side = "left",
                       #   tabPanel(
                       #     title = div(icon("chart-bar"), "MathScore"),
                       #     value = "caInput02",
                       #     withSpinner(DT::dataTableOutput('mathScore')),
                       #     br(),
                       #     br(),
                       #     uiOutput("msdb")
                       #    ),
                       #   tabPanel(
                       #     title = div(icon("image"), "Vafplot"),
                       #     value = "caInput03",
                       #     conditionalPanel(
                       #       condition = "input.plotOption == 'separate' ",
                       #       uiOutput("chooselistvaf")
                       #     ),
                       #     withSpinner(plotOutput("vaf",height = "100%")), 
                       #     uiOutput("vcdb")
                       #   ),
                       #   tabPanel(
                       #     title = div(icon("map"), "Mutsharedprivateplot"),
                       #     value = "caInput04",
                       #     withSpinner(plotOutput("mutSharedPrivatePlot",height = "100%")),   
                       #     uiOutput("mspdb")
                       #   ),
                       #   tabPanel(
                       #     title = div(icon("chart-bar"), "Stackplot"),
                       #     value = "caInput05",
                       #     withSpinner(plotOutput("stackplot",height = "100%")),
                       #     br(),
                       #     uiOutput("stkdb")
                       #   ),
                       #   tabPanel(
                       #     title = div(icon("box"), "Jaccardindex"),
                       #     value = "caInput06",
                       #     withSpinner(plotOutput("JaccardIndex",height = "100%")) ,
                       #     uiOutput("jidb")
                       #   )
                       # )
                     )
                   )
)

bodyclone <- tabItem('clone',
                     h2('Clonal analysis'),
                     fluidRow(
                       column(
                         width = 3,
                         box(
                           width = NULL,
                           conditionalPanel(
                             condition = "input.clt == 'c01'",
                             h3(strong("Parameter ")), 
                             sliderInput('width5',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
                             br(),
                             br(),
                             actionBttn('submit7',div(
                               strong("Start analysing"),align = 'center',
                               icon("hand-right", lib = "glyphicon")))
                           )
                         )
                       ),
                       column(
                         width = 9,
                         box(
                           width = NULL,
                           tabBox(
                             id = 'clt',
                             selected = 'c01',
                             side = 'left',
                             height = "100%",
                             width = "100%",
                             tabPanel(
                               value = 'c01',
                               title = div(icon("newspaper"), "Tumorcloneplot"),
                               plotOutput('cloneplot',height = "100%"),
                               br(),
                               uiOutput("clpdb")
                             )
                           )
                         )
                         # tabBox(
                         #   id = 'clt',
                         #   selected = 'c01',
                         #   side = 'left',
                         #   height = "100%",
                         #   width = "100%",
                         #   tabPanel(
                         #     value = 'c01',
                         #     title = div(icon("newspaper"), "Tumorcloneplot"),
                         #     withSpinner(plotOutput('cloneplot',height = "100%")),
                         #     br(),
                         #     uiOutput("clpdb")
                         #   )
                         #  )
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
                              conditionalPanel(
                                condition = "input.fat == 'F01'",
                                h3(strong("Parameter ")),
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
                                  strong("Start analysing"),align = 'center',
                                  icon("hand-right", lib = "glyphicon")))
                              ),
                              conditionalPanel(
                                condition = "input.fat == 'F02'",
                                h3(strong("Parameter ")),
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
                                  strong("Start analysing"),align = 'center',
                                  icon("hand-right", lib = "glyphicon")))
                              )
                            )
                          ),
                          column(
                            width = 9,
                            box(
                              width = NULL,
                              height = "100%",
                              tabBox(
                                id = 'fat',
                                side = 'left',
                                selected = 'F02',
                                width = "100%",
                                height = "100%",
                                tabPanel(
                                  title = div(icon("lightbulb"), "GO analysis"),
                                  value = 'F01',
                                  uiOutput("chooselist1"),
                                  plotOutput('GOplot',height = "100%",width = "100%"),
                                  br(),
                                  uiOutput("GOdb"),
                                  uiOutput('gotui')
                                ),
                                tabPanel(
                                  title = div(icon("microsoft"), "Pathway analysis"),
                                  value = 'F02',
                                  uiOutput("chooselist2"),
                                  plotOutput('Pathwayplot',height = "100%"),
                                  br(),
                                  uiOutput("Pathdb"),
                                  uiOutput('patht')
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
                             width = 3,
                             box(
                               width = NULL,
                               conditionalPanel(
                                 condition = "input.sgt == 'S01'",
                                 h3(strong("Parameter ")),
                                 fileInput('driverGenesFile1',label = div(style = "font-size:18px; font-weight:600; ", 'Upload driverGenesFile')), 
                                 numericInput('mutThreshold1', div(style = "font-size:18px; font-weight:600;  ", 'Mutation quantity threshold'), value = 50),
                                 selectInput("signaturesRef1", label = div(style = "font-size:18px; font-weight:600;  ", "Signautre reference"),
                                             choices = c(signatures.cosmic = "signatures.cosmic",
                                                         signatures.nature2013 = "signatures.nature2013"),
                                             selected = "signatures.cosmic"),
                                 sliderInput(inputId='widthsig1',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800, width = 500),
                                 sliderInput(inputId='heightsig1',label = div(style = "font-size:18px; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 600, width = 500), 
                                 br(),
                                 br(),
                                 actionBttn('submitSig1',div(
                                   strong("Start analysing"),align = 'center',
                                   icon("hand-right", lib = "glyphicon")))
                               ),
                               conditionalPanel(
                                 condition = "input.sgt == 'S02'",
                                 h3(strong("Parameter ")),
                                 fileInput('driverGenesFile2',label = div(style = "font-size:18px; font-weight:600; ", 'Upload driverGenesFile')), 
                                 numericInput('mutThreshold2', div(style = "font-size:18px; font-weight:600;  ", 'Mutation quantity threshold'), value = 50, step=10),
                                 selectInput("signaturesRef2", label = div(style = "font-size:18px; font-weight:600;  ", "Signautre reference"),
                                             choices = c(signatures.cosmic = "signatures.cosmic",
                                                         signatures.nature2013 = "signatures.nature2013"),
                                             selected = "signatures.cosmic"),
                                 # radioButtons(
                                 #   inputId = "sigplot", 
                                 #   label = div(style = "font-size:18px; font-weight:600; ", "Plot option"), 
                                 #   choiceNames = list(
                                 #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Signature probability"), 
                                 #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Branch-trunk")
                                 #   ),
                                 #   choiceValues = c("signaturesprob", "branchtrunk"),
                                 #   selected = "signaturesprob", 
                                 #   inline = TRUE), 
                                 numericInput('signiflevel', div(style = "font-size:18px; font-weight:600;  ", 'Significant level'), value = 0.05, min=0, max=1, step=0.1),
                                 sliderInput(inputId='widthsig2',label = div(style = "font-size:18px; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800, width = 500),
                                 sliderInput(inputId='heightsig2',label = div(style = "font-size:18px; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 600, width = 500), 
                                 br(),
                                 br(),
                                 actionBttn('submitSig2',div(
                                   strong("Start analysing"),align = 'center',
                                   icon("hand-right", lib = "glyphicon")))
                               )
                             )
                           ),
                           column(
                             width = 9,
                             box(
                               width = NULL,
                               tabBox(
                                 id = 'sgt',
                                 side = 'left',
                                 selected = 'S01',
                                 width = "100%",
                                 height = "100%",
                                 tabPanel(
                                   title = div(icon("newspaper"), "Signature"), 
                                   value = 'S01',
                                   plotOutput('sigOFAPlot1', height = "100%", width = "100%"),
                                   uiOutput("sigpdb1"),
                                   br(),
                                   uiOutput('sigOFATableUI1')
                                 ), 
                                 tabPanel(
                                   title = div(icon("image"), "Branch trunck"),
                                   value = 'S02',
                                   plotOutput('sigOFAPlot2', height = "100%", width = "100%"),
                                   uiOutput("sigpdb2"),
                                   br(),
                                   uiOutput('sigOFATableUI2')
                                 )
                               )
                             )
                             # tabBox(
                             #   id = 'sgt',
                             #   side = 'left',
                             #   selected = 'S01',
                             #   width = "100%",
                             #   height = "100%",
                             #   tabPanel(
                             #     title = div(icon("newspaper"), "Summary"), 
                             #     value = 'S01',
                             #     withSpinner(DT::dataTableOutput('sigOFA'))
                             #   ), 
                             #   tabPanel(
                             #     title = div(icon("image"), "Visulization"),
                             #     value = 'S02',
                             #     withSpinner(plotOutput('sigOFA2', height = "100%", width = "100%")),
                             #     br(),
                             #     uiOutput("sigpdb")
                             #   )
                             # )
                            )
                          )
                        )

bodySurvival <- tabItem('Survival',
                        h2('Phylotree visualiaztion'),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                              width = NULL,
                              h3(strong("Parameter ")), 
                              selectInput('phyloTreeType',div(style = "font-size:18px; font-weight:600; ", 'Type'),
                                          c( 'njtree','newick','beast','PAML'),
                                          selected = 'njtree'),
                              conditionalPanel(
                                condition = "input.phyloTreeType != 'njtree'",
                                fileInput(inputId = 'phylotree.dir', 
                                          label = div(style = "font-size:18px; font-weight:600; ", 'Phylotree file'), 
                                          width = 400)
                              ),
                              conditionalPanel(
                                condition = "input.phyloTreeType == 'njtree'",
                                checkboxInput('show.mutSig',div(style = "font-size:15px; font-weight:400; ", 'Show mutation signature'),value = TRUE),
                                checkboxInput('show.heatmap',div(style = "font-size:15px; font-weight:400; ", 'Show heatmap'),value = TRUE),
                                # radioButtons(
                                #   inputId = "sig.name", 
                                #   label = div(style = "font-size:18px; font-weight:600; ", "Signature name"), 
                                #   choiceNames = list(
                                #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Default"), 
                                #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Alias")
                                #   ),
                                #   choiceValues = c("default", "alias"),
                                #   selected = "alias", 
                                #   inline = TRUE),
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
                                )
                              ),
                              actionBttn('submit10',div(
                                strong("Start analysing"),align = 'center',
                                icon("hand-right", lib = "glyphicon")))
                            )
                          ),
                          column(
                            width = 9,
                            conditionalPanel(
                              condition = 'input.submit10',
                              box(
                                width = NULL,
                                height = "100%",
                                plotOutput("phylotree",height = 800),
                                br(),
                                uiOutput("phtdb")
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
            height:5px;
        }
          
        .progress-bar {
          background-image: linear-gradient(to right, #77C7FF, #3c8dbc ) !important;
          background-size: auto !important;
          font-size:0px;
          height:5px;
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
