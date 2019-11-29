options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(deconstructSigs))
suppressMessages(library(DT))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(MesKit))

#sider bar----

sidebar <- dashboardSidebar(
  width = 300,
  sidebarMenu(id="sidername",selected='home',
              menuItem(strong("Home"), tabName = "home", icon = icon("home")),
              menuItem(strong("Input Data"), tabName = "input", icon = icon("gear")),
              menuItem(strong("ITH evaluation"), tabName = "ITH", icon = icon("bar-chart")),
              menuItem(strong("Clonal analysis"), tabName = "clone", icon = icon("bar-chart")),
              menuItem(strong("Functional exploration"), tabName = "function", icon = icon("bar-chart")),
              menuItem(strong("Mutational signature analysis"), tabName = "signature", icon = icon("bar-chart")), 
              menuItem(strong("Phylogenetic tree visualization"), tabName = "Survival", icon = icon("tree"))
  )
)

#shinydashboar
bodyHome <- tabItem("home",
                    fluidRow(
                      box(
                        width = 12,
                        status = "info",
                        solidHeader = TRUE,
                        # title = strong("Wellcome to the MesKit reporter"),
                        title = div(strong("Introduction"),style = "font-size:27px; font-weight:500;"),
                        # h3(strong("Introduction")),
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
                        title = div(strong("Overview of MesKit package"),style = "font-size:2em; font-weight:500;"),
                        fluidRow(
                          column(
                            width = 7,
                            div(img(src = "images/result_viewer.png", width = 900,height = 720),style="text-align: left;")
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
                  fluidRow(
                    column(
                      width = 12, 
                      column(
                        width = 3,
                        box(
                          title = div(shiny::icon("gear"), "Input Section", inline =TRUE),
                          width = NULL,
                          conditionalPanel(
                            condition = "true",
                            fileInput(inputId = 'maf', 
                                      label = div(style = "font-size:1.5em; font-weight:600;",'MAF file',
                                                  tags$button(
                                                    Id = "iecontrol01",
                                                    type = "button",
                                                    class = "bttn-material-circle",
                                                    class = "btn action-button",
                                                    list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                    style = " background-position: center;padding:0;margin-bottom:7px;"
                                                  )
                                      ), 
                                      placeholder = "example data: HCC6046.maf", 
                                      width = 400),
                            fileInput(inputId = 'sampleInfo', 
                                      label = div(style = "font-size:1.5em; font-weight:600; ", 'Sample information File',
                                                  tags$button(
                                                    Id = "iecontrol02",
                                                    type = "button",
                                                    class = "bttn-material-circle",
                                                    class = "btn action-button",
                                                    list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                    style = " background-position: center;padding:0;margin-bottom:7px;"
                                                  )
                                      ), 
                                      placeholder = "example data: HCC6046.sampleInfo.txt", 
                                      width = 400),
                            
                            selectInput(inputId = 'mutType', 
                                        label = div(style = "font-size:1.5em; font-weight:600; ", 'Filter option'),
                                        choices = c(All = 'All',
                                                    nonSilent = 'nonSilent'), 
                                        selected = "All", width = 400), 
                            bsTooltip(id = "mutType", 
                                      title = "Choose whether use nonsilent list to filter variant classification.", 
                                      placement = "right", 
                                      trigger = "hover"), 
                            
                            conditionalPanel(
                              condition="input.mutType == 'nonSilent'", 
                              selectInput(inputId = 'mutNonSilent', 
                                          label = div(style = "font-size:1.5em; font-weight:600; ", 'Variant classification filter(Inclusive)'),
                                          choices = c(),
                                          multiple = TRUE, 
                                          width = 400), 
                              bsTooltip(id = "mutNonSilent", 
                                        title = "Select variant classification needed to be silent", 
                                        placement = "right", 
                                        trigger = "hover") 

                            ), 
                            
                            checkboxInput('useindel', 
                                          label = div(style = "font-size:15px; ", 'use indel'), 
                                          value = FALSE, 
                                          width = 400),
                            
                            textInput(inputId = "chrSilent", 
                                      label = div(style = "font-size:1.5em; font-weight:600; ", 'Chromosome filter(Exclusive)'), 
                                      value = "NULL",
                                      placeholder = "NULL"),
                            bsTooltip(id = "chrSilent", 
                                      title = "Choose chromosomes needed to be silent in this analysis", 
                                      placement = "right", 
                                      trigger = "hover"), 
                            
                            checkboxInput(inputId = 'useccf', label = div(style = "font-size:15px; ", 'use ccf'),value = FALSE, width = 200),
                            conditionalPanel(
                              condition = "input.useccf == true",
                              fileInput('ccf.cluster',label = div(style = "font-size:1.5em; font-weight:600; ", 'ccf.cluster file',
                                                                  tags$button(
                                                                    Id = "iecontrol03",
                                                                    type = "button",
                                                                    class = "bttn-material-circle",
                                                                    class = "btn action-button",
                                                                    list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                                    style = " background-position: center;padding:0;margin-bottom:7px;"
                                                                  )
                              ),
                              placeholder = "example data: HCC6046.cluster.tsv", width = 400),
                              fileInput('ccf.loci',label = div(style = "font-size:1.5em; font-weight:600; ", 'ccf.loci file',
                                                               tags$button(
                                                                 Id = "iecontrol04",
                                                                 type = "button",
                                                                 class = "bttn-material-circle",
                                                                 class = "btn action-button",
                                                                 list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                                 style = " background-position: center;padding:0;margin-bottom:7px;"
                                                               )
                              ),
                              placeholder = "example data: HCC6046.loci.tsv", width = 400)
                            ),
                            selectInput('ref', label = div(style = "font-size:1.5em; font-weight:600; ", 'Select genome reference'),
                                        choices = c('hg19','hg38'),selected = "hg19", width = 400)
                          ),
                          actionBttn('submit1',div(
                            strong("Upload data"),align = 'center'))
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
                          box(
                            width = NULL,box(
                              width = NULL,
                              div(strong("Maf Data Preview"),style = "font-size:27px; font-weight:500;"),
                              p("MAF files contain many fields of information about chromosome and gene mutations and their annotations. The following fields are highly recommended to be contained in the MAF files.",
                                style = "font-size:16px; font-weight:500;line-height:30px;"),
                              h4(strong("Data Preview:")),
                              DT::dataTableOutput('maftable', width = '100%')
                            ),
                            uiOutput("ie1"),
                            uiOutput("ie2"),
                            uiOutput("ie3"),
                            uiOutput("ie4")
                          )
                        )
                      )
                    )
                  )



bodyITH <- tabItem("ITH",
                   # h2('ITH evaluation'),
                   # fluidRow(
                   #   box(
                   #     width = 7,
                   #     title = div(strong("ITH evaluation"),style = "font-size:27px; font-weight:500;"),
                   #     status = "info",
                   #     solidHeader = TRUE,
                   #     p("MesKit offers several functions to estimate intra-tumoral heterogeneity with mutational data of bulk sequencing, including calculating MATH score identifying shared mutations and private mutations, clustering variant allele frequencies (VAF) etc.",
                   #       style = "font-size:20px; font-weight:500;line-height:40px;")
                   #   )
                   # ),
                   fluidRow(
                     column(
                       width = 3,
                       box(
                         width = NULL,
                         conditionalPanel(
                           condition = "input.tith == 'caInput00'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(),
                           tags$table(
                             tags$tr(id = "inline", 
                                     width = "100%",
                                     tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf:")),
                                     tags$td(width = "70%", textInput(inputId = "minvaf", value = 0.02, label = NULL)))
                           ), 
                           br(),
                           tags$table(
                             tags$tr(id = "inline",
                                     width = "100%",
                                     tags$td(width = "30%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Max vaf:")),
                                     tags$td(width = "70%", textInput(inputId = "maxvaf", value = 1.00, label = NULL)))
                           ), 
                           br(),
                           br(),
                           fluidRow(
                             column(
                               width = 9,
                               div(
                                 tags$button(
                                   id = "submit0", type = "button", class = "action-button bttn",
                                   class = "bttn-unite", class = paste0("bttn-md"),
                                   class = paste0("bttn-default"),
                                   list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                                 # tags$button(
                                 #   Id = "stop2",
                                 #   type = "button",
                                 #   class = "bttn-material-circle",
                                 #   class = "btn action-button",
                                 #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                 #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                 # )
                               )
                             )
                           )
                           # actionBttn('submit2',div(
                           #   strong("Start analysing"),align = 'center',
                           #   icon("hand-right", lib = "glyphicon")))
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput02'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(),
                           tags$table(
                             tags$tr(id = "inline", 
                                     width = "100%",
                                     tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf:")),
                                     tags$td(width = "70%", textInput(inputId = "minvaf", value = 0.02, label = NULL)))
                           ), 
                           br(),
                           tags$table(
                             tags$tr(id = "inline",
                                     width = "100%",
                                     tags$td(width = "30%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Max vaf:")),
                                     tags$td(width = "70%", textInput(inputId = "maxvaf", value = 1.00, label = NULL)))
                           ), 
                           br(),
                           br(),
                           fluidRow(
                             column(
                               width = 9,
                               div(
                                 tags$button(
                                   id = "submit2", type = "button", class = "action-button bttn",
                                   class = "bttn-unite", class = paste0("bttn-md"),
                                   class = paste0("bttn-default"),
                                   list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                                 # tags$button(
                                 #   Id = "stop2",
                                 #   type = "button",
                                 #   class = "bttn-material-circle",
                                 #   class = "btn action-button",
                                 #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                 #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                 # )
                               )
                             )
                           )
                           # actionBttn('submit2',div(
                           #   strong("Start analysing"),align = 'center',
                           #   icon("hand-right", lib = "glyphicon")))
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput03'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(),
                           selectInput("plotOption", label = div(style = "font-size:1.5em; font-weight:600;  ", "Plot option"),
                                       choices = c(
                                         Compare = "compare",
                                         Combine = "combine",
                                         Separate = "separate"
                                       ), selected = "compare",width = 300),
                           selectInput("themeOption", label = div(style = "font-size:1.5em; font-weight:600;  ", "Theme option"),
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
                           checkboxInput('showMATH',label = div(style = "font-size:15px; font-weight:400; ", 'Show MATH Score'),value = TRUE),
                           sliderInput('width1', label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'), min = 700,max = 1100, value = 850,width = 500),
                           br(),
                           br(),
                           fluidRow(
                             column(
                               width = 9,
                               div(
                                 tags$button(
                                   id = "submit3", type = "button", class = "action-button bttn",
                                   class = "bttn-unite", class = paste0("bttn-md"),
                                   class = paste0("bttn-default"),
                                   list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                                 # tags$button(
                                 #   Id = "stop3",
                                 #   type = "button",
                                 #   class = "bttn-material-circle",
                                 #   class = "btn action-button",
                                 #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                 #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                 # )
                               )
                             )
                           )
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput05'",
                           div(strong("Parameter"), style = "font-size:2em; font-weight:600;"),
                           br(),
                           selectInput(inputId = 'plotChoiceSpp', 
                                       label = div(style = "font-size:1.5em; font-weight:600; ", 
                                                   'Select'), 
                                       choices = c(
                                         sharedPrivatePlot = "sharedPrivatePlot",
                                         mutOncoTSG = "mutOncoTSG"
                                       ), selected = "sharedPrivatePlot", width = 300),
                          
                           conditionalPanel(
                             condition = "input.plotChoiceSpp == 'sharedPrivatePlot'",
                             checkboxInput('show.num1',label = div(style = "font-size:15px; font-weight:400; ", 'Show mutation number'),width = 200),
                             br(),
                             sliderInput('width2', label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
                             fluidRow(
                               column(
                                 width = 9,
                                 div(
                                   tags$button(
                                     id = "submit4", type = "button", class = "action-button bttn",
                                     class = "bttn-unite", class = paste0("bttn-md"),
                                     class = paste0("bttn-default"),
                                     list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                     style = "margin-bottom:0px;margin-right:0px;"
                                   )
                                   # tags$button(
                                   #   Id = "stop4",
                                   #   type = "button",
                                   #   class = "bttn-material-circle",
                                   #   class = "btn action-button",
                                   #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                   #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                   # )
                                 )
                               )
                             )
                           ),
                           
                           conditionalPanel(
                             condition = "input.plotChoiceSpp == 'mutOncoTSG'",
                           fileInput(inputId = 'oncogeneListFile', 
                                     label = div(style = "font-size:1.5em; font-weight:600; ", 'Oncogene list file'), 
                                     placeholder = "Defalut file: oncogene.list.txt", 
                                     width = 400),
                           fileInput(inputId = 'tsgListFile', 
                                     label = div(style = "font-size:1.5em; font-weight:600; ", 'TSG list file'), 
                                     placeholder = "Defalut file: TSG.list.txt", 
                                     width = 400),
                           checkboxInput('show.percentage',label = div(style = "font-size:15px; font-weight:400; ", 'Show Percentage'),value = T),
                           sliderInput('width3',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 600,max = 1100, value = 650,width = 500),
                           br(),
                           fluidRow(
                             column(
                               width = 9,
                               div(
                                 tags$button(
                                   id = "submit5", type = "button", class = "action-button bttn",
                                   class = "bttn-unite", class = paste0("bttn-md"),
                                   class = paste0("bttn-default"),
                                   list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                                 # tags$button(
                                 #   Id = "stop5",
                                 #   type = "button",
                                 #   class = "bttn-material-circle",
                                 #   class = "btn action-button",
                                 #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                 #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                 # )
                               )
                             )
                           )
                           )
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput06'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(), 
                           selectInput("JItype",h4(strong("Type")),
                                       choices = c(
                                         Lower = "lower",
                                         Upper = "upper",
                                         Full = "full"
                                       ), selected = "lower",width = 300),
                           sliderInput('width4',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'), min = 400, max = 1000, value = 500, width = 500),
                           br(),
                           br(),
                           fluidRow(
                             column(
                               width = 9,
                               div(
                                 tags$button(
                                   id = "submit6", type = "button", class = "action-button bttn",
                                   class = "bttn-unite", class = paste0("bttn-md"),
                                   class = paste0("bttn-default"),
                                   list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                                 # tags$button(
                                 #   Id = "stop6",
                                 #   type = "button",
                                 #   class = "bttn-material-circle",
                                 #   class = "btn action-button",
                                 #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                 #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                 # )
                               )
                             )
                           )
                           # actionBttn('submit6',div(
                           #   strong("Start analysing"),align = 'center',
                           #   icon("hand-right", lib = "glyphicon")))
                         )
                       )
                     ),
                     column(
                       width = 9,
                       box(
                         width = NULL,
                         div(strong("ITH evaluation"),style = "font-size:27px; font-weight:500;"),
                         p("MesKit offers several functions to estimate intra-tumoral heterogeneity (ITH) with mutational data of bulk sequencing, including calculating MATH score, clustering variant allele frequencies (VAF) etc, identifying shared mutations and private mutations and measures similarity between samples by Jaccard similarity coefficients.",
                           style = "font-size:20px; font-weight:500;line-height:40px;"),
                         tabBox(
                           id = 'tith',
                           height = "100%", 
                           width = "100%",
                           selected = "caInput00",
                           side = "left",
                           tabPanel(
                             title = div(icon("chart-bar"), "TMB"),
                             value = "caInput00",
                             DT::dataTableOutput('mathScoreTMB'),
                             br(),
                             br(),
                             uiOutput("msdbtmb")
                           ),
                           tabPanel(
                             title = div(icon("chart-bar"), "MATH Score"),
                             value = "caInput02",
                             DT::dataTableOutput('mathScore'),
                             br(),
                             br(),
                             uiOutput("msdb")
                           ),
                           tabPanel(
                             title = div(icon("image"), "Vaf clustering"),
                             value = "caInput03",
                             conditionalPanel(
                               condition = "input.plotOption == 'separate'",
                               uiOutput("chooselistvaf")
                             ),
                             div(plotOutput("vaf",height = "100%"),align = "center"),
                             uiOutput("vcdb")
                           ),
                           tabPanel(
                             title = div(icon("map"), "TrunkOrBranch summary"),
                             value = "caInput05",
                             conditionalPanel(
                               condition = "input.plotChoiceSpp == 'sharedPrivatePlot'",
                               div(plotOutput("mutSharedPrivatePlot",height = "100%"),align ="center"),  
                               br(),
                               uiOutput("mspdb")
                             ),
                             conditionalPanel(
                               condition = "input.plotChoiceSpp == 'mutOncoTSG'",
                               div(plotOutput("mutoncotsg",height = "100%"),align ="center"),
                               br(),
                               uiOutput("stkdb")
                             )
                           ),
                           tabPanel(
                             title = div(icon("box"), "Paired-samples similarity"),
                             value = "caInput06",
                             div(plotOutput("JaccardIndex",height = "100%"),align = "center") ,
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
                     # h2('Clonal analysis'),
                     # fluidRow(
                     #   box(
                     #     width = 7,
                     #     title = div(strong("Clonal analysis"),style = "font-size:27px; font-weight:500;"),
                     #     status = "info",
                     #     solidHeader = TRUE,
                     #     p("MesKit can decipher tumor clone distribution based on CCF (Cancer Cell Frequency) data generated by PyClone. For mutational data from multi-sample across different time points, MesKit infers subclonal relationship based on R package clonevol.",
                     #       style = "font-size:20px; font-weight:500;line-height:40px;")
                     #   )
                     # ),
                     fluidRow(
                       column(
                         width = 3,
                         box(
                           width = NULL,
                           conditionalPanel(
                             condition = "input.clt == 'c01'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             br(),
                             sliderInput('width5',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
                             br(),
                             br(),
                             fluidRow(
                               column(
                                 width = 9,
                                 div(
                                   tags$button(
                                     id = "submit7", type = "button", class = "action-button bttn",
                                     class = "bttn-unite", class = paste0("bttn-md"),
                                     class = paste0("bttn-default"),
                                     list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                     style = "margin-bottom:0px;margin-right:0px;"
                                   )
                                   # tags$button(
                                   #   Id = "stop7",
                                   #   type = "button",
                                   #   class = "bttn-material-circle",
                                   #   class = "btn action-button",
                                   #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                   #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                   # )
                                 )
                               )
                             )
                           ), 
                           conditionalPanel(
                             condition = "input.clt == 'c02'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             br(),
                             checkboxInput('showdensity', label = div(style = "font-size:15px; ", 'Show density'),value = FALSE,width = 400),
                             br(),
                             sliderInput('widthccfden',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
                             br(),
                             br(),
                             fluidRow(
                               column(
                                 width = 9,
                                 div(
                                   tags$button(
                                     id = "submitccfden", type = "button", class = "action-button bttn",
                                     class = "bttn-unite", class = paste0("bttn-md"),
                                     class = paste0("bttn-default"),
                                     list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                     style = "margin-bottom:0px;margin-right:0px;"
                                   )
                                   # tags$button(
                                   #   Id = "stop7",
                                   #   type = "button",
                                   #   class = "bttn-material-circle",
                                   #   class = "btn action-button",
                                   #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                   #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                   # )
                                 )
                               )
                             )
                           )
                         )
                       ),
                       column(
                         width = 9,
                         box(
                           width = NULL,
                           div(strong("Clonal analysis"),style = "font-size:27px; font-weight:500;"),
                           p("MesKit can decipher tumor clone distribution based on CCF (Cancer Cell Frequency) data generated by PyClone. Additionally, you can compare CCF of mutations identified across samples, which can be utilized to infer whether metastases are seeded in monoclonal manner or polyclonal manner from primary tumors.",
                             style = "font-size:20px; font-weight:500;line-height:40px;"),
                           tabBox(
                             id = 'clt',
                             selected = 'c01',
                             side = 'left',
                             height = "100%",
                             width = "100%",
                             tabPanel(
                               value = 'c01',
                               title = div(icon("newspaper"), "Subclonal plot"), 
                               div(plotOutput('cloneplot', height = "100%", width = "100%"), align = "center"),
                               br(),
                               uiOutput("clpdb")
                             ),
                             tabPanel(
                               value = 'c02',
                               title = div(icon("newspaper"), "CCF Density"), 
                               div(plotOutput('ccfdenplot', height = "100%", width = "100%"), align = "center"),
                               br(),
                               uiOutput("ccfdendb")
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
                        # h2('Functional analysis'),
                        # fluidRow(
                        #   box(
                        #     width = 7,
                        #     title = div(strong("Functional analysis"),style = "font-size:27px; font-weight:500;"),
                        #     status = "info",
                        #     solidHeader = TRUE,
                        #     p("MesKit supports GO and pathway enrichment analysis (KEGG and Reactome) both on whole-tree-level and branch-level of phylogenetic tree object .",
                        #       style = "font-size:20px; font-weight:500;line-height:40px;")
                        #   )
                        # ),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                              width = NULL,
                              conditionalPanel(
                                condition = "input.fat == 'F01'",
                                div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                br(),
# <<<<<<< HEAD
#                                 selectInput("plotType1", label = div(style = "font-size:18px; font-weight:600;  ", "Plot type"),
#                                             choices = c(
#                                               Bar = "bar",
#                                               Dot = "dot"
#                                             ),
#                                             selected = "bar",width = 300),
#                                 selectInput("GO.type", label = div(style = "font-size:18px; font-weight:600;  ", "GO type"),
#                                             choices = c(
#                                               MF = "MF",
#                                               BP = "BP",
#                                               CC = "CC"
#                                             ),
#                                             selected = "BP",width = 300),
#                                 selectInput("pAdjustMethod1", label = div(style = "font-size:18px; font-weight:600;  ", "Pval adjust method"),
#                                             choices = c(Holm = "holm",
#                                                         Hochberg = "hochberg",
#                                                         Hommel = "hommel",
#                                                         Bonferroni = "bonferroni",
#                                                         BH = "BH",
#                                                         BY = "BY",
#                                                         FDR = "fdr",
#                                                         None = "none"),
#                                             selected = "BH",width = 300),
# =======
#                                 
                                selectInput(inputId = "GO.type", 
                                            label = div(style = "font-size:1.5em; font-weight:600;  ", "GO.type"),
                                            choices = c(ALL = "ALL", 
                                                        BP = "BP",
                                                        MF = "MF", 
                                                        CC = "CC"),
                                            selected = "ALL"),
                                
                                selectInput(inputId = "plotType", 
                                            label = div(style = "font-size:1.5em; font-weight:600;  ", "plotType"),
                                            choices = c(dot = "dot",
                                                        bar = "bar"),
                                            selected = "dot"), 
                                
                                selectInput(inputId = "pAdjustMethod", 
                                            label = div(style = "font-size:1.5em; font-weight:600;  ", "pAdjustMethod"),
                                            choices = c(holm = "holm", 
                                                        hochberg = "hochberg",
                                                        hommel = "hommel", 
                                                        bonferroni = "bonferroni",
                                                        BH = "BH", 
                                                        BY = "BY",
                                                        fdr = "fdr", 
                                                        none = "none"),
                                            selected = "BH"), 

                                tags$table(
                                  tags$tr(id = "inline", 
                                          width = "100%",
                                          tags$td(width = "20%", div(style = "font-size:1.5em; font-weight:600;  ", "Pval:")),
                                          tags$td(width = "70%", textInput(inputId = "pval1", value = 0.05, label = NULL)))
                                ), 
                                br(),
                                tags$table(
                                  tags$tr(id = "inline",
                                          width = "100%",
                                          tags$td(width = "20%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Qval:")),
                                          tags$td(width = "70%", textInput(inputId = "qval1", value =  0.20, label = NULL)))
                                ),
                                br(),
                                numericInput(inputId = "showCategory", 
                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "showCategory"), 
                                             value = 5),

                                sliderInput('width6',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800),
                                sliderInput('height6',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 600, value = 560),
                                br(),
                                fluidRow(
                                  column(
                                    width = 9,
                                    div(
                                      tags$button(
                                        id = "submit8", type = "button", class = "action-button bttn",
                                        class = "bttn-unite", class = paste0("bttn-md"),
                                        class = paste0("bttn-default"),
                                        list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                        style = "margin-bottom:0px;margin-right:0px;"
                                      )
                                      # tags$button(
                                      #   Id = "stop8",
                                      #   type = "button",
                                      #   class = "bttn-material-circle",
                                      #   class = "btn action-button",
                                      #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                      #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                      # )
                                    )
                                  )
                                )
                                # actionBttn('submit8',div(
                                #   strong("Start analysing"),align = 'center',
                                #   icon("hand-right", lib = "glyphicon")))
                              ),
                              conditionalPanel(
                                condition = "input.fat == 'F02'",
                                div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                br(),
# <<<<<<< HEAD
#                                 selectInput("plotType2", label = div(style = "font-size:18px; font-weight:600;  ", "Plot type"),
#                                             choices = c(
#                                               Bar = "bar",
#                                               Dot = "dot"
#                                             ),
#                                             selected = "dot",width = 300),
#                                 selectInput("pathway.type", label = div(style = "font-size:18px; font-weight:600;  ", "Pathway type"),
#                                             choices = c(
#                                               KEGG = "KEGG",
#                                               Reactome = "Reactome"
#                                             ),
#                                             selected = "BP",width = 300),
#                                 selectInput("pAdjustMethod2", label = div(style = "font-size:18px; font-weight:600;  ", "Pval adjust method"),
#                                             choices = c(Holm = "holm",
#                                                         Hochberg = "hochberg",
#                                                         Hommel = "hommel",
#                                                         Bonferroni = "bonferroni",
#                                                         BH = "BH",
#                                                         BY = "BY",
#                                                         FDR = "fdr",
#                                                         None = "none"),
#                                             selected = "BH",width = 300),
# =======
#                                 
                                
                                selectInput(inputId = "pathway.type", 
                                            label = div(style = "font-size:1.5em; font-weight:600;  ", "Pathway.type"),
                                            choices = c(KEGG = "KEGG", 
                                                        Reactome = "Reactome"),
                                            selected = "KEGG"),
                                
                                selectInput(inputId = "pathplotType", 
                                            label = div(style = "font-size:1.5em; font-weight:600;  ", "plotType"),
                                            choices = c(dot = "dot",
                                                        bar = "bar"),
                                            selected = "dot"), 
                                
                                
                                selectInput(inputId = "pathpAdjustMethod", 
                                            label = div(style = "font-size:1.5em; font-weight:600;  ", "pAdjustMethod"),
                                            choices = c(holm = "holm", 
                                                        hochberg = "hochberg",
                                                        hommel = "hommel", 
                                                        bonferroni = "bonferroni",
                                                        BH = "BH", 
                                                        BY = "BY",
                                                        fdr = "fdr", 
                                                        none = "none"),
                                            selected = "BH"), 
                                tags$table(
                                  tags$tr(id = "inline", 
                                          width = "100%",
                                          tags$td(width = "20%", div(style = "font-size:1.5em; font-weight:600; ", "Pval:")),
                                          tags$td(width = "70%", textInput(inputId = "pval2", value = 0.05, label = NULL)))
                                ), 
                                br(),
                                tags$table(
                                  tags$tr(id = "inline",
                                          width = "100%",
                                          tags$td(width = "20%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Qval:")),
                                          tags$td(width = "70%", textInput(inputId = "qval2", value =  0.20, label = NULL)))
                                ),
                                br(),
                                # tags$table(
                                #   tags$tr(id = "inline",
                                #           width = "100%",
                                #           tags$td(width = "50%", tags$div(style = "font-size:18px; font-weight:600; ", "Show category:")),
                                #           tags$td(width = "50%", textInput(inputId = "showCategory2", value =  5, label = NULL)))
                                # ),
                                # br(),
                                
                                numericInput(inputId = "pathshowCategory", 
                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "showCategory"), 
                                             value = 5),
                                
                                sliderInput('width7',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800),
                                sliderInput('height7',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 600, value = 560),
                                br(),
                                fluidRow(
                                  column(
                                    width = 9,
                                    div(
                                      tags$button(
                                        id = "submit9", type = "button", class = "action-button bttn",
                                        class = "bttn-unite", class = paste0("bttn-md"),
                                        class = paste0("bttn-default"),
                                        list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                        style = "margin-bottom:0px;margin-right:0px;"
                                      )
                                      # tags$button(
                                      #   Id = "stop9",
                                      #   type = "button",
                                      #   class = "bttn-material-circle",
                                      #   class = "btn action-button",
                                      #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                      #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                      # )
                                    )
                                  )
                                )
                                # actionBttn('submit9',div(
                                #   strong("Start analysing"),align = 'center',
                                #   icon("hand-right", lib = "glyphicon")))
                              )
                            )
                          ),
                          column(
                            width = 9,
                            box(
                              width = NULL,
                              height = "100%",
                              div(strong("Functional exploration"),style = "font-size:27px; font-weight:500;"),
                              p("MesKit supports GO and pathway enrichment analysis (KEGG/Reactome) both on tree-level and branch-level of phylogenetic tree objects.",
                                style = "font-size:20px; font-weight:500;line-height:40px;"),
                              tabBox(
                                id = 'fat',
                                side = 'left',
                                selected = 'F01',
                                width = "100%",
                                height = "100%",
                                tabPanel(
                                  title = div(icon("lightbulb"), "GO analysis"),
                                  value = 'F01',
                                  uiOutput("chooselist1"),
                                  div(plotOutput('GOplot',height = "100%",width = "100%"),align = "center"),
                                  br(),
                                  uiOutput("GOdb"),
                                  br(),
                                  br(),
                                  uiOutput('gotui')
                                ),
                                tabPanel(
                                  title = div(icon("microsoft"), "Pathway analysis"),
                                  value = 'F02',
                                  uiOutput("chooselist2"),
                                  div(plotOutput('Pathwayplot',height = "100%"),align = "center"),
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
                         # h2('Mutational signature analysis'),
                         # fluidRow(
                         #   box(
                         #     width = 7,
                         #     title = div(strong("Mutational signature analysis"),style = "font-size:27px; font-weight:500;"),
                         #     status = "info",
                         #     solidHeader = TRUE,
                         #     p("MesKit integrates mutational signature analysis by implementing R package deconstructSig, identifying potential signatures which could be attributed to known mutational processes for each branch/trunk of the NJtree object.",
                         #       style = "font-size:20px; font-weight:500;line-height:40px;")
                         #   )
                         # ),
                         fluidRow(
                           column(
                             width = 3,
                             box(
                               width = NULL,
                               conditionalPanel(
                                 condition = "input.sgt == 'S01'",
                                 div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                 br(),
                                 checkboxInput(inputId="oncogeneMapping", label = div(style = "font-size:1.5em; font-weight:600; ", 'Oncogene mapping'), value = FALSE),
                                 conditionalPanel(
                                   condition = "input.oncogeneMapping == true",
                                   fileInput(inputId = 'driverGenesFile', 
                                             label = div(style = "font-size:1.5em; font-weight:600; ", 'Oncogene list'),
                                             placeholder = "Default file: putative_driver_genes.txt", 
                                             width = 400)
                                 ), 
                                 numericInput('mutThreshold', div(style = "font-size:1.5em; font-weight:600;  ", 'Mutation quantity threshold'), value = 50),
                                 selectInput("signaturesRef", label = div(style = "font-size:1.5em; font-weight:600;  ", "Signautre reference"),
                                             choices = c(cosmic = "cosmic",
                                                         nature2013 = "nature2013"),
                                             selected = "cosmic"),
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(
                                     width = 9,
                                     div(
                                       tags$button(
                                         id = "submitSig", type = "button", class = "action-button bttn",
                                         class = "bttn-unite", class = paste0("bttn-md"),
                                         class = paste0("bttn-default"),
                                         list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                         style = "margin-bottom:0px;margin-right:0px;"
                                       )
                                       # tags$button(
                                       #   Id = "stopSig",
                                       #   type = "button",
                                       #   class = "bttn-material-circle",
                                       #   class = "btn action-button",
                                       #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                       #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                       # )
                                     )
                                   )
                                 )
                                 # actionBttn('submitSig',div(
                                 #   strong("Start analysing"),align = 'center',
                                 #   icon("hand-right", lib = "glyphicon")))
                               ),
                               conditionalPanel(
                                 condition = "input.sgt == 'S02'",
                                 div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                 br(),
                                 fileInput('driverGenesFile1',label = div(style = "font-size:1.5em; font-weight:600; ", 'Upload driverGenesFile')), 
                                 numericInput('mutThreshold1', div(style = "font-size:1.5em; font-weight:600;  ", 'Mutation quantity threshold'), value = 50),
                                 selectInput("signaturesRef1", label = div(style = "font-size:1.5em; font-weight:600;  ", "Signautre reference"),
                                             choices = c(cosmic = "cosmic",
                                                         nature2013 = "nature2013"),
                                             selected = "cosmic"),
                                 sliderInput(inputId='widthsig1',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800, width = 500),
                                 sliderInput(inputId='heightsig1',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(
                                     width = 9,
                                     div(
                                       tags$button(
                                         id = "submitSig1", type = "button", class = "action-button bttn",
                                         class = "bttn-unite", class = paste0("bttn-md"),
                                         class = paste0("bttn-default"),
                                         list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                         style = "margin-bottom:0px;margin-right:0px;"
                                       )
                                       # tags$button(
                                       #   Id = "stopSig1",
                                       #   type = "button",
                                       #   class = "bttn-material-circle",
                                       #   class = "btn action-button",
                                       #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                       #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                       # )
                                     )
                                   )
                                 )
                                 # actionBttn('submitSig1',div(
                                 #   strong("Start analysing"),align = 'center',
                                 #   icon("hand-right", lib = "glyphicon")))
                               ),
                               conditionalPanel(
                                 condition = "input.sgt == 'S03'",
                                 div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                 br(),
                                 fileInput('driverGenesFile2',label = div(style = "font-size:1.5em; font-weight:600; ", 'Upload driverGenesFile')), 
                                 numericInput('mutThreshold2', div(style = "font-size:1.5em; font-weight:600;  ", 'Mutation quantity threshold'), value = 50, step=10),
                                 selectInput("signaturesRef2", label = div(style = "font-size:1.5em; font-weight:600;  ", "Signautre reference"),
                                             choices = c(cosmic = "cosmic",
                                                         nature2013 = "nature2013"),
                                             selected = "cosmic"),
                                 # radioButtons(
                                 #   inputId = "sigplot", 
                                 #   label = div(style = "font-size:1.5em; font-weight:600; ", "Plot option"), 
                                 #   choiceNames = list(
                                 #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Signature probability"), 
                                 #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Branch-trunk")
                                 #   ),
                                 #   choiceValues = c("signaturesprob", "branchtrunk"),
                                 #   selected = "signaturesprob", 
                                 #   inline = TRUE), 
                                 numericInput('signiflevel', div(style = "font-size:1.5em; font-weight:600;  ", 'Significant level'), value = 0.05, min=0, max=1, step=0.1),
                                 sliderInput(inputId='widthsig2',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800, width = 500),
                                 sliderInput(inputId='heightsig2',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(
                                     width = 9,
                                     div(
                                       tags$button(
                                         id = "submitSig2", type = "button", class = "action-button bttn",
                                         class = "bttn-unite", class = paste0("bttn-md"),
                                         class = paste0("bttn-default"),
                                         list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                         style = "margin-bottom:0px;margin-right:0px;"
                                       )
                                       # tags$button(
                                       #   Id = "stopSig2",
                                       #   type = "button",
                                       #   class = "bttn-material-circle",
                                       #   class = "btn action-button",
                                       #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                       #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                       # )
                                     )
                                   )
                                 )
                                 # actionBttn('submitSig2',div(
                                 #   strong("Start analysing"),align = 'center',
                                 #   icon("hand-right", lib = "glyphicon")))
                               )
                             )
                           ),
                           column(
                             width = 9,
                             box(
                               width = NULL,
                               div(strong("Mutational signature analysis"),style = "font-size:27px; font-weight:500;"),
                               p("MesKit integrates mutational signature analysis by implementing R package deconstructSig, identifying potential signatures which could be attributed to known mutational processes for each branch/trunk of the NJtree object.",
                                 style = "font-size:20px; font-weight:500;line-height:40px;"),
                               tabBox(
                                 id = 'sgt',
                                 side = 'left',
                                 selected = 'S01',
                                 width = "100%",
                                 height = "100%",
                                 tabPanel(
                                   title = div(icon("newspaper"), "Summary"), 
                                   value = 'S01',
                                   DT::dataTableOutput('sigOFAt',width = "100%"),
                                   br(),
                                   uiOutput("sigpdb")
                                 ),
                                 tabPanel(
                                   title = div(icon("microsoft"), "Signature"), 
                                   value = 'S02',
                                   div(plotOutput('sigOFAPlot1', height = "100%", width = "100%"),align = "center"),
                                   uiOutput("sigpdb1"),
                                   br(),
                                   uiOutput('sigOFATableUI1')
                                 ), 
                                 tabPanel(
                                   title = div(icon("image"), "Mutational trunkOrBranch plot"),
                                   value = 'S03',
                                   div(plotOutput('sigOFAPlot2', height = "100%", width = "100%"),align = "center"),
                                   uiOutput("sigpdb2")
                                   # br(),
                                   # uiOutput('sigOFATableUI2')
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
                        # h2('Phylotree visualiaztion'),
                        # fluidRow(
                        #   box(
                        #     width = 7,
                        #     title = div(strong("Phylotree visualiaztion"),style = "font-size:27px; font-weight:500;"),
                        #     status = "info",
                        #     solidHeader = TRUE,
                        #     p("Phylogenetic tree becomes more widely used in depicting evolutionary relationships among tumors. Here, MesKit is able to plot tree-like phylogentic", 
                        #        br(),"trees with mutational signature information, along with a CCF heatmap (if CCF data is available) ",
                        #       style = "font-size:20px; font-weight:500;line-height:40px;")
                        #   )
                        #   # infoBoxOutput("progressBox2")
                        # ),
                        fluidRow(
                          column(
                            width = 3,
                            box(
                              width = NULL,
                              div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                              # selectInput('phyloTreeType',div(style = "font-size:1.5em; font-weight:600; ", 'Type'),
                              #             c( 'njtree','newick','beast','PAML'),
                              #             selected = 'njtree'),
                              # conditionalPanel(
                              #   condition = "input.phyloTreeType != 'njtree'",
                              #   fileInput(inputId = 'phylotree.dir', 
                              #             label = div(style = "font-size:1.5em; font-weight:600; ", 'Phylotree file'), 
                              #             width = 400)
                              # ),
                              checkboxInput('showmutSig',div(style = "font-size:15px; font-weight:400; ", 'Show mutation signature'),value = TRUE),
                              checkboxInput('showheatmap',div(style = "font-size:15px; font-weight:400; ", 'Show heatmap'),value = TRUE),
                                # radioButtons(
                                #   inputId = "sig.name", 
                                #   label = div(style = "font-size:1.5em; font-weight:600; ", "Signature name"), 
                                #   choiceNames = list(
                                #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Default"), 
                                #     tags$span(style = "font-size:14.5px; font-weight:400; ", "Alias")
                                #   ),
                                #   choiceValues = c("default", "alias"),
                                #   selected = "alias", 
                                #   inline = TRUE),
                              conditionalPanel(
                                condition = "input.showheatmap == true",
                                radioButtons(
                                  inputId = "heatmap.type",
                                  label = div(style = "font-size:1.5em; font-weight:600; ", "Heatmap type"),
                                  choiceNames = list(
                                    tags$span(style = "font-size:14.5px; font-weight:400; ", "Binary"), 
                                    tags$span(style = "font-size:14.5px; font-weight:400; ", "CCF")
                                  ),
                                  choiceValues = c("binary", "CCF"),
                                  selected = "binary", 
                                  inline = TRUE
                                )
                              ), 
                              fluidRow(
                                column(
                                  width = 9,
                                  div(
                                    tags$button(
                                      id = "submit10", type = "button", class = "action-button bttn",
                                      class = "bttn-unite", class = paste0("bttn-md"),
                                      class = paste0("bttn-default"),
                                      list(strong("Start analysing"),icon("hand-right", lib = "glyphicon")),
                                      style = "margin-bottom:0px;margin-right:0px;"
                                    )
                                    # tags$button(
                                    #   Id = "stop10",
                                    #   type = "button",
                                    #   class = "bttn-material-circle",
                                    #   class = "btn action-button",
                                    #   list(tags$img(src = "image/stop.png",width = "40px",height = "40px")),
                                    #   style = " background-position: center;padding:0;margin-bottom:7px;"
                                    # )
                                  )
                                )
                              )
                              # actionBttn('submit10',div(
                              #   strong("Start analysing"),align = 'center',
                              #   icon("hand-right", lib = "glyphicon")))
                            )
                          ),
                          column(
                            width = 9,
                            box(
                              width = 13,
                              div(strong("Phylogenetic tree visualization"),style = "font-size:27px; font-weight:500;"),
                              p("Phylogenetic tree becomes more widely used in depicting evolutionary relationships among tumors. Here, MesKit is able to plot tree-like phylogentic trees with mutational signature information, along with a binary mutation heatmap or CCF heatmap (if CCF data is available).", 
                                br(),"trees with mutational signature information, along with a CCF heatmap (if CCF data is available) ",
                                style = "font-size:20px; font-weight:500;line-height:40px;"),
                              conditionalPanel(
                                condition = 'input.submit10',
                                width = NULL,
                                height = "100%",
                                div(plotOutput("phylotree",height = 600,width = 950),align = "center"),
                                br(),
                                uiOutput("phtdb")
                              )
                            )
                          )
                        )
)




#Main function----
dbHeader <- dashboardHeader(title = "MesKit", titleWidth = 300, 
                            tags$li(class = "dropdown", actionLink(inputId = "help", label = div(style = "font-size:15px; font-weight:400; ", "Help"))), 
                            tags$li(class = "dropdown", actionLink(inputId = "contact", label = div(style = "font-size:15px; font-weight:400; ", "Contact"))))
dbHeader$children[[2]]$children <-  tags$a(href='https://github.com/Niinleslie/MesKit',
                                           tags$img(src='image/logo1.jpg',height='65',width='250'))
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
      # tags$head(tags$style(
      #   type="text/css", 
      #   ".progress.shiny-file-input-progress {
      #       height:5px;
      #   }
      #     
      #   .progress-bar {
      #     background-image: linear-gradient(to right, #77C7FF, #3c8dbc ) !important;
      #     background-size: auto !important;
      #     font-size:0px;
      #     height:5px;
      #   }"
      # )),
      
      tags$script(HTML('
      $(document).ready(function() {
        $("header").find("nav").append(\'<span class="textnvbar"> MesKit: Analysis and visualize multi-sample whole-exome sequencing data</span>\');
      })
     ')), 
      
      
      tags$head(
        tags$style(type="text/css", "#inline label{ display: table-cell; text-align: centers; vertical-align: middle; width=400; } 
                #inline .form-group { display: table-row; width=400; }"),
        tags$style(HTML("
        .shiny-output-error-validation {
                              color: brown;
                         }
                         .shiny-notification {
                              height: 200px;
                              width: 600px;
                              position:fixed;
                              font-size: 30px;
                              top: calc(50% - 100px);
                              left: calc(50% + 100px);
                         }
                         .dt-right {
    text-align: justify !important;
                         }


table.dataTable tbody th, table.dataTable tbody td {
    padding: 10px 1.5em !important;
}

.tooltip {
    min-width: 15em !important;
}

                         # .shiny-notification-close {
                         #      float: right;/*image size adjust  */
                         #      font-weight: bold;
                         #      font-size: 30px;
                         #      bottom: 9px;
                         #      position: relative;
                         #      padding-left: 4px;
                         #      color: #444;
                         #      cursor: default;
                         #  }
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
        tags$link(rel = "stylesheet", type = "text/css", href = "main.css")
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
