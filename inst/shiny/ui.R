options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(DT))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
# suppressMessages(library(shinycssloaders))
# suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(MesKit))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))
#suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

#sider bar----

sidebar <- dashboardSidebar(
  width = 300,
  sidebarMenu(id="sidername",selected='home',
              menuItem(strong("Home"), tabName = "home", icon = icon("home")),              
              menuItem(strong("Input Data"), tabName = "input", icon = icon("gear")),
              menuItem(strong("Alterational landscape"), tabName = "AL", icon = icon("bar-chart")),
              menuItem(strong("ITH evaluation"), tabName = "ITH", icon = icon("bar-chart")),
              menuItem(strong("Metastatic routes inference"), tabName = "clone", icon = icon("bar-chart")),
              # menuItem(strong("Functional exploration"), tabName = "function", icon = icon("bar-chart")),
              menuItem(strong("Phylogenetic tree visualization"), tabName = "tree", icon = icon("tree")) 
              # menuItem(strong("Phylogenetic tree visualization"), tabName = "Survival", icon = icon("tree"))
  )
)

#shinydashboar
bodyHome <- tabItem("home",
                    fluidRow(
                      box(
                        width = 12,
                        status = "info",
                        solidHeader = TRUE,
                        title = div(strong("Introduction"),style = "font-size:27px; font-weight:500;"),
                        p("Cancer develops as a result of the accumulation of genetic aberrations, which promotes the generation of distinct subpopulations of tumor cells and shapes intra-tumor heterogeneity (ITH). ITH is involved in tumor growth, progression, invasion, and metastasis, presenting one of the most significant barriers to accurate diagnoses and effective treatments of cancers. Therefore, dissecting and interpreting ITH of tumor dynamics is one of the major tasks in cancer research. Here, we present MesKit, an R-based package, to provide commonly used analysis and visualization modules in MRS study.",
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
                                div(img(src = "image/MesKit_workflow.png", width = "90%",height = "72%"),
                                    style="text-align: center;float:left;margin:0;padding:0")
                          ),
                          column(
                              width = 5,
                              div(
                                  h3(strong("With this MesKit Shiny APP:")),
                                  p("- Identify clonal/subclonal somatic mutations",br(),
                                    "- Quantify heterogeneity within or between tumors from the same patient",br(),
                                    "- Infer metastases clonal origin",br(),
                                    "- Perform mutational signature analysis",br(),
                                    "- Construct and visualize phylogenetic trees",
                                    style = "font-size:18px; font-weight:500;line-height:50px"),
                                  style = "text-align: left;float:left;padding-left:0px;margin:0px"
                              )
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
                          div(shiny::icon("gear"), strong("Input Section"), inline =TRUE,style = "font-size:27px; font-weight:500;"),
                          br(),
                          width = NULL,
                            fileInput(inputId = 'mafFile', 
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
                            # bsTooltip(id = "maf",
                            #           title = "Upload maf data",
                            #           placement = "right",
                            #           trigger = "hover"),
                            # checkboxInput(inputId = 'useccf', 
                            #               label = div(style = "font-size:1.5em; font-weight:600;position: relative;padding-left:15px", 'use ccf'),value = FALSE, width = 200),
                            # bsTooltip(id = "useccf",
                            #           title = "Click if provide CCF data.",
                            #           placement = "top",
                            #           trigger = "hover"),
                          fileInput('ccfFile',label = div(style = "font-size:1.5em; font-weight:600; ", 'CCF file',
                                                          tags$button(
                                                              Id = "iecontrol02",
                                                              type = "button",
                                                              class = "bttn-material-circle",
                                                              class = "btn action-button",
                                                              list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                              style = " background-position: center;padding:0;margin-bottom:7px;"
                                                          )
                          ),
                          placeholder = "example data: HCC6046.CCF.txt", width = 400),
                            # conditionalPanel(
                            #     condition = "input.useccf == true",
                            #     fileInput('ccfFile',label = div(style = "font-size:1.5em; font-weight:600; ", 'CCF file',
                            #                                     tags$button(
                            #                                         Id = "iecontrol03",
                            #                                         type = "button",
                            #                                         class = "bttn-material-circle",
                            #                                         class = "btn action-button",
                            #                                         list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                            #                                         style = " background-position: center;padding:0;margin-bottom:7px;"
                            #                                     )
                            #     ),
                            #     placeholder = "example data: HCC6046.CCF.txt", width = 400)
                            # ),
                            # selectInput(inputId = "method", label = div(style = "font-size:1.5em; font-weight:600;  ", "Tree construction approach"),
                            #             choices = c(
                            #                 "Neighbor joining" = "NJ",
                            #                 "Maximum parsimony" = "MP",
                            #                 "Maximum likelihood" = "ML",
                            #                 "FASTME.bal" = "FASTME.bal",
                            #                 "FASTME.ols" = "FASTME.ols"
                            #             ), selected = "NJ",width = 300),
                            # bsTooltip(id = "method",
                            #           title = "Approach to construct phylogenetic trees.",
                            #           placement = "top",
                            #           trigger = "hover"),
                            # selectInput(inputId = 'mutType', 
                            #             label = div(style = "font-size:1.5em; font-weight:600; ", 'Filter option'),
                            #             choices = c(All = 'All',
                            #                         nonSilent = 'nonSilent'), 
                            #             selected = "All", width = 400), 
                            # bsTooltip(id = "mutType", 
                            #           title = "Choose whether use nonsilent list to filter variant classification.", 
                            #           placement = "top", 
                            #           trigger = "hover"), 
                            # conditionalPanel(
                            #   condition="input.mutType == 'nonSilent'", 
                            #   selectInput(inputId = 'mutNonSilent', 
                            #               label = div(style = "font-size:1.5em; font-weight:600; ", 'Variant classification filter(Inclusive)'),
                            #               choices = c(),
                            #               multiple = TRUE, 
                            #               width = 400), 
                            #   bsTooltip(id = "mutNonSilent", 
                            #             title = "Select variant classification needed to be silent", 
                            #             placement = "top", 
                            #             trigger = "hover") 
                            # 
                            # ), 
                            # 
                            # checkboxInput('useindel', 
                            #               label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'use indel'), 
                            #               value = FALSE, 
                            #               width = 400),
                            # bsTooltip(id = "useindel",
                            #           title = "Whether to use INDELs besides somatic SNVs",
                            #           placement = "top",
                            #           trigger = "hover"),
                            # 
                            # textInput(inputId = "chrSilent", 
                            #           label = div(style = "font-size:1.5em; font-weight:600; ", 'Chromosome filter(Exclusive)'), 
                            #           value = "NULL",
                            #           placeholder = "NULL"),
                            # bsTooltip(id = "chrSilent",
                            #           title = "Choose chromosomes needed to be silent in this analysis",
                            #           placement = "right",
                            #           trigger = "hover"),
                            selectInput('ref', label = div(style = "font-size:1.5em; font-weight:600; ", 'Genome reference'),
                                        choices = c('hg19','hg38'),selected = "hg19", width = 400),
                        bsTooltip(id = "ref",
                                    title = "Human reference genome versions of hg19 or hg38 by UCSC",
                                    placement = "top",
                                    trigger = "hover"),
                          actionBttn('submit1',div(
                            strong("Upload data"),align = 'center'))
                        )
                      ), 
                        column(
                          width = 9, 
                          box(
                            width = NULL,
                            uiOutput("datapreview"),
                              # DT::dataTableOutput('maftable', width = '100%')
                            uiOutput("ie1"),
                            uiOutput("ie2")
                          )
                        )
                      )
                    )
                  )

bodyITH <- tabItem("ITH",
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
                           bsTooltip(id = "minvaf",
                                     title = "The minimum value of vaf",
                                     placement = "top",
                                     trigger = "hover"),
                           br(),
                           tags$table(
                             tags$tr(id = "inline",
                                     width = "100%",
                                     tags$td(width = "30%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Max vaf:")),
                                     tags$td(width = "70%", textInput(inputId = "maxvaf", value = 1.00, label = NULL)))
                           ), 
                           bsTooltip(id = "maxvaf",
                                     title = "The maximum value of vaf",
                                     placement = "top",
                                     trigger = "hover"),
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
                                   list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                               )
                             )
                           )
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
                           bsTooltip(id = "minvaf",
                                     title = "The minimum value of vaf",
                                     placement = "top",
                                     trigger = "hover"),
                           br(),
                           tags$table(
                             tags$tr(id = "inline",
                                     width = "100%",
                                     tags$td(width = "30%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Max vaf:")),
                                     tags$td(width = "70%", textInput(inputId = "maxvaf", value = 1.00, label = NULL)))
                           ), 
                           bsTooltip(id = "maxvaf",
                                     title = "The maximum value of vaf",
                                     placement = "top",
                                     trigger = "hover"),
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
                                   list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                               )
                             )
                           )
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput03'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(),
                           checkboxInput('withinTumor_vafcluster',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Use seg file'),value = FALSE),
                           bsTooltip(id = "withinTumor_vafcluster",
                                     title = "cluster vaf within tumors",
                                     placement = "top",
                                     trigger = "hover"),
                           # selectInput("plotOption", label = div(style = "font-size:1.5em; font-weight:600;  ", "Plot option"),
                           #             choices = c(
                           #               Compare = "compare",
                           #               Combine = "combine",
                           #               Separate = "separate"
                           #             ), selected = "compare",width = 300),
                           # bsTooltip(id = "plotOption",
                           #           title = "Three options for displaying clustering of VAFs in samples.",
                           #           placement = "top",
                           #           trigger = "hover"),
                           # selectInput("themeOption", label = div(style = "font-size:1.5em; font-weight:600;  ", "Theme option"),
                           #             choices = c(NPG = "npg",
                           #                         AAAS = "aaas",
                           #                         NEJM = "nejm",
                           #                         Lancet = "lancet",
                           #                         JAMA = "jama",
                           #                         JCO = "jco",
                           #                         UCSCGB = "ucscgb",
                           #                         D3 = "d3",
                           #                         LocusZoom = "locuszoom",
                           #                         IGV = "igv",
                           #                         UChicago = "uchicago",
                           #                         'Star Trek' = "startrek",
                           #                         'Tron Legacy' = 'tron',
                           #                         Futurama = "futurama",
                           #                         'Rick and Morty' = 'rickandmorty',
                           #                         'The Simpsons' = 'simpsons',
                           #                         GSEA = 'gsea'),
                           #             selected = "aaas",width = 300),
                           # bsTooltip(id = "themeOption",
                           #           title = "Select a theme palette from ggsci",
                           #           placement = "top",
                           #           trigger = "hover"),
                           # checkboxInput('showMATH',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Show MATH score'),value = TRUE),
                           # bsTooltip(id = "showMATH",
                           #           title = "Whether to show MATH Scores in the plot",
                           #           placement = "top",
                           #           trigger = "hover"),
                           # checkboxInput('useseg',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Use seg file'),value = FALSE),
                           # bsTooltip(id = "useseg",
                           #           title = "Upload copy number file to ignore variants in copy number altered regions",
                           #           placement = "top",
                           #           trigger = "hover"),
                           # conditionalPanel(
                           #     condition = "input.useseg == true",
                           #     fileInput(inputId = 'segFile', 
                           #               label = div(style = "font-size:1.5em; font-weight:600; ", 'Seg file'),
                           #               placeholder = "Default file: HCC6046.seg.txt", 
                           #               width = 400)
                           # ),
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
                                   list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                               )
                             )
                           )
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput05'",
                           div(strong("Parameter"), style = "font-size:2em; font-weight:600;"),
                           br(),
                           # selectInput(inputId = 'plotChoiceSpp', 
                           #             label = div(style = "font-size:1.5em; font-weight:600; ", 
                           #                         'Select'), 
                           #             choices = c(
                           #               sharedPrivatePlot = "sharedPrivatePlot",
                           #               mutOncoTSG = "mutOncoTSG"
                           #             ), selected = "sharedPrivatePlot", width = 300),
                           # bsTooltip(id = "plotChoiceSpp",
                           #           title = "Plot choice",
                           #           placement = "top",
                           #           trigger = "hover"),
                           checkboxInput('showNum1',label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Show mutation number'),width = 500),
                           bsTooltip(id = "showNum1",
                                     title = "Show mutational number",
                                     placement = "top",
                                     trigger = "hover"),
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
                                           list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                           style = "margin-bottom:0px;margin-right:0px;"
                                       )
                                   )
                               )
                           )
                           # conditionalPanel(
                           #   condition = "input.plotChoiceSpp == 'sharedPrivatePlot'",
                           #   checkboxInput('showNum1',label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Show mutation number'),width = 500),
                           #   bsTooltip(id = "showNum1",
                           #             title = "Show mutational number",
                           #             placement = "top",
                           #             trigger = "hover"),
                           #   br(),
                           #   sliderInput('width2', label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
                           #   fluidRow(
                           #     column(
                           #       width = 9,
                           #       div(
                           #         tags$button(
                           #           id = "submit4", type = "button", class = "action-button bttn",
                           #           class = "bttn-unite", class = paste0("bttn-md"),
                           #           class = paste0("bttn-default"),
                           #           list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                           #           style = "margin-bottom:0px;margin-right:0px;"
                           #         )
                           #       )
                           #     )
                           #   )
                           # )
                           # 
                           # conditionalPanel(
                           #   condition = "input.plotChoiceSpp == 'mutOncoTSG'",
                           # fileInput(inputId = 'oncogeneListFile', 
                           #           label = div(style = "font-size:1.5em; font-weight:600; ", 'Oncogene list file'), 
                           #           placeholder = "Defalut file: oncogene.list.txt", 
                           #           width = 400),
                           # fileInput(inputId = 'tsgListFile', 
                           #           label = div(style = "font-size:1.5em; font-weight:600; ", 'TSG list file'), 
                           #           placeholder = "Defalut file: TSG.list.txt", 
                           #           width = 400),
                           # checkboxInput('showPercentage',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Show Percentage'),value = T),
                           # bsTooltip(id = "showPercentage",
                           #           title = "Show the number of each mutations in the stack plot",
                           #           placement = "top",
                           #           trigger = "hover"),
                           # sliderInput('width3',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 600,max = 1100, value = 650,width = 500),
                           # br(),
                           # fluidRow(
                           #   column(
                           #     width = 9,
                           #     div(
                           #       tags$button(
                           #         id = "submit5", type = "button", class = "action-button bttn",
                           #         class = "bttn-unite", class = paste0("bttn-md"),
                           #         class = paste0("bttn-default"),
                           #         list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                           #         style = "margin-bottom:0px;margin-right:0px;"
                           #       )
                           #     )
                           #   )
                           # )
                           # )
                         ),
                         conditionalPanel(
                           condition = "input.tith == 'caInput06'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(), 
                           selectInput("JItype",div(style = "font-size:1.5em; font-weight:600; ", 'Type'),
                                       choices = c(
                                         Lower = "lower",
                                         Upper = "upper",
                                         Full = "full"
                                       ), selected = "lower",width = 300),
                           bsTooltip(id = "JItype",
                                     title = "Display full matrix, lower triangular or upper triangular matrix",
                                     placement = "top",
                                     trigger = "hover"),
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
                                   list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                   style = "margin-bottom:0px;margin-right:0px;"
                                 )
                               )
                             )
                           )
                         ),
                         conditionalPanel(
                             condition = "input.tith == 'caInput_ccfauc'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf:")),
                                         tags$td(width = "70%", textInput(inputId = "minccf_ccfauc", value = 0, label = NULL)))
                             ), 
                             bsTooltip(id = "minccf_ccfauc",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('withinTumor_ccfauc',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within type'),
                                           width = 500),
                             bsTooltip(id = "withinTumor_ccfauc",
                                       title = "calculate ccf within type",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             fluidRow(
                                 column(
                                     width = 9,
                                     div(
                                         tags$button(
                                             id = "submit_ccfauc", type = "button", class = "action-button bttn",
                                             class = "bttn-unite", class = paste0("bttn-md"),
                                             class = paste0("bttn-default"),
                                             list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                             style = "margin-bottom:0px;margin-right:0px;"
                                         )
                                     )
                                 )
                             )
                         ),
                         conditionalPanel(
                             condition = "input.tith == 'caInput_calfst'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf:")),
                                         tags$td(width = "70%", textInput(inputId = "minvaf_calfst", value = 0.08, label = NULL)))
                             ), 
                             bsTooltip(id = "minvaf_calfst",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('withinTumor_calfst',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within type'),
                                           width = 500),
                             bsTooltip(id = "withinTumor_calfst",
                                       title = "calculate ccf within type",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             fluidRow(
                                 column(
                                     width = 9,
                                     div(
                                         tags$button(
                                             id = "submit_calfst", type = "button", class = "action-button bttn",
                                             class = "bttn-unite", class = paste0("bttn-md"),
                                             class = paste0("bttn-default"),
                                             list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                             style = "margin-bottom:0px;margin-right:0px;"
                                         )
                                     )
                                 )
                             )
                         ),
                         conditionalPanel(
                             condition = "input.tith == 'caInput_calneidist'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf:")),
                                         tags$td(width = "70%", textInput(inputId = "minvaf_calneidist", value = 0.08, label = NULL)))
                             ), 
                             bsTooltip(id = "minvaf_calneidist",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('withinTumor_calneidist',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within type'),
                                           width = 500),
                             bsTooltip(id = "withinTumor_calneidist",
                                       title = "calculate ccf within type",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             fluidRow(
                                 column(
                                     width = 9,
                                     div(
                                         tags$button(
                                             id = "submit_calneidist", type = "button", class = "action-button bttn",
                                             class = "bttn-unite", class = paste0("bttn-md"),
                                             class = paste0("bttn-default"),
                                             list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                             style = "margin-bottom:0px;margin-right:0px;"
                                         )
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
                         div(strong("ITH evaluation"),style = "font-size:27px; font-weight:500;"),
                         p("MesKit offers several functions to estimate intra-tumoral heterogeneity (ITH) with mutational data of bulk sequencing, including calculating MATH score, clustering variant allele frequencies (VAF), identifying shared mutations/private mutations and measuring similarity between samples by Jaccard similarity coefficients.",
                           style = "font-size:20px; font-weight:500;line-height:40px;"),
                         tabBox(
                           id = 'tith',
                           height = "100%", 
                           width = "100%",
                           selected = "caInput02",
                           side = "left",
                           # tabPanel(
                           #   title = div(icon("chart-bar"), "TMB"),
                           #   value = "caInput00",
                           #   uiOutput('warningMessage01'),
                           #   DT::dataTableOutput('mathScoreTMB'),
                           #   br(),
                           #   br(),
                           #   uiOutput("msdbtmb")
                           # ),
                           tabPanel(
                             title = div(icon("chart-bar"), "MATH score"),
                             value = "caInput02",
                             uiOutput('warningMessage02'),
                             DT::dataTableOutput('mathScore'),
                             br(),
                             br(),
                             uiOutput("msdb")
                           ),
                           tabPanel(
                             title = div(icon("image"), "VAF clustering"),
                             value = "caInput03",
                             uiOutput("vafcluster.patientlist"),
                             uiOutput("vafcluster.samplelist"),
                             uiOutput('warningMessage03'),
                             div(plotOutput("vaf",height = "100%"),align = "center"),
                             uiOutput("vcdb")
                           ),
                           # tabPanel(
                           #   title = div(icon("map"), "TrunkOrBranch summary"),
                           #   value = "caInput05",
                           #   conditionalPanel(
                           #     condition = "input.plotChoiceSpp == 'sharedPrivatePlot'",
                           #     uiOutput('warningMessage04'),
                           #     div(plotOutput("mutSharedPrivatePlot",height = "100%"),align ="center"),  
                           #     br(),
                           #     uiOutput("mspdb")
                           #   ),
                           #   conditionalPanel(
                           #     condition = "input.plotChoiceSpp == 'mutOncoTSG'",
                           #     uiOutput('warningMessage05'),
                           #     div(plotOutput("mutoncotsg",height = "100%"),align ="center"),
                           #     br(),
                           #     uiOutput("stkdb")
                           #   )
                           # ),
                           # tabPanel(
                           #   title = div(icon("box"), "Paired-samples similarity"),
                           #   value = "caInput06",
                           #   uiOutput('warningMessage06'),
                           #   div(plotOutput("JaccardIndex",height = "100%"),align = "center") ,
                           #   uiOutput("jidb")
                           # ),
                           tabPanel(
                               title = div(icon("box"), "ccfAUC"),
                               value = "caInput_ccfauc",
                               uiOutput('warningMessage_ccfauc'),
                               uiOutput('ccfauc.patientlist'),
                               div(plotOutput("ccfauc_plot",height = "100%"),align = "center") ,
                               uiOutput("ccfauc_db_ui"),
                               uiOutput("ccfauc_table_ui")
                           ),
                           tabPanel(
                               title = div(icon("box"), "calFst"),
                               value = "caInput_calfst",
                               uiOutput('warningMessage_calfst'),
                               uiOutput('calfst.patientlist'),
                               div(plotOutput("calfst_plot",height = "100%"),align = "center") ,
                               uiOutput("calfst_db_ui"),
                               uiOutput("calfst_avg_table_ui"),
                               uiOutput("calfst_pair_table_ui")
                           ),
                           tabPanel(
                               title = div(icon("box"), "calNeiDist"),
                               value = "caInput_calneidist",
                               uiOutput('calneidist.patientlist'),
                               uiOutput('warningMessage_calneidist'),
                               div(plotOutput("calneidist_plot",height = "100%"),align = "center") ,
                               uiOutput("calneidist_db_ui"),
                               uiOutput("calneidist_avg_table_ui"),
                               uiOutput("calneidist_pair_table_ui")
                           )
                         )
                       )
                     )
                   )
)

bodyAL <- tabItem("AL",
                  fluidRow(
                      column(
                          width = 3,
                          box(
                              width = NULL,
                              conditionalPanel(
                                  condition = "input.al_tabbox == 'pannel_classifymut'",
                                  div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                  br(),
                                  fluidRow(
                                      column(
                                          width = 9,
                                          div(
                                              tags$button(
                                                  id = "submit_classifymut", type = "button", class = "action-button bttn",
                                                  class = "bttn-unite", class = paste0("bttn-md"),
                                                  class = paste0("bttn-default"),
                                                  list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                                  style = "margin-bottom:0px;margin-right:0px;"
                                              )
                                          )
                                      )
                                  )
                              ),
                              conditionalPanel(
                                  condition = "input.al_tabbox == 'pannel_plotcna'",
                                  div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                  br(),
                                  fileInput(inputId = 'segFile_plotcna', 
                                            label = div(style = "font-size:1.5em; font-weight:600; ", 'Seg file'),
                                            placeholder = "Default file: HCC6046.seg.tsv", 
                                            width = 400),
                                  fluidRow(
                                      column(
                                          width = 9,
                                          div(
                                              tags$button(
                                                  id = "submit_plotcna", type = "button", class = "action-button bttn",
                                                  class = "bttn-unite", class = paste0("bttn-md"),
                                                  class = paste0("bttn-default"),
                                                  list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                                  style = "margin-bottom:0px;margin-right:0px;"
                                              )
                                          )
                                      )
                                  )
                              ),
                          )
                      ),
                      column(
                          width = 9,
                          box(
                              width = NULL,
                              div(strong("Alterational Landscape"),style = "font-size:27px; font-weight:500;"),
                              p("",
                                style = "font-size:20px; font-weight:500;line-height:40px;"),
                              tabBox(
                                  id = 'al_tabbox',
                                  selected = 'pannel_classifymut',
                                  side = 'left',
                                  height = "100%",
                                  width = "100%",
                                  tabPanel(
                                      value = 'pannel_classifymut',
                                      title = div(icon("newspaper"), "Mutational landscape"), 
                                      uiOutput('warningMessage_classifymut'),
                                      div(plotOutput('classifymut_plot', height = "100%", width = "100%"), align = "center"),
                                      br(),
                                      uiOutput("classifymut_db_ui"),
                                      uiOutput("classifymut_table_ui")
                                  ),
                                  tabPanel(
                                      value = 'pannel_plotcna',
                                      title = div(icon("newspaper"), "CNA profile"), 
                                      uiOutput('warningMessage_plotcna'),
                                      div(plotOutput('plotcna_plot', height = "100%", width = "100%"), align = "center"),
                                      br(),
                                      uiOutput("plotcna_db_ui"),
                                      uiOutput("plotcna_table_ui")
                                  )
                              )
                          )
                      )
                )
)

bodyclone <- tabItem('clone',
                     fluidRow(
                       column(
                         width = 3,
                         box(
                           width = NULL,
                           conditionalPanel(
                             condition = "input.clt == 'clone_compareccf'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             br(),
                             # checkboxInput('showdensity', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Show density'),value = FALSE,width = 400),
                             # bsTooltip(id = "showdensity",
                             #           title = "If TRUE, the default, perform density estimation",
                             #           placement = "top",
                             #           trigger = "hover"),
                             # sliderInput('widthccfden',label = div(style = "font-size:1.5em; font-weight:600;", 'Image width'),min = 700,max = 1100, value = 850,width = 500),
                             # br(),
                             # br(),
                             fluidRow(
                               column(
                                 width = 9,
                                 div(
                                   tags$button(
                                     id = "submit_compareccf", type = "button", class = "action-button bttn",
                                     class = "bttn-unite", class = paste0("bttn-md"),
                                     class = paste0("bttn-default"),
                                     list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                     style = "margin-bottom:0px;margin-right:0px;"
                                   )
                                 )
                               )
                             )
                           ),
                           conditionalPanel(
                               condition = "input.clt == 'clone_comparejsi'",
                               div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf:")),
                                           tags$td(width = "70%", textInput(inputId = "minvaf_comparejsi", value = 0.08, label = NULL)))
                               ), 
                               bsTooltip(id = "minvaf_comparejsi",
                                         title = "The minimum value of vaf",
                                         placement = "top",
                                         trigger = "hover"),
                               checkboxInput('pairByTumor_comparejsi',
                                             value = FALSE,
                                             label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Pair by type'),
                                             width = 500),
                               bsTooltip(id = "pairByTumor_comparejsi",
                                         title = "calculate JSI within type",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               fluidRow(
                                   column(
                                       width = 9,
                                       div(
                                           tags$button(
                                               id = "submit_comparejsi", type = "button", class = "action-button bttn",
                                               class = "bttn-unite", class = paste0("bttn-md"),
                                               class = paste0("bttn-default"),
                                               list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                               style = "margin-bottom:0px;margin-right:0px;"
                                           )
                                       )
                                   )
                               )
                           ),
                           conditionalPanel(
                               condition = "input.clt == 'clone_testneutral'",
                               div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                           tags$td(width = "70%", textInput(inputId = "minvaf_testneutral", value = 0.1, label = NULL)))
                               ), 
                               bsTooltip(id = "minvaf_testneutral",
                                         title = "The minimum value of vaf",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Max vaf: ")),
                                           tags$td(width = "70%", textInput(inputId = "maxvaf_testneutral", value = 0.3, label = NULL)))
                               ), 
                               bsTooltip(id = "maxvaf_testneutral",
                                         title = "The maximum value of vaf",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               fluidRow(
                                   column(
                                       width = 9,
                                       div(
                                           tags$button(
                                               id = "submit_testneutral", type = "button", class = "action-button bttn",
                                               class = "bttn-unite", class = paste0("bttn-md"),
                                               class = paste0("bttn-default"),
                                               list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                               style = "margin-bottom:0px;margin-right:0px;"
                                           )
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
                             selected = 'clone_compareccf',
                             side = 'left',
                             height = "100%",
                             width = "100%",
                             tabPanel(
                               value = 'clone_compareccf',
                               title = div(icon("newspaper"), "compareCCF"), 
                               uiOutput('warningMessage_compareccf'),
                               uiOutput('compareccf.patientlist'),
                               uiOutput('compareccf.samplelist'), 
                               br(),
                               uiOutput("compareccf_table_ui")
                             ),
                             tabPanel(
                                 title = div(icon("box"), "compareJSI"),
                                 value = "clone_comparejsi",
                                 uiOutput('comparejsi.patientlist'),
                                 uiOutput('warningMessage_comparejsi'),
                                 div(plotOutput("comparejsi_plot",height = "100%"),align = "center") ,
                                 uiOutput("comparejsi_db_ui"),
                                 uiOutput("comparejsi_avg_table_ui"),
                                 uiOutput("comparejsi_pair_table_ui")
                             ),
                             tabPanel(
                                 title = div(icon("box"), "testNeutral"),
                                 value = "clone_testneutral",
                                 uiOutput('testneutral.patientlist'),
                                 uiOutput('testneutral.samplelist'),
                                 uiOutput('warningMessage_testneutral'),
                                 div(plotOutput("testneutral_plot",height = "100%"),align = "center") ,
                                 uiOutput("testneutral_db_ui"),
                                 uiOutput("testneutral_table_ui")
                             )
                           )
                         )
                       )
                     )
)    


# bodyfunction <- tabItem('function',
#                         fluidRow(
#                           column(
#                             width = 3,
#                             box(
#                               width = NULL,
#                               conditionalPanel(
#                                 condition = "input.fat == 'F01'",
#                                 div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
#                                 br(),
#                                 checkboxInput(inputId="driverGenesMapping1", label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Driver genes mapping'), value = FALSE),
#                                 bsTooltip(id = "driverGenesMapping1",
#                                           title = 'The file with driver gene list.',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 conditionalPanel(
#                                     condition = "input.driverGenesMapping1 == true",
#                                     fileInput(inputId = 'driverGenesFile1', 
#                                               label = div(style = "font-size:1.5em; font-weight:600; ", 'Driver genes list'),
#                                               placeholder = "Default file: putative_driver_genes.txt", 
#                                               width = 400)
#                                 ),
#                                 selectInput(inputId = "GOtype", 
#                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "GO type"),
#                                             choices = c(ALL = "ALL", 
#                                                         BP = "BP",
#                                                         MF = "MF", 
#                                                         CC = "CC"),
#                                             selected = "BP"),
#                                 bsTooltip(id = "GOtype",
#                                           title = "One of BP, MF, and CC sub-ontologies, or ALL for pooling 3 GO sub-ontologies.",
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 
#                                 selectInput(inputId = "plotType", 
#                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "Plot type"),
#                                             choices = c(dot = "dot",
#                                                         bar = "bar"),
#                                             selected = "dot"),
#                                 bsTooltip(id = "plotType",
#                                           title = "One of dot, bar",
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 
#                                 selectInput(inputId = "pAdjustMethod", 
#                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "pAdjustMethod"),
#                                             choices = c(holm = "holm", 
#                                                         hochberg = "hochberg",
#                                                         hommel = "hommel", 
#                                                         bonferroni = "bonferroni",
#                                                         BH = "BH", 
#                                                         BY = "BY",
#                                                         fdr = "fdr", 
#                                                         none = "none"),
#                                             selected = "BH"), 
#                                 bsTooltip(id = "pAdjustMethod",
#                                           title = "Method to adjust P value.",
#                                           placement = "top",
#                                           trigger = "hover"),
# 
#                                 tags$table(
#                                   tags$tr(id = "inline", 
#                                           width = "100%",
#                                           tags$td(width = "20%", div(style = "font-size:1.5em; font-weight:600;  ", "P-value:")),
#                                           tags$td(width = "70%", textInput(inputId = "pval1", value = 0.05, label = NULL)))
#                                 ), 
#                                 bsTooltip(id = "pval1",
#                                           title = "Cutoff value of pvalue. Default pval=0.05",
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 br(),
#                                 tags$table(
#                                   tags$tr(id = "inline",
#                                           width = "100%",
#                                           tags$td(width = "20%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Q-value:")),
#                                           tags$td(width = "70%", textInput(inputId = "qval1", value =  0.20, label = NULL)))
#                                 ),
#                                 bsTooltip(id = "qval1",
#                                           title = "Cutoff value of qvalue. Default qval=0.20",
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 br(),
#                                 numericInput(inputId = "showCategory", 
#                                              label = div(style = "font-size:1.5em; font-weight:600;  ", "Show category"), 
#                                              value = 5),
#                                 bsTooltip(id = "showCategory",
#                                           title = "Category numbers",
#                                           placement = "top",
#                                           trigger = "hover"),
# 
#                                 sliderInput('width6',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800),
#                                 sliderInput('height6',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 600, value = 560),
#                                 br(),
#                                 fluidRow(
#                                   column(
#                                     width = 9,
#                                     div(
#                                       tags$button(
#                                         id = "submit8", type = "button", class = "action-button bttn",
#                                         class = "bttn-unite", class = paste0("bttn-md"),
#                                         class = paste0("bttn-default"),
#                                         list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
#                                         style = "margin-bottom:0px;margin-right:0px;"
#                                       )
#                                     )
#                                   )
#                                 )
#                               ),
#                               conditionalPanel(
#                                 condition = "input.fat == 'F02'",
#                                 div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
#                                 br(),
#                                 checkboxInput(inputId="driverGenesMapping2", label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Driver genes mapping'), value = FALSE),
#                                 bsTooltip(id = "driverGenesMapping2",
#                                           title = 'The file with driver gene list.',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 conditionalPanel(
#                                     condition = "input.driverGenesMapping2 == true",
#                                     fileInput(inputId = 'driverGenesFile2', 
#                                               label = div(style = "font-size:1.5em; font-weight:600; ", 'Driver genes list'),
#                                               placeholder = "Default file: putative_driver_genes.txt", 
#                                               width = 400)
#                                 ),
#                                 selectInput(inputId = "pathwaytype", 
#                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "Pathway type"),
#                                             choices = c(KEGG = "KEGG", 
#                                                         Reactome = "Reactome"),
#                                             selected = "KEGG"),
#                                 bsTooltip(id = "pathwaytype",
#                                           title = 'One of "KEGG" or "Reactome". Default type="KEGG"',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 
#                                 selectInput(inputId = "pathplotType", 
#                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "Plot type"),
#                                             choices = c(dot = "dot",
#                                                         bar = "bar"),
#                                             selected = "dot"), 
#                                 bsTooltip(id = "pathplotType",
#                                           title = 'One of "dot", "bar"',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 selectInput(inputId = "pathpAdjustMethod", 
#                                             label = div(style = "font-size:1.5em; font-weight:600;  ", "pAdjustMethod"),
#                                             choices = c(holm = "holm", 
#                                                         hochberg = "hochberg",
#                                                         hommel = "hommel", 
#                                                         bonferroni = "bonferroni",
#                                                         BH = "BH", 
#                                                         BY = "BY",
#                                                         fdr = "fdr", 
#                                                         none = "none"),
#                                             selected = "BH"), 
#                                 bsTooltip(id = "pathpAdjustMethod",
#                                           title = 'One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none".',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 tags$table(
#                                   tags$tr(id = "inline", 
#                                           width = "100%",
#                                           tags$td(width = "20%", div(style = "font-size:1.5em; font-weight:600; ", "P-value:")),
#                                           tags$td(width = "70%", textInput(inputId = "pval2", value = 0.05, label = NULL)))
#                                 ), 
#                                 bsTooltip(id = "pval2",
#                                           title = 'Cutoff value of pvalue. Default pval=0.05',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 br(),
#                                 tags$table(
#                                   tags$tr(id = "inline",
#                                           width = "100%",
#                                           tags$td(width = "20%", tags$div(style = "font-size:1.5em; font-weight:600; ", "Q-value:")),
#                                           tags$td(width = "70%", textInput(inputId = "qval2", value =  0.20, label = NULL)))
#                                 ),
#                                 bsTooltip(id = "qval2",
#                                           title = 'Cutoff value of qvalue. Default qval=0.20',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 br(),
#                                 numericInput(inputId = "pathshowCategory", 
#                                              label = div(style = "font-size:1.5em; font-weight:600;  ", "Show category"), 
#                                              value = 5),
#                                 bsTooltip(id = "pathshowCategory",
#                                           title = 'Category numbers',
#                                           placement = "top",
#                                           trigger = "hover"),
#                                 
#                                 sliderInput('width7',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800),
#                                 sliderInput('height7',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 600, value = 560),
#                                 br(),
#                                 fluidRow(
#                                   column(
#                                     width = 9,
#                                     div(
#                                       tags$button(
#                                         id = "submit9", type = "button", class = "action-button bttn",
#                                         class = "bttn-unite", class = paste0("bttn-md"),
#                                         class = paste0("bttn-default"),
#                                         list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
#                                         style = "margin-bottom:0px;margin-right:0px;"
#                                       )
#                                     )
#                                   )
#                                 )
#                               )
#                             )
#                           ),
#                           column(
#                             width = 9,
#                             box(
#                               width = NULL,
#                               height = "100%",
#                               div(strong("Functional exploration"),style = "font-size:27px; font-weight:500;"),
#                               p("MesKit supports GO and pathway enrichment analysis (KEGG/Reactome) both on tree-level and branch-level of phylogenetic tree objects.",
#                                 style = "font-size:20px; font-weight:500;line-height:40px;"),
#                               tabBox(
#                                 id = 'fat',
#                                 side = 'left',
#                                 selected = 'F01',
#                                 width = "100%",
#                                 height = "100%",
#                                 tabPanel(
#                                   title = div(icon("lightbulb"), "GO analysis"),
#                                   value = 'F01',
#                                   uiOutput("chooselist1"),
#                                   uiOutput('warningMessage09'),
#                                   div(plotOutput('GOplot',height = "100%",width = "100%"),align = "center"),
#                                   br(),
#                                   uiOutput("GOdb"),
#                                   br(),
#                                   br(),
#                                   uiOutput('gotui')
#                                 ),
#                                 tabPanel(
#                                   title = div(icon("microsoft"), "Pathway analysis"),
#                                   value = 'F02',
#                                   uiOutput("chooselist2"),
#                                   uiOutput('warningMessage10'),
#                                   div(plotOutput('Pathwayplot',height = "100%"),align = "center"),
#                                   br(),
#                                   uiOutput("Pathdb"),
#                                   uiOutput('patht')
#                                 )
#                               )
#                             )
#                           )
#                         )
# )



bodytree <- tabItem('tree',
                         fluidRow(
                           column(
                             width = 3,
                             box(
                               width = NULL,
                               conditionalPanel(
                                   condition = "input.sgt == 'S_phylotree'",
                                   div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                   checkboxInput('showmutSig',div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Show mutation signature'),value = TRUE),
                                   bsTooltip(id = "showmutSig",
                                             title = 'Whether to show mutation signatures on branch or trunk',
                                             placement = "top",
                                             trigger = "hover"),
                                   checkboxInput('showheatmap',div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Show heatmap'),value = TRUE),
                                   bsTooltip(id = "showheatmap",
                                             title = 'Whether to show heatmap of somatic mutations. Default is TRUE.',
                                             placement = "top",
                                             trigger = "hover"),
                                   checkboxInput('showbootstrap',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Show bootstrap value'),value = TRUE),
                                   bsTooltip(id = "showbootstrap",
                                             title = 'Whether to add bootstrap value on internal nodes.',
                                             placement = "top",
                                             trigger = "hover"),
                                   conditionalPanel(
                                       condition = "input.showbootstrap == true",
                                       checkboxInput('usebox',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Use box'),value = TRUE),
                                       bsTooltip(id = "usebox",
                                                 title = 'Whether to add box arround bootstrap value on tree. ',
                                                 placement = "top",
                                                 trigger = "hover")
                                   ), 
                                   fluidRow(
                                       column(
                                           width = 9,
                                           div(
                                               tags$button(
                                                   id = "submit_phylotree", type = "button", class = "action-button bttn",
                                                   class = "bttn-unite", class = paste0("bttn-md"),
                                                   class = paste0("bttn-default"),
                                                   list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                                   style = "margin-bottom:0px;margin-right:0px;"
                                               )
                                           )
                                       )
                                   )
                               ), 
                               conditionalPanel(
                                 condition = "input.sgt == 'S_treemutsig'",
                                 div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                 br(),
                                 fileInput('driverGenesFile1',label = div(style = "font-size:1.5em; font-weight:600; ", 'Upload driver genes file')), 
                                 numericInput('mutThreshold1', div(style = "font-size:1.5em; font-weight:600;  ", 'Mutation quantity threshold'), value = 15),
                                 bsTooltip(id = "mutThreshold1",
                                           title = 'The threshold for the variants in a branch. Default 15.". Option: "nature2013". ',
                                           placement = "top",
                                           trigger = "hover"),
                                 selectInput("signaturesRef1", label = div(style = "font-size:1.5em; font-weight:600;  ", "Signautre reference"),
                                             choices = c("cosmic_v2",
                                                         "nature2013",
                                                         "geome_cosmic_v3",
                                                         "exome_cosmic_v3"),
                                             selected = "cosmic"),
                                 bsTooltip(id = "signaturesRef1",
                                           title = 'The parameter used for deconstructSig. Default "cosmic". Option: "nature2013". ". Option: "nature2013". ',
                                           placement = "top",
                                           trigger = "hover"),
                                 sliderInput(inputId='width_treemutsig',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800, width = 500),
                                 sliderInput(inputId='height_treemutsig',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(
                                     width = 9,
                                     div(
                                       tags$button(
                                         id = "submit_treemutsig", type = "button", class = "action-button bttn",
                                         class = "bttn-unite", class = paste0("bttn-md"),
                                         class = paste0("bttn-default"),
                                         list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                         style = "margin-bottom:0px;margin-right:0px;"
                                       )
                                     )
                                   )
                                 )
                               ), 
                               conditionalPanel(
                                 condition = "input.sgt == 'S_muttrunkbranch'",
                                 div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                 br(),
                                 fileInput('driverGenesFile2',label = div(style = "font-size:1.5em; font-weight:600; ", 'Upload driverGenesFile')), 
                                 numericInput('mutThreshold2', div(style = "font-size:1.5em; font-weight:600;  ", 'Mutation quantity threshold'), value = 15, step=10),
                                 bsTooltip(id = "mutThreshold2",
                                           title = 'the threshold for the variants in a branch. Default 15.',
                                           placement = "top",
                                           trigger = "hover"),
                                 selectInput("signaturesRef2", label = div(style = "font-size:1.5em; font-weight:600;  ", "Signautre reference"),
                                             choices = c("cosmic_v2",
                                                         "nature2013",
                                                         "geome_cosmic_v3",
                                                         "exome_cosmic_v3"),
                                             selected = "cosmic"),
                                 bsTooltip(id = "signaturesRef2",
                                           title = 'The parameter used for deconstructSig. Default "cosmic". Option: "nature2013".',
                                           placement = "top",
                                           trigger = "hover"),
                                 numericInput('conflevel', div(style = "font-size:1.5em; font-weight:600;  ", 'Confidence level'), value = 0.95, min=0, max=1, step=0.05),
                                 bsTooltip(id = "conflevel",
                                           title = 'Confidence level of the interval for wilcox.test. Default: 0.95. Option: on the scale of 0 to 1.',
                                           placement = "top",
                                           trigger = "hover"),
                                 sliderInput(inputId='width_muttrunkbranch',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 470, width = 500),
                                 sliderInput(inputId='height_muttrunkbranch',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 470, width = 500), 
                                 br(),
                                 br(),
                                 fluidRow(
                                   column(
                                     width = 9,
                                     div(
                                       tags$button(
                                         id = "submit_muttrunkbranch", type = "button", class = "action-button bttn",
                                         class = "bttn-unite", class = paste0("bttn-md"),
                                         class = paste0("bttn-default"),
                                         list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
                                         style = "margin-bottom:0px;margin-right:0px;"
                                       )
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
                               div(strong("Mutational signature analysis"),style = "font-size:27px; font-weight:500;"),
                               p("MesKit integrates mutational signature analysis by implementing R package deconstructSig, identifying potential signatures which could be attributed to known mutational processes for each branch/trunk of the NJtree object.",
                                 style = "font-size:20px; font-weight:500;line-height:40px;"),
                               tabBox(
                                 id = 'sgt',
                                 side = 'left',
                                 selected = 'S_phylotree',
                                 width = "100%",
                                 height = "100%",
                                 tabPanel(
                                     title = div(icon("newspaper"), "Plot phylotree"),
                                     value = 'S_phylotree',
                                     uiOutput("phylotree.patientlist"),
                                     uiOutput('warningMessage15'),
                                     div(plotOutput("phylotree_plot",height = 600,width = 900),align = "center"),
                                     br(),
                                     uiOutput("phylotree_downloadbutton_ui")
                                 ),
                                 tabPanel(
                                     title = div(icon("image"), "Mutational trunkOrBranch plot"),
                                     value = 'S_muttrunkbranch',
                                     # uiOutput('warningMessage_muttrunkbranch'),
                                     uiOutput("muttrunkbranch.patientlist"),
                                     div(plotOutput('muttrunkbranch_plot', height = "100%", width = "100%"),align = "center"),
                                     uiOutput("muttrunkbranch_download_button_ui"),
                                     br(),
                                     uiOutput('muttrunkbranch_table_ui')
                                 ), 
                                 tabPanel(
                                   title = div(icon("microsoft"), "Signature plot"), 
                                   value = 'S_treemutsig',
                                   # uiOutput('warningMessage_treemutsig'),
                                   uiOutput("treemutsig.patientlist"),
                                   div(plotOutput('treemutsig_plot', height = "100%", width = "100%"),align = "center"),
                                   uiOutput("treemutsig_download_button_ui"),
                                   br(),
                                   uiOutput('treemutsig_table_ui')
                                 )
                               )
                             )
                           )
                         )
)

# bodySurvival <- tabItem('Survival',
#                         fluidRow(
#                           column(
#                             width = 3,
#                             box(
#                               width = NULL,
#                               div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
#                               checkboxInput('showmutSig',div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Show mutation signature'),value = TRUE),
#                               bsTooltip(id = "showmutSig",
#                                         title = 'Whether to show mutation signatures on branch or trunk',
#                                         placement = "top",
#                                         trigger = "hover"),
#                               checkboxInput('showheatmap',div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Show heatmap'),value = TRUE),
#                               bsTooltip(id = "showheatmap",
#                                         title = 'Whether to show heatmap of somatic mutations. Default is TRUE.',
#                                         placement = "top",
#                                         trigger = "hover"),
#                               # conditionalPanel(
#                               #     condition = "input.showheatmap == true",
#                               #     radioButtons(
#                               #         inputId = "heatmaptype",
#                               #         label = div(style = "font-size:1.5em; font-weight:600;  ", "Heatmap type"),
#                               #         choiceNames = list(
#                               #             tags$span(style = "font-size:1.5em; font-weight:600; ", "Binary"), 
#                               #             tags$span(style = "font-size:1.5em; font-weight:600; ", "CCF")
#                               #         ),
#                               #         choiceValues = c("binary", "CCF"),
#                               #         selected = "binary", 
#                               #         inline = TRUE
#                               #     ),
#                               #     bsTooltip(id = "heatmaptype",
#                               #               title = ' "binary" (default) for printing a binary heatmap of mutations; or "CCF" for printing a cancer cell frequency (CCF) heatmap.',
#                               #               placement = "top",
#                               #               trigger = "hover")
#                               # ), 
#                               checkboxInput('showbootstrap',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Show bootstrap value'),value = TRUE),
#                               bsTooltip(id = "showbootstrap",
#                                         title = 'Whether to add bootstrap value on internal nodes.',
#                                         placement = "top",
#                                         trigger = "hover"),
#                               conditionalPanel(
#                                   condition = "input.showbootstrap == true",
#                                   checkboxInput('usebox',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Use box'),value = TRUE),
#                                   bsTooltip(id = "usebox",
#                                             title = 'Whether to add box arround bootstrap value on tree. ',
#                                             placement = "top",
#                                             trigger = "hover")
#                               ), 
#                               fluidRow(
#                                 column(
#                                   width = 9,
#                                   div(
#                                     tags$button(
#                                       id = "submit_phylotree", type = "button", class = "action-button bttn",
#                                       class = "bttn-unite", class = paste0("bttn-md"),
#                                       class = paste0("bttn-default"),
#                                       list(strong("Start analysis"),icon("hand-right", lib = "glyphicon")),
#                                       style = "margin-bottom:0px;margin-right:0px;"
#                                     )
#                                   )
#                                 )
#                               )
#                             )
#                           ),
#                           column(
#                             width = 9,
#                             box(
#                               width = 13,
#                               div(strong("Phylogenetic tree visualization"),style = "font-size:27px; font-weight:500;"),
#                               p("Phylogenetic tree becomes more widely used in depicting evolutionary relationships among tumors. Here, MesKit is able to plot tree-like phylogentic trees with mutational signature information, along with a binary mutation heatmap or CCF heatmap (if CCF data is available).", 
#                                 br(),"trees with mutational signature information, along with a CCF heatmap (if CCF data is available) ",
#                                 style = "font-size:20px; font-weight:500;line-height:40px;"),
#                               conditionalPanel(
#                                 condition = 'submit_phylotree',
#                                 width = NULL,
#                                 height = "100%",
#                                 uiOutput("phylotree.patientlist"),
#                                 uiOutput('warningMessage15'),
#                                 div(plotOutput("phylotree_plot",height = 600,width = 900),align = "center"),
#                                 br(),
#                                 uiOutput("phylotree_downloadbutton_ui")
#                               )
#                             )
#                           )
#                         )
# )




#Main function----
dbHeader <- dashboardHeader(title = "MesKit", titleWidth = 300, 
                            # tags$li(class = "dropdown", actionLink(inputId = "help", label = div(style = "font-size:15px; font-weight:400; ", "Help"))), 
                            # tags$li(class = "dropdown", actionLink(inputId = "contact", label = div(style = "font-size:15px; font-weight:400; ", "Contact"))),
                            dropdownMenu(
                                type = "notifications", 
                                icon = icon("question-circle"),
                                badgeStatus = NULL,
                                headerText = "Help:",
                                notificationItem("MesKit github page", icon = icon("file"),
                                                 href = "https://github.com/Niinleslie/MesKit")
                            ),
                            dropdownMenu(
                                type = "notifications", 
                                icon = icon("envelope"),
                                badgeStatus = NULL,
                                headerText = "",
                                tags$li(p("Mengni Liu, liumn5@mail2.sysu.edu.cn")),
                                tags$li(p("Chengwei Wang, wangchw8@outlook.com")),
                                tags$li(p("Jianyu Chen, chenjy327@mail2.sysu.edu.cn")),
                                tags$li(p("Xin Wang, wangx555@mail2.sysu.edu.cn"))
                                # notificationItem("Mengni Liu, liumn5@mail2.sysu.edu.cn, Sun Yat-sen university", icon = icon("user"),href = "liumn5@mail2.sysu.edu.cn"),
                                # notificationItem("Chengwei Wang, wangchw8@outlook.com, Sun Yat-sen university", icon = icon("user"),href = "wangchw8@outlook.com"),
                                # notificationItem("Jianyu Chen, chenjy327@mail2.sysu.edu.cn, Sun Yat-sen university", icon = icon("user"),href = "chenjy327@mail2.sysu.edu.cn"),
                                # notificationItem("Xin Wang, wangx555@mail2.sysu.edu.cn, Sun Yat-sen university", icon = icon("user"),href = "wangx555@mail2.sysu.edu.cn")
                            )
                            )
dbHeader$children[[2]]$children <-  tags$a(href='https://github.com/Niinleslie/MesKit',
                                           tags$img(src='image/logo.jpg',height='65',width='250'))
shinyUI(
  dashboardPage(
    skin = "blue",
    header=dbHeader ,
    sidebar=sidebar,
    body=dashboardBody(
      ## add text behind the sidebar (design error)
      tags$head(tags$style(HTML(
        "/* logo */
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
        font-family: 'Helvetica Neue',Helvetica,Arial,sans-serif;
        padding: 0 15px;
        overflow: hidden;
        color: white;
        }
        .checkbox { /* checkbox is a div class*/
        line-height: 25px;}
       input[type='checkbox']{ 
        width: 23px; /*Desired width*/
        height: 23px; /*Desired height*/
        line-height: 25px; 
      }
      span { 
            line-height: 30px; 
        }
        
        "))),
      
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
        $("header").find("nav").append(\'<span class="textnvbar"> MesKit: dissect cancer evolution from multi-region derived tumor biopsies</span>\');
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
                         .shiny-output-error-validation {
                            color: green;
                            font-size:27px; 
                            font-weight:500;
                         }
                         .main-sidebar { font-size: 16px; }


table.dataTable tbody th, table.dataTable tbody td {
    padding: 10px 1.5em !important;
}

.tooltip {
    min-width: 15em !important;
}

.progress-message, .progress-detail {
    display: block !important;
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
        bodyAL,
        bodyITH,
        bodyclone,
        # bodyfunction,
        bodytree
        # bodySurvival
      )
    )
  )
)
