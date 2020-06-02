options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(DT))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))
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
                           condition = "input.tith == 'ith_mathscore'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(),
                           checkboxInput('mathscore_withintumor',
                                         value = FALSE,
                                         label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within tumor'),
                                         width = 500),
                           bsTooltip(id = "mathscore_withintumor",
                                     title = "calculate math score within tumor",
                                     placement = "top",
                                     trigger = "hover"),
                           tags$table(
                             tags$tr(id = "inline", 
                                     width = "100%",
                                     tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                     tags$td(width = "70%", textInput(inputId = "mathscore_minvaf", value = 0.02, label = NULL)))
                           ), 
                           bsTooltip(id = "mathscore_minvaf",
                                     title = "The minimum value of vaf",
                                     placement = "top",
                                     trigger = "hover"),
                           br(),
                           fluidRow(
                             column(
                               width = 9,
                               div(
                                 tags$button(
                                   id = "submit_mathscore", type = "button", class = "action-button bttn",
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
                           condition = "input.tith == 'ith_vafcluster'",
                           div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                           br(),
                           tags$table(
                               tags$tr(id = "inline", 
                                       width = "100%",
                                       tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                       tags$td(width = "70%", textInput(inputId = "vafcluster_minvaf", value = 0.02, label = NULL)))
                           ), 
                           bsTooltip(id = "vafcluster_minvaf",
                                     title = "The minimum value of vaf",
                                     placement = "top",
                                     trigger = "hover"),
                           br(),
                           tags$table(
                               tags$tr(id = "inline", 
                                       width = "100%",
                                       tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Max vaf: ")),
                                       tags$td(width = "70%", textInput(inputId = "vafcluster_maxvaf", value = 1, label = NULL)))
                           ), 
                           bsTooltip(id = "vafcluster_maxvaf",
                                     title = "The maximum value of vaf",
                                     placement = "top",
                                     trigger = "hover"),
                           br(),
                           checkboxInput('vafcluster_withintumor',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Within tumor'),value = FALSE),
                           bsTooltip(id = "vafcluster_withintumor",
                                     title = "cluster vaf within tumors",
                                     placement = "top",
                                     trigger = "hover"),
                           fileInput(inputId = 'vafcluster_segfile', 
                                     label = div(style = "font-size:1.5em; font-weight:600; ", 'Segment file'),
                                     width = 400),
                          sliderInput(inputId='vafcluster_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 470, width = 500),
                          sliderInput(inputId='vafcluster_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 470, width = 500), 
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
                             condition = "input.tith == 'ith_ccfauc'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf: ")),
                                         tags$td(width = "70%", textInput(inputId = "ccfauc_minccf", value = 0, label = NULL)))
                             ), 
                             bsTooltip(id = "ccfauc_minccf",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             checkboxInput('ccfauc_withintumor',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within tumor'),
                                           width = 500),
                             bsTooltip(id = "ccfauc_withintumor",
                                       title = "calculate ccf within type",
                                       placement = "top",
                                       trigger = "hover"),
                             sliderInput(inputId='ccfauc_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 747, width = 500),
                             sliderInput(inputId='ccfauc_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
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
                             condition = "input.tith == 'ith_calfst'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                         tags$td(width = "70%", textInput(inputId = "calfst_minvaf", value = 0.02, label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_minvaf",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Min total depth: ")),
                                         tags$td(width = "50%", textInput(inputId = "calfst_mintotaldepth", value = 2, label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_mintotaldepth",
                                       title = "The minimum total allele depth for filtering variants. Default: 2.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             checkboxInput('calfst_withinTumor',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within Tumor'),
                                           width = 500),
                             bsTooltip(id = "calfst_withinTumor",
                                       title = "calculate fst within tumor",
                                       placement = "top",
                                       trigger = "hover"),
                             
                             checkboxInput('calfst_usecircle',
                                           value = TRUE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Use circle'),
                                           width = 500),
                             bsTooltip(id = "calfst_usecircle",
                                       title = "Logical (Default:TRUE). Whether to use circle as visualization method of correlation matrix",
                                       placement = "top",
                                       trigger = "hover"),
                             
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Title: ")),
                                         tags$td(width = "70%", textInput(inputId = "calfst_title", value = NULL, label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_title",
                                       title = "The title of the plot.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.cex: ")),
                                         tags$td(width = "70%", textInput(inputId = "calfst_numbercex", value = 8, label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_numbercex",
                                       title = "The size of text shown in correlation plot. Default 8.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),

                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.col: ")),
                                         tags$td(width = "70%", textInput(inputId = "calfst_numbercol", value = "#C77960", label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_numbercol",
                                       title = "The color of text shown in correlation plot. Default #C77960.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             sliderInput(inputId='calfst_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 560, width = 500),
                             sliderInput(inputId='calfst_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                             
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
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                         tags$td(width = "70%", textInput(inputId = "calneidist_minccf", value = 0.08, label = NULL)))
                             ), 
                             br(),
                             bsTooltip(id = "calneidist_minccf",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('calneidist_withintumor',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within Tumor'),
                                           width = 500),
                             bsTooltip(id = "calneidist_withintumor",
                                       title = "calculate ccf within type",
                                       placement = "top",
                                       trigger = "hover"),
                             
                             checkboxInput('calneidist_usecircle',
                                           value = TRUE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Use circle'),
                                           width = 500),
                             bsTooltip(id = "calneidist_usecircle",
                                       title = "Logical (Default:TRUE). Whether to use circle as visualization method of correlation matrix",
                                       placement = "top",
                                       trigger = "hover"),
                             
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Title: ")),
                                         tags$td(width = "70%", textInput(inputId = "calneidist_title", value = NULL, label = NULL)))
                             ), 
                             bsTooltip(id = "calneidist_title",
                                       title = "The title of the plot.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.cex: ")),
                                         tags$td(width = "70%", textInput(inputId = "calneidist_numbercex", value = 8, label = NULL)))
                             ), 
                             bsTooltip(id = "calneidist_numbercex",
                                       title = "The size of text shown in correlation plot. Default 8.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.col: ")),
                                         tags$td(width = "70%", textInput(inputId = "calneidist_numbercol", value = "#C77960", label = NULL)))
                             ), 
                             bsTooltip(id = "calneidist_numbercol",
                                       title = "The color of text shown in correlation plot. Default #C77960.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             
                             sliderInput(inputId='calneidist_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 560, width = 500),
                             sliderInput(inputId='calneidist_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                             
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
                           selected = "ith_mathscore",
                           side = "left",
                           tabPanel(
                             title = div(icon("chart-bar"), "MATH score"),
                             value = "ith_mathscore",
                             DT::dataTableOutput('mathScore'),
                             br(),
                             br(),
                             uiOutput("msdb")
                           ),
                           tabPanel(
                             title = div(icon("image"), "VAF clustering"),
                             value = "ith_vafcluster",
                             uiOutput("vafcluster.patientlist"),
                             uiOutput("vafcluster.samplelist"),
                             div(plotOutput("vaf",height = "100%", width = "100%"),align = "center"),
                             uiOutput("vcdb"),
                             uiOutput("vafcluster_table_ui")
                           ),
                           tabPanel(
                               title = div(icon("box"), "ccfAUC"),
                               value = "ith_ccfauc",
                               uiOutput('warningMessage_ccfauc'),
                               uiOutput('ccfauc.patientlist'),
                               div(plotOutput("ccfauc_plot",height = "100%", width = "100%"),align = "center") ,
                               uiOutput("ccfauc_db_ui"),
                               uiOutput("ccfauc_table_ui")
                           ),
                           tabPanel(
                               title = div(icon("box"), "calFst"),
                               value = "ith_calfst",
                               uiOutput('warningMessage_calfst'),
                               uiOutput('calfst.patientlist'),
                               div(plotOutput("calfst_plot",height = "100%", width = "100%"),align = "center") ,
                               uiOutput("calfst_db_ui"),
                               # uiOutput("calfst_avg_table_ui"),
                               uiOutput("calfst_pair_table_ui")
                           ),
                           tabPanel(
                               title = div(icon("box"), "calNeiDist"),
                               value = "caInput_calneidist",
                               uiOutput('calneidist.patientlist'),
                               uiOutput('warningMessage_calneidist'),
                               div(plotOutput("calneidist_plot",height = "100%", width = "100%"),align = "center") ,
                               uiOutput("calneidist_db_ui"),
                               # uiOutput("calneidist_avg_table_ui"),
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
                                  condition = "input.al_tabbox == 'pannel_plotmutprofile'",
                                  div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                  br(),
                                  uiOutput("plotmutprofile.patientlist"),
                                  selectInput("plotmutprofile_class", label = div(style = "font-size:1.5em; font-weight:600;  ", "Class"),
                                              choices = c("SP","CS","SPCS"),
                                              selected = "SP"),
                                  bsTooltip(id = "plotmutprofile_class",
                                            title = 'The class which would be represented, default is "SP" (Shared pattern: Public/Shared/Private),other options: "CS" (Clonal status: Clonal/Subclonl) and "SPCS".',
                                            placement = "top",
                                            trigger = "hover"),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Top genes count: ")),
                                              tags$td(width = "50%", textInput(inputId = "plotmutprofile_topGenesCount", value = 10, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotmutprofile_topGenesCount",
                                            title = "The number of genes print, default is 10",
                                            placement = "top",
                                            trigger = "hover"),
                                  checkboxInput('plotmutprofile_classByTumor', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Class by tumor'),value = FALSE,width = 400),
                                  bsTooltip(id = "plotmutprofile_classByTumor",
                                            title = "FALSE(Default). Define shared pattern of mutations based on tumor types (TRUE) or samples (FALSE)",
                                            placement = "top",
                                            trigger = "hover"),
                                  checkboxInput('plotmutprofile_remove_empty_columns', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Remove empty columns'),value = TRUE,width = 400),
                                  bsTooltip(id = "plotmutprofile_remove_empty_columns",
                                            title = "Whether remove the samples without alterations",
                                            placement = "top",
                                            trigger = "hover"),
                                  checkboxInput('plotmutprofile_remove_empty_rows', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Remove empty rows'),value = TRUE,width = 400),
                                  bsTooltip(id = "plotmutprofile_remove_empty_rows",
                                            title = "Whether remove the genes without alterations.",
                                            placement = "top",
                                            trigger = "hover"),
                                  checkboxInput('plotmutprofile_showColnames', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Show column names'),value = TRUE,width = 400),
                                  bsTooltip(id = "plotmutprofile_showColnames",
                                            title = "TRUE(Default). Show sample names of columns.",
                                            placement = "top",
                                            trigger = "hover"),
                                  sliderInput('plotmutprofile_width', label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'), min = 400,max = 1200, value = 900,width = 500),
                                  sliderInput('plotmutprofile_height', label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'), min = 400,max = 1200, value = 900,width = 500),
                                  fluidRow(
                                      column(
                                          width = 9,
                                          div(
                                              tags$button(
                                                  id = "submit_plotmutprofile", type = "button", class = "action-button bttn",
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
                                  fileInput(inputId = 'plotcna_segfile', 
                                            label = div(style = "font-size:1.5em; font-weight:600; ", 'Segment file'),
                                            width = 400),
                                  fileInput(inputId = 'plotcna_gisticAmpGenesFile', 
                                            label = div(style = "font-size:1.5em; font-weight:600; ", 'Gistic Amplification genes file'),
                                            width = 400),
                                  fileInput(inputId = 'plotcna_gisticDelGenesFile', 
                                            label = div(style = "font-size:1.5em; font-weight:600; ", 'Gistic deletion genes file'),
                                            width = 400),
                                  fileInput(inputId = 'plotcna_gisticAllLesionsFile', 
                                            label = div(style = "font-size:1.5em; font-weight:600; ", 'Gistic all lesions file'),
                                            width = 400),
                                  selectInput('plotcna_refBuild', label = div(style = "font-size:1.5em; font-weight:600; ", 'Genome reference'),
                                              choices = c('hg18','hg19','hg38'),selected = "hg19", width = 400),
                                  bsTooltip(id = "plotcna_refBuild",
                                            title = "Human reference genome versions of hg18, hg19 or hg38 by UCSC. Default: hg19.",
                                            placement = "top",
                                            trigger = "hover"),
                                  checkboxInput('plotcna_showGISTICgene', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Show GISTIC gene'),value = FALSE,width = 400),
                                  bsTooltip(id = "plotcna_showGISTICgene",
                                            title = "Whether GISTIC gene in seg.Default FALSE.",
                                            placement = "top",
                                            trigger = "hover"),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Gistic qvalue: ")),
                                              tags$td(width = "50%", textInput(inputId = "plotcna_gisticqval", value = 0.25, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_gisticqval",
                                            title = "The threshold of gistic Q value. Default is 0.25",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  uiOutput("plotcna.patientlist"),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Sample text size: ")),
                                              tags$td(width = "50%", textInput(inputId = "plotcna_sampletextsize", value = 11, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_sampletextsize",
                                            title = "Size of sample name.Default 11.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Legend text size: ")),
                                              tags$td(width = "50%", textInput(inputId = "plotcna_legendtextsize", value = 9, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_legendtextsize",
                                            title = "Size of legend text.Default 9.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Legend title size: ")),
                                              tags$td(width = "50%", textInput(inputId = "plotcna_legendtitlesize", value = 11, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_legendtitlesize",
                                            title = "Size of legend title.Default 11.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Sample bar size: ")),
                                              tags$td(width = "50%", textInput(inputId = "plotcna_samplebarheight", value = 0.5, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_samplebarheight",
                                            title = "Bar height of each sample .Default 0.5.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "70%", div(style = "font-size:1.5em; font-weight:600; ", "Chromosome bar size: ")),
                                              tags$td(width = "30%", textInput(inputId = "plotcna_chrombarheight", value = 0.5, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_chrombarheight",
                                            title = "Bar height of each chromosome .Default 0.5.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  sliderInput('plotcna_width', label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'), min = 400,max = 1200, value = 800,width = 500),
                                  sliderInput('plotcna_height', label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'), min = 400,max = 1200, value = 800,width = 500),
                                  br(),
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
                                  selected = 'pannel_plotmutprofile',
                                  side = 'left',
                                  height = "100%",
                                  width = "100%",
                                  tabPanel(
                                      value = 'pannel_plotmutprofile',
                                      title = div(icon("newspaper"), "Mutational landscape"), 
                                      uiOutput('warningMessage_plotmutprofile'),
                                      div(plotOutput('plotmutprofile_plot', height = "100%", width = "100%"), align = "center"),
                                      br(),
                                      uiOutput("plotmutprofile_download_button_ui")
                                  ),
                                  tabPanel(
                                      value = 'pannel_plotcna',
                                      title = div(icon("newspaper"), "CNA profile"), 
                                      uiOutput('warningMessage_plotcna'),
                                      div(plotOutput('plotcna_plot', height = "100%", width = "100%"), align = "center"),
                                      br(),
                                      uiOutput("plotcna_download_button_ui"),
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
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf: ")),
                                         tags$td(width = "70%", textInput(inputId = "compareccf_minccf", value = 0, label = NULL)))
                             ), 
                             br(),
                             bsTooltip(id = "compareccf_minccf",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('compareccf_pairbytumor',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Pair by tumor'),
                                           width = 500),
                             bsTooltip(id = "compareccf_pairbytumor",
                                       title = "Compare CCF by tumor",
                                       placement = "top",
                                       trigger = "hover"),
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
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf: ")),
                                           tags$td(width = "70%", textInput(inputId = "comparejsi_minccf", value = 0, label = NULL)))
                               ), 
                               bsTooltip(id = "comparejsi_minccf",
                                         title = "The minimum value of ccf",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               checkboxInput('comparejsi_pairbytumor',
                                             value = FALSE,
                                             label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Pair by tumor'),
                                             width = 500),
                               bsTooltip(id = "comparejsi_pairbytumor",
                                         title = "calculate JSI by tumor",
                                         placement = "top",
                                         trigger = "hover"),
                               
                               checkboxInput('comparejsi_usecircle',
                                             value = TRUE,
                                             label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Use circle'),
                                             width = 500),
                               bsTooltip(id = "comparejsi_usecircle",
                                         title = "Logical (Default:TRUE). Whether to use circle as visualization method of correlation matrix",
                                         placement = "top",
                                         trigger = "hover"),
                               
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Title: ")),
                                           tags$td(width = "70%", textInput(inputId = "comparejsi_title", value = NULL, label = NULL)))
                               ), 
                               bsTooltip(id = "comparejsi_title",
                                         title = "The title of the plot.",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.cex: ")),
                                           tags$td(width = "70%", textInput(inputId = "comparejsi_numbercex", value = 8, label = NULL)))
                               ), 
                               bsTooltip(id = "comparejsi_numbercex",
                                         title = "The size of text shown in correlation plot. Default 8.",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.col: ")),
                                           tags$td(width = "70%", textInput(inputId = "comparejsi_numbercol", value = "#C77960", label = NULL)))
                               ), 
                               bsTooltip(id = "comparejsi_numbercol",
                                         title = "The color of text shown in correlation plot. Default #C77960.",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               
                               sliderInput(inputId='comparejsi_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 560, width = 500),
                               sliderInput(inputId='comparejsi_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                               
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
                               checkboxInput('testneutral_withintumor',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Within Tumor'),value = FALSE),
                               bsTooltip(id = "testneutral_withintumor",
                                         title = 'Test neutral within tumros in each patients,default is FALSE.',
                                         placement = "top",
                                         trigger = "hover"),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Min total depth: ")),
                                           tags$td(width = "50%", textInput(inputId = "testneutral_mintotaldepth", value = 2, label = NULL)))
                               ), 
                               bsTooltip(id = "testneutral_mintotaldepth",
                                         title = "The minimun total depth of coverage. Defalut: 2",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                           tags$td(width = "70%", textInput(inputId = "testneutral_minvaf", value = 0.1, label = NULL)))
                               ), 
                               bsTooltip(id = "testneutral_minvaf",
                                         title = "The minimum value of vaf",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Max vaf: ")),
                                           tags$td(width = "70%", textInput(inputId = "testneutral_maxvaf", value = 0.3, label = NULL)))
                               ), 
                               bsTooltip(id = "testneutral_maxvaf",
                                         title = "The maximum value of vaf",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "R2 threshold: ")),
                                           tags$td(width = "50%", textInput(inputId = "testneutral_R2threshold", value = 0.98, label = NULL)))
                               ), 
                               bsTooltip(id = "testneutral_R2threshold",
                                         title = "The threshod of R2 to decide whether a tumor follows neutral evolution. Default: 0.98",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Min mut count: ")),
                                           tags$td(width = "50%", textInput(inputId = "testneutral_minmutcount", value = 20, label = NULL)))
                               ), 
                               bsTooltip(id = "testneutral_minmutcount",
                                         title = "The minimun number of subclonal mutations used to fit model. Default: 20",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               sliderInput(inputId='testneutral_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 560, width = 500),
                               sliderInput(inputId='testneutral_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                               
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
                                 div(plotOutput("comparejsi_plot",height = "100%", width = "100%"),align = "center") ,
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
                                 div(plotOutput("testneutral_plot",height = "100%", width = "100%"),align = "center") ,
                                 uiOutput("testneutral_db_ui"),
                                 uiOutput("testneutral_table_ui")
                             )
                           )
                         )
                       )
                     )
)    





bodytree <- tabItem('tree',
                         fluidRow(
                           column(
                             width = 3,
                             box(
                               width = NULL,
                               conditionalPanel(
                                   condition = "input.sgt == 'S_getphylotree'",
                                   div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                   br(),
                                   selectInput("getphylotree_method", label = div(style = "font-size:1.5em; font-weight:600;  ", "Tree construct method"),
                                               choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                                               selected = "NJ"),
                                   bsTooltip(id = "getphylotree_method",
                                             title = "Approach to construct phylogenetic trees",
                                             placement = "top",
                                             trigger = "hover"),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                               tags$td(width = "70%", textInput(inputId = "getphylotree_minvaf", value = 0.02, label = NULL)))
                                   ), 
                                   bsTooltip(id = "getphylotree_minvaf",
                                             title = "The minimum value of vaf",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf: ")),
                                               tags$td(width = "70%", textInput(inputId = "getphylotree_minccf", value = 0, label = NULL)))
                                   ), 
                                   bsTooltip(id = "getphylotree_minccf",
                                             title = "The minimum value of ccf",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "70%", div(style = "font-size:1.5em; font-weight:600; ", "Bootstrap rep num: ")),
                                               tags$td(width = "30%", textInput(inputId = "getphylotree_bootstraprepnum", value = 100, label = NULL)))
                                   ), 
                                   bsTooltip(id = "getphylotree_bootstraprepnum",
                                             title = "Bootstrap iterations. Default 100.",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   actionBttn('submit_getphylotree',div(
                                       strong("Get phylotree"),align = 'center'))
                               ), 
                               conditionalPanel(
                                   condition = "input.sgt == 'S_plotphylotree'",
                                   div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                   selectInput("plotphylotree_branchcol", label = div(style = "font-size:1.5em; font-weight:600;  ", "Branch color"),
                                               choices = c("mutType",
                                                           "mutSig",
                                                           "NULL"),
                                               selected = "mutType"),
                                   bsTooltip(id = "plotphylotree_branchcol",
                                             title = "Specify the colors of branches (Default: mutType). Other options: 'mutSig' for coloring branches by branch mutation signature;",
                                             placement = "top",
                                             trigger = "hover"),
                                   
                                   conditionalPanel(
                                       condition = "input.plotphylotree_branchcol == 'mutSig'",
                                       selectInput("plotphylotree_signatureref", label = div(style = "font-size:1.5em; font-weight:600;  ", "Signautre reference"),
                                                   choices = c("cosmic_v2",
                                                               "nature2013",
                                                               "exome_cosmic_v3"),
                                                   selected = "cosmic_v2"),
                                       bsTooltip(id = "plotphylotree_signatureref",
                                                 title = 'signature reference',
                                                 placement = "top",
                                                 trigger = "hover"),
                                       numericInput('plotphylotree_minmutcount', div(style = "font-size:1.5em; font-weight:600;  ", 'Minimal mutation number'), value = 15),
                                       bsTooltip(id = "plotphylotree_minmutcount",
                                                 title = 'The threshold for the variants in a branch. Default 15.',
                                                 placement = "top",
                                                 trigger = "hover"),
                                   ),
                                   checkboxInput('plotphylotree_showbootstrap',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Show bootstrap value'),value = TRUE),
                                   bsTooltip(id = "plotphylotree_showbootstrap",
                                             title = 'Whether to add bootstrap value on internal nodes.',
                                             placement = "top",
                                             trigger = "hover"),
                                   conditionalPanel(
                                       condition = "input.plotphylotree_showbootstrap == true",
                                       checkboxInput('plotphylotree_usebox',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Use box'),value = TRUE),
                                       bsTooltip(id = "plotphylotree_usebox",
                                                 title = 'Whether to add box arround bootstrap value on tree. ',
                                                 placement = "top",
                                                 trigger = "hover")
                                   ), 
                                   numericInput('plotphylotree_minratio', div(style = "font-size:1.5em; font-weight:600;  ", 'Min ratio'),
                                                value = 0.05, min = 0.05, max = 1, step = 0.05),
                                   bsTooltip(id = "plotphylotree_minratio",
                                             title = "Double (Default:1/20). If min.ratio is not NULL,all edge length of a phylogenetic tree should be greater than min.ratio*the longest edge length.If not, the edge length will be reset as min.ratio*longest edge length.",
                                             placement = "top",
                                             trigger = "hover"),
                                   sliderInput(inputId='plotphylotree_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 900,max = 1200, value = 1000, width = 500),
                                   sliderInput(inputId='plotphylotree_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 700,max = 1000, value = 900, width = 500),
                                   # 
                                   fluidRow(
                                       column(
                                           width = 9,
                                           div(
                                               tags$button(
                                                   id = "submit_plotphylotree", type = "button", class = "action-button bttn",
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
                                 checkboxInput('treemutsig_withintumor',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Within tumor'),value = FALSE),
                                 bsTooltip(id = "treemutsig_withintumor",
                                           title = 'Exploring signatures within tumor. Default: FALSE.',
                                           placement = "top",
                                           trigger = "hover"),
                                 selectInput("treemutsig_signatureref", label = div(style = "font-size:1.5em; font-weight:600;  ", "Signautre reference"),
                                             choices = c("cosmic_v2",
                                                         "nature2013",
                                                         "exome_cosmic_v3"),
                                             selected = "cosmic"),
                                 bsTooltip(id = "treemutsig_signatureref",
                                           title = 'signature reference',
                                           placement = "top",
                                           trigger = "hover"),
                                 numericInput('treemutsig_minmutcount', div(style = "font-size:1.5em; font-weight:600;  ", 'Minimal mutation number'), value = 15),
                                 bsTooltip(id = "treemutsig_minmutcount",
                                           title = 'The threshold for the variants in a branch. Default 15.',
                                           placement = "top",
                                           trigger = "hover"),
                                 
                                 numericInput('treemutsig_signaturecutoff', div(style = "font-size:1.5em; font-weight:600;  ", 'Signature cutoff'), min = 0,max = 1,step = 0.05,value = 0.1),
                                 bsTooltip(id = "treemutsig_signaturecutoff",
                                           title = 'The threshold for the variants in a branch. Default 15.',
                                           placement = "top",
                                           trigger = "hover"),
                                 
                                 selectInput("treemutsig_mode", label = div(style = "font-size:1.5em; font-weight:600;  ", "Mode"),
                                             choices = c('NULL',
                                                         'Original',
                                                         'Reconstructed',
                                                         'Difference'),
                                             selected = "NULL"),
                                 bsTooltip(id = "treemutsig_mode",
                                           title = "type of mutation spectrum.Default: NULL. Options:'Original','Reconstructed' or 'Difference'",
                                           placement = "top",
                                           trigger = "hover"),
                                 
                                 sliderInput(inputId='treemutsig_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 800, width = 500),
                                 sliderInput(inputId='treemutsig_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
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
                                 checkboxInput('muttrunkbranch_ct',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'CT'),value = FALSE),
                                 bsTooltip(id = "muttrunkbranch_ct",
                                           title = 'Distinction between C>T at CpG and C>T at other sites, Default FALSE',
                                           placement = "top",
                                           trigger = "hover"),
                                 numericInput(inputId = "muttrunkbranch_pvalue",
                                           label = div(style = "font-size:1.5em; font-weight:600; ", 'P-value'),
                                           min = 0,max = 1,step = 0.05,
                                           value = 0.05),
                                 bsTooltip(id = "muttrunkbranch_pvalue",
                                           title = "Confidence level of the interval for Fisher test. Default: 0.05.",
                                           placement = "right",
                                           trigger = "hover"),
                                 sliderInput(inputId='muttrunkbranch_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 470, width = 500),
                                 sliderInput(inputId='muttrunkbranch_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 470, width = 500), 
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
                                 selected = 'S_getphylotree',
                                 width = "100%",
                                 height = "100%",
                                 tabPanel(
                                     title = div(icon("newspaper"), "Get phyloTree"),
                                     value = 'S_getphylotree',
                                     # uiOutput("getphylotree.patientlist"),
                                     # uiOutput("getphylotree.slotlist"),
                                     # uiOutput("getphylotree.slotlist"),
                                 ),
                                 tabPanel(
                                     title = div(icon("newspaper"), "Plot phylotree"),
                                     value = 'S_plotphylotree',
                                     uiOutput("phylotree.patientlist"),
                                     div(plotOutput("phylotree_plot",height = "100%",width = "100%"),align = "center"),
                                     br(),
                                     uiOutput("phylotree_downloadbutton_ui")
                                 ),
                                 tabPanel(
                                     title = div(icon("microsoft"), "Signature plot"), 
                                     value = 'S_treemutsig',
                                     # uiOutput('warningMessage_treemutsig'),
                                     uiOutput("treemutsig.patientlist"),
                                     div(plotOutput('treemutsig_plot', height = "100%", width = "100%"),align = "center"),
                                     uiOutput("treemutsig_download_button_ui"),
                                     # br(),
                                     # uiOutput('treemutsig_table_ui')
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
                                 )
                               )
                             )
                           )
                         )
)




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
