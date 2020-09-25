
#required packages
suppressMessages(library(shiny))
suppressMessages(library(DT))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))
suppressMessages(library(shinyBS))
suppressMessages(library(MesKit))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))

#sider bar----

sidebar <- dashboardSidebar(
  width = 300,
  sidebarMenu(id="sidername",selected='home',
              menuItem(strong("Home"), tabName = "home", icon = icon("home")),              
              menuItem(strong("Input data"), tabName = "input", icon = icon("gear")),
              menuItem(strong("Mutational landscape"), tabName = "AL", icon = icon("bar-chart")),
              menuItem(strong("ITH evaluation"), tabName = "ITH", icon = icon("bar-chart")),
              menuItem(strong("Metastatic routes inference"), tabName = "clone", icon = icon("bar-chart")),
              menuItem(strong("PhyloTree-based analysis"), tabName = "tree", icon = icon("tree")) 
  )
)

#shinydashboar
bodyHome <- tabItem("home",
                    fluidRow(
                      box(
                        width = 12,
                        status = "info",
                        solidHeader = TRUE,
                        title = div(strong("Introduction"),style = "font-size:2em; font-weight:500;"),
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
                        title = div(strong("Overview of MesKit"),style = "font-size:2em; font-weight:500;"),
                        fluidRow(
                          column(
                              width = 6,
                                div(img(src = "image/MesKit_workflow.png", width = "90%",height = "72%"),
                                    style="text-align: center;float:left;margin:0;padding:0")
                          ),
                          column(
                              width = 5,
                              div(
                                  br(),
                                  br(),
                                  br(),
                                  br(),
                                  h3(strong("With this MesKit Shiny APP:")),
                                  p("- Identify clonal/subclonal somatic mutations",br(),
                                    "- Quantify heterogeneity within or between tumors from the same patient",br(),
                                    "- Infer metastases clonal origin",br(),
                                    "- Perform mutational signature analysis",br(),
                                    "- Construct and visualize phylogenetic trees",
                                    style = "font-size:18px; font-weight:500;line-height:50px;"),
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
                                      placeholder = "Example: HCC_LDC.maf", 
                                      width = 400),
                          fileInput(inputId = 'clinFile', 
                                    label = div(style = "font-size:1.5em; font-weight:600;",'Clinical file',
                                                tags$button(
                                                  Id = "iecontrol_clin",
                                                  type = "button",
                                                  class = "bttn-material-circle",
                                                  class = "btn action-button",
                                                  list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                  style = " background-position: center;padding:0;margin-bottom:7px;"
                                                )
                                    ), 
                                    placeholder = "Example: HCC_LDC.clin.txt", 
                                    width = 400),
                          checkboxInput('useccffile', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px",
                                                                  'CCF file',
                                                                  tags$button(
                                                                    Id = "iecontrol02",
                                                                    type = "button",
                                                                    class = "bttn-material-circle",
                                                                    class = "btn action-button",
                                                                    list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                                    style = " background-position: center;padding:0;margin-bottom:7px;"
                                                                  )),
                                        value = FALSE,width = 400),
                          bsTooltip(id = "useccffile",
                                    title = "CCF file of somatic mutations. Default NULL.",
                                    placement = "top",
                                    trigger = "hover"),
                          conditionalPanel(
                            condition = "input.useccffile == true",
                            fileInput('ccfFile',label = '',placeholder = "Example: HCC_LDC.ccf.tsv", width = 400)
                          ), 
                            selectInput('ref', label = div(style = "font-size:1.5em; font-weight:600; ", 'Genome reference'),
                                        choices = c('hg18','hg19','hg38'),selected = "hg19", width = 400),
                        bsTooltip(id = "ref",
                                    title = "Human reference genome versions of hg18,hg19 or hg38 by UCSC",
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
                            uiOutput("ie2"),
                            uiOutput("ie_clin"),
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
                           uiOutput("mathscore_patientid_ui"),
                           checkboxInput('mathscore_withintumor',
                                         value = FALSE,
                                         label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within tumor'),
                                         width = 500),
                           bsTooltip(id = "mathscore_withintumor",
                                     title = "Calculate math score within tumors in each patients",
                                     placement = "top",
                                     trigger = "hover"),
                           checkboxInput('mathscore_useadjvaf', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use adjusted vaf'),value = FALSE, width = 400),
                           bsTooltip(id = "mathscore_useadjvaf",
                                     title = "Let VAF = VAF_adj, Tumor_Average_VAF = Tumor_Average_VAF_adj.Default: FALSE. ",
                                     placement = "top",
                                     trigger = "hover"),
                           checkboxInput('mathscore_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                           bsTooltip(id = "mathscore_usetumorsamplelabel",
                                     title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                     placement = "top",
                                     trigger = "hover"),
                           tags$table(
                             tags$tr(id = "inline", 
                                     width = "100%",
                                     tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                     tags$td(width = "60%", textInput(inputId = "mathscore_minvaf", value = 0.02, label = NULL)))
                           ), 
                           bsTooltip(id = "mathscore_minvaf",
                                     title = "The minimum VAF for filtering variants. Default: 0.02. ",
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
                           uiOutput("vafcluster_patientid_ui"),
                           tags$table(
                               tags$tr(id = "inline", 
                                       width = "100%",
                                       tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                       tags$td(width = "60%", textInput(inputId = "vafcluster_minvaf", value = 0.02, label = NULL)))
                           ), 
                           bsTooltip(id = "vafcluster_minvaf",
                                     title = "The minimum value of VAF. Default: 0. Option: on the scale of 0 to 1.",
                                     placement = "top",
                                     trigger = "hover"),
                           br(),
                           tags$table(
                               tags$tr(id = "inline", 
                                       width = "100%",
                                       tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Max vaf")),
                                       tags$td(width = "60%", textInput(inputId = "vafcluster_maxvaf", value = 1, label = NULL)))
                           ), 
                           bsTooltip(id = "vafcluster_maxvaf",
                                     title = "The maximum value of VAF. Default: 0. Option: on the scale of 0 to 1.",
                                     placement = "top",
                                     trigger = "hover"),
                           br(),
                           checkboxInput('vafcluster_withintumor',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Within tumor'),value = FALSE),
                           bsTooltip(id = "vafcluster_withintumor",
                                     title = "Cluster VAF within tumors in each patients,default is FALSE.",
                                     placement = "top",
                                     trigger = "hover"),
                           checkboxInput('vafcluster_useadjvaf', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use adjusted vaf'),value = FALSE, width = 400),
                           bsTooltip(id = "vafcluster_useadjvaf",
                                     title = "Let VAF = VAF_adj, Tumor_Average_VAF = Tumor_Average_VAF_adj.Default: FALSE. ",
                                     placement = "top",
                                     trigger = "hover"),
                           checkboxInput('vafcluster_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                           bsTooltip(id = "vafcluster_usetumorsamplelabel",
                                     title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
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
                             uiOutput("ccfauc_patientid_ui"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                         tags$td(width = "60%", textInput(inputId = "ccfauc_minccf", value = 0, label = NULL)))
                             ), 
                             bsTooltip(id = "ccfauc_minccf",
                                       title = "The minimum value of CCF. Default: 0.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             checkboxInput('ccfauc_withintumor',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within tumor'),
                                           width = 500),
                             bsTooltip(id = "ccfauc_withintumor",
                                       title = "Calculate AUC within tumors in each patients, default is FALSE.",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('ccfauc_useadjvaf', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use adjusted vaf'),value = FALSE, width = 400),
                             bsTooltip(id = "ccfauc_useadjvaf",
                                       title = "Let VAF = VAF_adj, Tumor_Average_VAF = Tumor_Average_VAF_adj.Default: FALSE. ",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('ccfauc_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                             bsTooltip(id = "ccfauc_usetumorsamplelabel",
                                       title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
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
                             br(),
                             uiOutput("calfst_patientid_ui"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                         tags$td(width = "60%", textInput(inputId = "calfst_minvaf", value = 0.02, label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_minvaf",
                                       title = "Specify the minimum VAF_adj, default is 0.02.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Min total depth")),
                                         tags$td(width = "40%", textInput(inputId = "calfst_mintotaldepth", value = 2, label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_mintotaldepth",
                                       title = "The minimum total allele depth for filtering variants. Default: 2.",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('calfst_withinTumor',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within Tumor'),
                                           width = 500),
                             bsTooltip(id = "calfst_withinTumor",
                                       title = "Calculate fst within types in each patients,default is FALSE.",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('calfst_useadjvaf', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use adjusted vaf'),value = FALSE, width = 400),
                             bsTooltip(id = "calfst_useadjvaf",
                                       title = "Let VAF = VAF_adj, Tumor_Average_VAF = Tumor_Average_VAF_adj.Default: FALSE. ",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('calfst_usecircle',
                                           value = TRUE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Use circle'),
                                           width = 500),
                             checkboxInput('calfst_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                             bsTooltip(id = "calfst_usetumorsamplelabel",
                                       title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                       placement = "top",
                                       trigger = "hover"),
                             bsTooltip(id = "calfst_usecircle",
                                       title = "Logical (Default:TRUE). Whether to use circle as visualization method of correlation matrix",
                                       placement = "top",
                                       trigger = "hover"),
                             
                             # tags$table(
                             #     tags$tr(id = "inline", 
                             #             width = "100%",
                             #             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Title: ")),
                             #             tags$td(width = "70%", textInput(inputId = "calfst_title", value = NULL, label = NULL)))
                             # ), 
                             # bsTooltip(id = "calfst_title",
                             #           title = "The title of the plot.",
                             #           placement = "top",
                             #           trigger = "hover"),
                             # br(),
                             
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.cex")),
                                         tags$td(width = "50%", textInput(inputId = "calfst_numbercex", value = 8, label = NULL)))
                             ), 
                             bsTooltip(id = "calfst_numbercex",
                                       title = "The size of text shown in correlation plot. Default 8.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),

                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.col")),
                                         tags$td(width = "50%", textInput(inputId = "calfst_numbercol", value = "#C77960", label = NULL)))
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
                             condition = "input.tith == 'ith_mutheatmap'",
                             div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                             br(),
                             uiOutput("mutheatmap_patientid_ui"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                         tags$td(width = "60%", textInput(inputId = "mutheatmap_minvaf", value = 0.02, label = NULL)))
                             ), 
                             bsTooltip(id = "mutheatmap_minvaf",
                                       title = "The minimum value of VAF. Default: 0.02. Option: on the scale of 0 to 1.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                         tags$td(width = "60%", textInput(inputId = "mutheatmap_minccf", value = 0, label = NULL)))
                             ), 
                             bsTooltip(id = "mutheatmap_minccf",
                                       title = "The minimum value of CCF. Default: 0.02. Option: on the scale of 0 to 1.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             checkboxInput('mutheatmap_useccf',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Use ccf'),
                                           width = 500),
                             bsTooltip(id = "mutheatmap_useccf",
                                       title = "Logical. If FALSE (default), print a binary heatmap of mutations. Otherwise, print a cancer cell frequency (CCF) heatmap.",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('mutheatmap_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                             bsTooltip(id = "mutheatmap_usetumorsamplelabel",
                                       title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                       placement = "top",
                                       trigger = "hover"),
                             # bsTooltip(id = "mutheatmap_useccf",
                             #           title = "Logical. If FALSE (default), print a binary heatmap of mutations. Otherwise, print a cancer cell frequency (CCF) heatmap.",
                             #           placement = "top",
                             #           trigger = "hover"),
                             # fileInput(inputId = 'mutheatmap_genelist', 
                             #           label = div(style = "font-size:1.5em; font-weight:600; ", 'Gene list file'),
                             #           placeholder = "Default: IntOGen-DriverGenes_HC.tsv",
                             #           width = 400),
                             uiOutput("mutheatmap_parameters_ui"),
                             # checkboxInput('mutheatmap_plotgenelist',
                             #               value = FALSE,
                             #               label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Plot gene list'),
                             #               width = 500),
                             # bsTooltip(id = "mutheatmap_plotgenelist",
                             #           title = "If TRUE, plot heatmap with genes on geneList when geneList is not NULL.Default FALSE.",
                             #           placement = "top",
                             #           trigger = "hover"),
                             # 
                             # checkboxInput('mutheatmap_showgene',
                             #               value = FALSE,
                             #               label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Show gene'),
                             #               width = 500),
                             # bsTooltip(id = "mutheatmap_showgene",
                             #           title = "Show the name of genes next to the heatmap.Default FALSE.",
                             #           placement = "top",
                             #           trigger = "hover"),
                             # 
                             # checkboxInput('mutheatmap_showgenelist',
                             #               value = TRUE,
                             #               label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Show gene list'),
                             #               width = 500),
                             # bsTooltip(id = "mutheatmap_showgenelist",
                             #           title = "Show the names of gene on the geneList.Default TRUE.",
                             #           placement = "top",
                             #           trigger = "hover"),
                             
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "60%", div(style = "font-size:1.5em; font-weight:600; ", "Mutation threshold")),
                                         tags$td(width = "30%", textInput(inputId = "mutheatmap_mutthreshold", value = 150, label = NULL)))
                             ), 
                             bsTooltip(id = "mutheatmap_mutthreshold",
                                       title = "show.gene and show.geneList will be FALSE when patient have more mutations than threshold.Default is 150.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Sample text size")),
                                         tags$td(width = "40%", textInput(inputId = "mutheatmap_sampletextsize", value = 9, label = NULL)))
                             ), 
                             bsTooltip(id = "mutheatmap_sampletextsize",
                                       title = "Size of sample name.Default 9.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Legend title size")),
                                         tags$td(width = "40%", textInput(inputId = "mutheatmap_legendtitlesize", value = 10, label = NULL)))
                             ), 
                             bsTooltip(id = "mutheatmap_legendtitlesize",
                                       title = "Size of legend title.Default 9.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             sliderInput(inputId='mutheatmap_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 400,max = 1000, value = 700, width = 500),
                             sliderInput(inputId='mutheatmap_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 560, width = 500), 
                             
                             br(),
                             fluidRow(
                                 column(
                                     width = 9,
                                     div(
                                         tags$button(
                                             id = "submit_mutheatmap", type = "button", class = "action-button bttn",
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
                             uiOutput("calneidist_patientid_ui"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                         tags$td(width = "60%", textInput(inputId = "calneidist_minccf", value = 0.08, label = NULL)))
                             ), 
                             br(),
                             bsTooltip(id = "calneidist_minccf",
                                       title = "Specify the minimum CCF, default is 0.08.",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('calneidist_withintumor',
                                           value = FALSE,
                                           label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Within Tumor'),
                                           width = 500),
                             bsTooltip(id = "calneidist_withintumor",
                                       title = "Calculate fst within tumors in each patients,default is FALSE.",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('calneidist_useadjvaf', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use adjusted vaf'),value = FALSE, width = 400),
                             bsTooltip(id = "calneidist_useadjvaf",
                                       title = "Let VAF = VAF_adj, Tumor_Average_VAF = Tumor_Average_VAF_adj.Default: FALSE. ",
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
                             checkboxInput('calneidist_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                             bsTooltip(id = "calneidist_usetumorsamplelabel",
                                       title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                       placement = "top",
                                       trigger = "hover"),
                             
                             # tags$table(
                             #     tags$tr(id = "inline", 
                             #             width = "100%",
                             #             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Title: ")),
                             #             tags$td(width = "70%", textInput(inputId = "calneidist_title", value = NULL, label = NULL)))
                             # ), 
                             # bsTooltip(id = "calneidist_title",
                             #           title = "The title of the plot.",
                             #           placement = "top",
                             #           trigger = "hover"),
                             # br(),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.cex")),
                                         tags$td(width = "50%", textInput(inputId = "calneidist_numbercex", value = 8, label = NULL)))
                             ), 
                             bsTooltip(id = "calneidist_numbercex",
                                       title = "The size of text shown in correlation plot. Default 8.",
                                       placement = "top",
                                       trigger = "hover"),
                             br(),
                             
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.col")),
                                         tags$td(width = "50%", textInput(inputId = "calneidist_numbercol", value = "#C77960", label = NULL)))
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
                         p("Understanding the origin and development of intra-tumor heterogeneity is clinically important, which has the potential to yield insights to guide therapeutic strategies. MesKit has integrated several metrics to estimate ITH within region or between regions borrowed from published research and population genetics.",
                           style = "font-size:20px; font-weight:500;line-height:40px;"),
                         tabBox(
                           id = 'tith',
                           height = "100%", 
                           width = "100%",
                           selected = "ith_mathscore",
                           side = "left",
                           tabPanel(
                             title = div(icon("chart-bar"), "MATH score", style = "font-size:1.5em; font-weight:600; "),
                             value = "ith_mathscore",
                             DT::dataTableOutput('mathScore'),
                             br(),
                             br(),
                             uiOutput("msdb")
                           ),
                           tabPanel(
                             title = div(icon("image"), "VAF clustering", style = "font-size:1.5em; font-weight:600; "),
                             value = "ith_vafcluster",
                             uiOutput("vafcluster.patientlist"),
                             uiOutput("vafcluster_table_ui"),
                             uiOutput("vafcluster.samplelist"),
                             div(plotOutput("vaf",height = "100%", width = "100%"),align = "left"),
                             br(),
                             uiOutput("vcdb")
                           ),
                           tabPanel(
                               title = div(icon("image"), "AUC of CCF", style = "font-size:1.5em; font-weight:600; "),
                               value = "ith_ccfauc",
                               uiOutput('ccfauc.patientlist'),
                               uiOutput("ccfauc_table_ui"),
                               div(plotOutput("ccfauc_plot",height = "100%", width = "100%"),align = "left") ,
                               uiOutput("ccfauc_db_ui")
                           ),
                           tabPanel(
                               title = div(icon("image"), "Fixation index", style = "font-size:1.5em; font-weight:600; "),
                               value = "ith_calfst",
                               uiOutput('calfst.patientlist'),
                               uiOutput("calfst_pair_table_ui"),
                               div(plotOutput("calfst_plot",height = "100%", width = "100%"),align = "left") ,
                               uiOutput("calfst_db_ui")
                               # uiOutput("calfst_avg_table_ui"),
                           ),
                           tabPanel(
                               title = div(icon("image"), "Nei's distance", style = "font-size:1.5em; font-weight:600; "),
                               value = "caInput_calneidist",
                               uiOutput('calneidist.patientlist'),
                               
                               uiOutput("calneidist_pair_table_ui"),
                               
                               div(plotOutput("calneidist_plot",height = "100%", width = "100%"),align = "left") ,
                               uiOutput("calneidist_db_ui")
                               # uiOutput("calneidist_avg_table_ui"),
                           ),
                           tabPanel(
                               title = div(icon("image"), "Heatmap", style = "font-size:1.5em; font-weight:600; "),
                               value = "ith_mutheatmap",
                               uiOutput('mutheatmap.patientlist'),
                               div(plotOutput("mutheatmap_plot",height = "100%", width = "100%"),align = "left") ,
                               uiOutput("mutheatmap_db_ui")
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
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Top genes count")),
                                              tags$td(width = "40%", textInput(inputId = "plotmutprofile_topGenesCount", value = 10, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotmutprofile_topGenesCount",
                                            title = "The number of genes print, default is 10",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  checkboxInput('plotmutprofile_usegenelist', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Gene list'),value = FALSE,width = 400),
                                  bsTooltip(id = "plotmutprofile_usegenelist",
                                            title = "A list of genes to restrict the analysis. Default NULL.",
                                            placement = "top",
                                            trigger = "hover"),
                                  conditionalPanel(
                                    condition = "input.plotmutprofile_usegenelist == true",
                                    fileInput(inputId = 'plotmutprofile_genelist', 
                                              label = div(style = "font-size:1.5em; font-weight:600; ", 'Gene list file'),
                                              placeholder = "Example: IntOGen-DriverGenes_HC.tsv",
                                              width = 400)
                                    ), 
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
                                  checkboxInput('plotmutprofile_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                                  bsTooltip(id = "plotmutprofile_useadjvaf",
                                            title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
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
                                            label = div(style = "font-size:1.5em; font-weight:600;",'Segment file',
                                                        tags$button(
                                                          Id = "iecontrol_seg",
                                                          type = "button",
                                                          class = "bttn-material-circle",
                                                          class = "btn action-button",
                                                          list(tags$img(src = "image/button.png",width = "22px",height = "22px")),
                                                          style = " background-position: center;padding:0;margin-bottom:7px;"
                                                        )
                                            ), 
                                            placeholder = "Example: HCC_LDC.seg.txt", 
                                            width = 400),
                                  # fileInput(inputId = 'plotcna_segfile', 
                                  #           label = div(style = "font-size:1.5em; font-weight:600; ",
                                  #                       'Segment file'),
                                  #           width = 400),
                                  checkboxInput('plotmutprofile_usegisticAmpGenes', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Gistic amplification genes'),value = FALSE,width = 400),
                                  bsTooltip(id = "plotmutprofile_usegisticAmpGenes",
                                            title = "Amplification Genes file generated by GISTIC.",
                                            placement = "top",
                                            trigger = "hover"),
                                  conditionalPanel(
                                    condition = "input.plotmutprofile_usegisticAmpGenes == true",
                                    fileInput(inputId = 'plotcna_gisticAmpGenesFile', 
                                              label = '',
                                              placeholder = "Example: LIHC_amp_genes.conf_99.txt",
                                              width = 400)
                                  ), 
                                  checkboxInput('plotmutprofile_usegisticDelGenes', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Gistic deletion genes'),value = FALSE,width = 400),
                                  bsTooltip(id = "plotmutprofile_usegisticDelGenes",
                                            title = "Deletion Genes file generated by GISTIC.",
                                            placement = "top",
                                            trigger = "hover"),
                                  conditionalPanel(
                                    condition = "input.plotmutprofile_usegisticDelGenes == true",
                                    fileInput(inputId = 'plotcna_gisticDelGenesFile', 
                                              label = '',
                                              placeholder = "Example: LIHC_del_genes.conf_99.txt",
                                              width = 400)
                                  ), 
                                  checkboxInput('plotmutprofile_usegisticAllLesions', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Gistic all lesions'),value = FALSE,width = 400),
                                  bsTooltip(id = "plotmutprofile_usegisticAllLesions",
                                            title = "Information of all lesions generated by GISTIC.",
                                            placement = "top",
                                            trigger = "hover"),
                                  conditionalPanel(
                                    condition = "input.plotmutprofile_usegisticAllLesions == true",
                                    fileInput(inputId = 'plotcna_gisticAllLesionsFile', 
                                              label = '',
                                              placeholder = "Example: LIHC_all_lesions.conf_99.txt",
                                              
                                              width = 400)
                                  ), 
                                  selectInput('plotcna_refBuild', label = div(style = "font-size:1.5em; font-weight:600; ", 'Genome reference'),
                                              choices = c('hg18','hg19','hg38'),selected = "hg19", width = 400),
                                  bsTooltip(id = "plotcna_refBuild",
                                            title = "Human reference genome versions of hg18, hg19 or hg38 by UCSC. Default: hg19.",
                                            placement = "top",
                                            trigger = "hover"),
                                  uiOutput("plotcna_gistic_parameters_ui"),
                                  br(),
                                  uiOutput("plotcna.patientlist"),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Sample text size")),
                                              tags$td(width = "40%", textInput(inputId = "plotcna_sampletextsize", value = 11, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_sampletextsize",
                                            title = "Size of sample name.Default 11.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Legend text size")),
                                              tags$td(width = "40%", textInput(inputId = "plotcna_legendtextsize", value = 9, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_legendtextsize",
                                            title = "Size of legend text.Default 9.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Legend title size")),
                                              tags$td(width = "40%", textInput(inputId = "plotcna_legendtitlesize", value = 11, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_legendtitlesize",
                                            title = "Size of legend title.Default 11.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Sample bar size")),
                                              tags$td(width = "40%", textInput(inputId = "plotcna_samplebarheight", value = 0.5, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_samplebarheight",
                                            title = "Bar height of each sample.Default 0.5.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  tags$table(
                                      tags$tr(id = "inline", 
                                              width = "100%",
                                              tags$td(width = "70%", div(style = "font-size:1.5em; font-weight:600; ", "Chromosome bar size")),
                                              tags$td(width = "30%", textInput(inputId = "plotcna_chrombarheight", value = 0.5, label = NULL)))
                                  ), 
                                  bsTooltip(id = "plotcna_chrombarheight",
                                            title = "Bar height of each chromosome.Default 0.5.",
                                            placement = "top",
                                            trigger = "hover"),
                                  br(),
                                  checkboxInput('plotcna_showrownames',label = div(style = "font-size:1.5em; font-weight:600;padding-left:15px ", 'Show row names'),value = TRUE),
                                  checkboxInput('plotcna_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                                  bsTooltip(id = "plotcna_usetumorsamplelabel",
                                            title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                            placement = "top",
                                            trigger = "hover"),
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
                              div(strong("Mutational landscape"),style = "font-size:27px; font-weight:500;"),
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
                                      title = div(icon("image"), "Mutational profile",style = "font-size:1.5em; font-weight:600; "), 
                                      div(plotOutput('plotmutprofile_plot', height = "100%", width = "100%"), align = "left"),
                                      br(),
                                      uiOutput("plotmutprofile_download_button_ui")
                                  ),
                                  tabPanel(
                                      value = 'pannel_plotcna',
                                      title = div(icon("image"), "CNA profile", style = "font-size:1.5em; font-weight:600; "),
                                      uiOutput("ie_seg"),
                                      uiOutput("plotcna_table_ui"),
                                      div(plotOutput('plotcna_plot', height = "100%", width = "100%"), align = "left"),
                                      br(),
                                      uiOutput("plotcna_download_button_ui")
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
                             uiOutput("compareccf_patientid_ui"),
                             tags$table(
                                 tags$tr(id = "inline", 
                                         width = "100%",
                                         tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                         tags$td(width = "60%", textInput(inputId = "compareccf_minccf", value = 0, label = NULL)))
                             ), 
                             br(),
                             bsTooltip(id = "compareccf_minccf",
                                       title = "The minimum value of ccf",
                                       placement = "top",
                                       trigger = "hover"),
                             checkboxInput('compareccf_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                             bsTooltip(id = "compareccf_usetumorsamplelabel",
                                       title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                       placement = "top",
                                       trigger = "hover"),
                             uiOutput("compareccf_pairbytumor_ui"),
                             # checkboxInput('compareccf_pairbytumor',
                             #               value = FALSE,
                             #               label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Pair by tumor'),
                             #               width = 500),
                             # bsTooltip(id = "compareccf_pairbytumor",
                             #           title = "Compare CCF by tumor",
                             #           placement = "top",
                             #           trigger = "hover"),
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
                               uiOutput("comparejsi_patientid_ui"),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                           tags$td(width = "60%", textInput(inputId = "comparejsi_minccf", value = 0, label = NULL)))
                               ), 
                               bsTooltip(id = "comparejsi_minccf",
                                         title = "The minimum value of ccf",
                                         placement = "top",
                                         trigger = "hover"),
                               checkboxInput('comparejsi_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                               bsTooltip(id = "comparejsi_usetumorsamplelabel",
                                         title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               uiOutput("comparejsi_pairbytumor_ui"),
                               # checkboxInput('comparejsi_pairbytumor',
                               #               value = FALSE,
                               #               label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Pair by tumor'),
                               #               width = 500),
                               # bsTooltip(id = "comparejsi_pairbytumor",
                               #           title = "calculate JSI by tumor",
                               #           placement = "top",
                               #           trigger = "hover"),
                               
                               checkboxInput('comparejsi_usecircle',
                                             value = TRUE,
                                             label = div(style = "font-size:1.5em; font-weight:600; padding-left:12px", 'Use circle'),
                                             width = 500),
                               bsTooltip(id = "comparejsi_usecircle",
                                         title = "Logical (Default:TRUE). Whether to use circle as visualization method of correlation matrix",
                                         placement = "top",
                                         trigger = "hover"),
                               
                               # tags$table(
                               #     tags$tr(id = "inline", 
                               #             width = "100%",
                               #             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Title: ")),
                               #             tags$td(width = "70%", textInput(inputId = "comparejsi_title", value = NULL, label = NULL)))
                               # ), 
                               # bsTooltip(id = "comparejsi_title",
                               #           title = "The title of the plot.",
                               #           placement = "top",
                               #           trigger = "hover"),
                               # br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.cex")),
                                           tags$td(width = "40%", textInput(inputId = "comparejsi_numbercex", value = 8, label = NULL)))
                               ), 
                               bsTooltip(id = "comparejsi_numbercex",
                                         title = "The size of text shown in correlation plot. Default 8.",
                                         placement = "top",
                                         trigger = "hover"),
                               br(),
                               tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "number.col")),
                                           tags$td(width = "40%", textInput(inputId = "comparejsi_numbercol", value = "#C77960", label = NULL)))
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
                               br(),
                               uiOutput("testneutral_patientid_ui"),
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
                           div(strong("Metastatic routes inference"),style = "font-size:27px; font-weight:500;"),
                           p("Since metastasis is the ultimate cause of death for most patients, it is particularly important to gain a systematic understanding of how tumor disseminates and the scale of ongoing parallel evolution in metastatic and primary site. Here, we provide two functions to help distinguish monoclonal from polyclonal seeding. ",
                             style = "font-size:20px; font-weight:500;line-height:40px;"),
                           tabBox(
                             id = 'clt',
                             selected = 'clone_compareccf',
                             side = 'left',
                             height = "100%",
                             width = "100%",
                             tabPanel(
                               value = 'clone_compareccf',
                               title = div(icon("chart-bar"), "CCF comparison", style = "font-size:1.5em; font-weight:600; "), 
                               uiOutput('compareccf.patientlist'),
                               uiOutput('compareccf.samplelist'), 
                               br(),
                               uiOutput("compareccf_table_ui")
                             ),
                             tabPanel(
                                 title = div(icon("image"), "Jaccard similarity index", style = "font-size:1.5em; font-weight:600; "),
                                 value = "clone_comparejsi",
                                 uiOutput('comparejsi.patientlist'),
                                 uiOutput("comparejsi_pair_table_ui"),
                                 br(),
                                 div(plotOutput("comparejsi_plot",height = "100%", width = "100%"),align = "left") ,
                                 uiOutput("comparejsi_db_ui")
                             )
                             # tabPanel(
                             #     title = div(icon("box"), "testNeutral"),
                             #     value = "clone_testneutral",
                             #     uiOutput('testneutral.patientlist'),
                             #     uiOutput("testneutral_table_ui"),
                             #     br(),
                             #     uiOutput('testneutral.samplelist'),
                             #     uiOutput('warningMessage_testneutral'),
                             #     div(plotOutput("testneutral_plot",height = "100%", width = "100%"),align = "center") ,
                             #     uiOutput("testneutral_db_ui")
                             # )
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
                               # conditionalPanel(
                               #     condition = "input.sgt == 'S_getphylotree'",
                               #     div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                               #     br(),
                               #     selectInput("getphylotree_method", label = div(style = "font-size:1.5em; font-weight:600;  ", "Tree construct method"),
                               #                 choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                               #                 selected = "NJ"),
                               #     bsTooltip(id = "getphylotree_method",
                               #               title = "Approach to construct phylogenetic trees",
                               #               placement = "top",
                               #               trigger = "hover"),
                               #     tags$table(
                               #         tags$tr(id = "inline", 
                               #                 width = "100%",
                               #                 tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                               #                 tags$td(width = "70%", textInput(inputId = "getphylotree_minvaf", value = 0.02, label = NULL)))
                               #     ), 
                               #     bsTooltip(id = "getphylotree_minvaf",
                               #               title = "The minimum value of vaf",
                               #               placement = "top",
                               #               trigger = "hover"),
                               #     br(),
                               #     tags$table(
                               #         tags$tr(id = "inline", 
                               #                 width = "100%",
                               #                 tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf: ")),
                               #                 tags$td(width = "70%", textInput(inputId = "getphylotree_minccf", value = 0, label = NULL)))
                               #     ), 
                               #     bsTooltip(id = "getphylotree_minccf",
                               #               title = "The minimum value of ccf",
                               #               placement = "top",
                               #               trigger = "hover"),
                               #     br(),
                               #     tags$table(
                               #         tags$tr(id = "inline", 
                               #                 width = "100%",
                               #                 tags$td(width = "60%", div(style = "font-size:1.5em; font-weight:600; ", "Boostrap repetitions: ")),
                               #                 tags$td(width = "40%", textInput(inputId = "getphylotree_bootstraprepnum", value = 100, label = NULL)))
                               #     ), 
                               #     bsTooltip(id = "getphylotree_bootstraprepnum",
                               #               title = "Bootstrap iterations. Default 100.",
                               #               placement = "top",
                               #               trigger = "hover"),
                               #     br(),
                               #     actionBttn('submit_getphylotree',div(
                               #         strong("Get phylotree"),align = 'center'))
                               # ), 
                               conditionalPanel(
                                   condition = "input.sgt == 'S_plotphylotree'",
                                   div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                   br(),
                                   uiOutput("plotphylotree_patientid_ui"),
                                   selectInput("plotphylotree_getphylotree_method", label = div(style = "font-size:1.5em; font-weight:600;  ", "Tree construct method"),
                                               choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                                               selected = "NJ"),
                                   bsTooltip(id = "plotphylotree_getphylotree_method",
                                             title = "Approach to construct phylogenetic trees",
                                             placement = "top",
                                             trigger = "hover"),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                               tags$td(width = "60%", textInput(inputId = "plotphylotree_getphylotree_minvaf", value = 0.02, label = NULL)))
                                   ), 
                                   bsTooltip(id = "plotphylotree_getphylotree_minvaf",
                                             title = "The minimum value of vaf",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                               tags$td(width = "60%", textInput(inputId = "plotphylotree_getphylotree_minccf", value = 0, label = NULL)))
                                   ), 
                                   bsTooltip(id = "plotphylotree_getphylotree_minccf",
                                             title = "The minimum value of ccf",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Boostrap repetitions")),
                                               tags$td(width = "30%", textInput(inputId = "plotphylotree_getphylotree_bootstraprepnum", value = 100, label = NULL)))
                                   ), 
                                   bsTooltip(id = "plotphylotree_getphylotree_bootstraprepnum",
                                             title = "Bootstrap iterations. Default 100.",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
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
                                       tags$table(
                                         tags$tr(id = "inline", 
                                                 width = "100%",
                                                 tags$td(width = "60%", div(style = "font-size:1.5em; font-weight:600; ", "Minimal mutation number")),
                                                 tags$td(width = "20%", textInput(inputId = "plotphylotree_minmutcount", value = 15, label = NULL)))
                                       ), 
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
                                   checkboxInput('plotphylotree_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                                   bsTooltip(id = "plotphylotree_usetumorsamplelabel",
                                             title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                             placement = "top",
                                             trigger = "hover"),
                                   # conditionalPanel(
                                   #     condition = "input.plotphylotree_showbootstrap == true",
                                   #     checkboxInput('plotphylotree_usebox',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'Use box'),value = TRUE),
                                   #     bsTooltip(id = "plotphylotree_usebox",
                                   #               title = 'Whether to add box arround bootstrap value on tree. ',
                                   #               placement = "top",
                                   #               trigger = "hover")
                                   # ), 
                                   tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ratio")),
                                             tags$td(width = "60%", textInput(inputId = "plotphylotree_minratio", value = 0.05, label = NULL)))
                                   ), 
                                   bsTooltip(id = "plotphylotree_minratio",
                                             title = "Double (Default:1/20). If min.ratio is not NULL,all edge length of a phylogenetic tree should be greater than min.ratio*the longest edge length.If not, the edge length will be reset as min.ratio*longest edge length.",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   sliderInput(inputId='plotphylotree_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 300,max = 1200, value = 500, width = 500),
                                   sliderInput(inputId='plotphylotree_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 300,max = 1200, value = 500, width = 500),
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
                                   condition = "input.sgt == 'S_comparetree'",
                                   div(strong("Parameter"),style = "font-size:2em; font-weight:600;"),
                                   br(),
                                   uiOutput("comparetree.patientlist"),
                                   selectInput("comparetree_getphylotree_method1",
                                               label = div(style = "font-size:1.5em; font-weight:600;  ",
                                                           "Tree1 construct method"),
                                               choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                                               selected = "NJ"),
                                   bsTooltip(id = "comparetree_getphylotree_method1",
                                             title = "Approach to construct phylogenetic trees",
                                             placement = "top",
                                             trigger = "hover"),
                                   selectInput("comparetree_getphylotree_method2", 
                                               label = div(style = "font-size:1.5em; font-weight:600;  ",
                                                           "Tree2 construct method"),
                                               choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                                               selected = "MP"),
                                   bsTooltip(id = "comparetree_getphylotree_method2",
                                             title = "Approach to construct phylogenetic trees",
                                             placement = "top",
                                             trigger = "hover"),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                               tags$td(width = "60%", textInput(inputId = "comparetree_getphylotree_minvaf", value = 0.02, label = NULL)))
                                   ), 
                                   bsTooltip(id = "comparetree_getphylotree_minvaf",
                                             title = "The minimum value of vaf",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                               tags$td(width = "60%", textInput(inputId = "comparetree_getphylotree_minccf", value = 0, label = NULL)))
                                   ), 
                                   bsTooltip(id = "comparetree_getphylotree_minccf",
                                             title = "The minimum value of ccf",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   tags$table(
                                       tags$tr(id = "inline", 
                                               width = "100%",
                                               tags$td(width = "50%", 
                                                       div(style = "font-size:1.5em; font-weight:600; ",
                                                           "Boostrap repetitions")),
                                               tags$td(width = "30%", 
                                                       textInput(inputId = "comparetree_getphylotree_bootstraprepnum", 
                                                                 value = 100, label = NULL)))
                                   ), 
                                   bsTooltip(id = "comparetree_getphylotree_bootstraprepnum",
                                             title = "Bootstrap iterations. Default 100.",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   # div(strong("Parameter(phyloTree2)"),style = "font-size:1.6em; font-weight:600;"),
                                   # selectInput("comparetree_getphylotree_method2", label = div(style = "font-size:1.5em; font-weight:600;  ", "Tree construct method"),
                                   #             choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                                   #             selected = "MP"),
                                   # bsTooltip(id = "comparetree_getphylotree_method2",
                                   #           title = "Approach to construct phylogenetic trees",
                                   #           placement = "top",
                                   #           trigger = "hover"),
                                   # tags$table(
                                   #     tags$tr(id = "inline", 
                                   #             width = "100%",
                                   #             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf: ")),
                                   #             tags$td(width = "70%", textInput(inputId = "comparetree_getphylotree_minvaf2", value = 0.02, label = NULL)))
                                   # ), 
                                   # bsTooltip(id = "comparetree_getphylotree_minvaf2",
                                   #           title = "The minimum value of vaf",
                                   #           placement = "top",
                                   #           trigger = "hover"),
                                   # br(),
                                   # tags$table(
                                   #     tags$tr(id = "inline", 
                                   #             width = "100%",
                                   #             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf: ")),
                                   #             tags$td(width = "70%", textInput(inputId = "comparetree_getphylotree_minccf2", value = 0, label = NULL)))
                                   # ), 
                                   # bsTooltip(id = "comparetree_getphylotree_minccf2",
                                   #           title = "The minimum value of ccf",
                                   #           placement = "top",
                                   #           trigger = "hover"),
                                   # br(),
                                   # tags$table(
                                   #     tags$tr(id = "inline", 
                                   #             width = "100%",
                                   #             tags$td(width = "60%", div(style = "font-size:1.5em; font-weight:600; ", "Boostrap repetitions: ")),
                                   #             tags$td(width = "40%", textInput(inputId = "comparetree_getphylotree_bootstraprepnum2", value = 100, label = NULL)))
                                   # ), 
                                   # bsTooltip(id = "comparetree_getphylotree_bootstraprepnum2",
                                   #           title = "Bootstrap iterations. Default 100.",
                                   #           placement = "top",
                                   #           trigger = "hover"),
                                   # br(),
                                   # div(strong("Parameter(compareTree)"),style = "font-size:1.6em; font-weight:600;"),
                                   checkboxInput('comparetree_showbootstrap',
                                                 div(style = "font-size:1.5em; font-weight:600; padding-left:15px ",
                                                     'Show bootstrap value'),value = FALSE),
                                   bsTooltip(id = "comparetree_showbootstrap",
                                             title = 'Whether to add bootstrap value on internal nodes.',
                                             placement = "top",
                                             trigger = "hover"),
                                   checkboxInput('comparetree_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                                   bsTooltip(id = "comparetree_usetumorsamplelabel",
                                             title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
                                             placement = "top",
                                             trigger = "hover"),
                                   tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ratio")),
                                             tags$td(width = "60%", textInput(inputId = "comparetree_minratio", value = 0.05, label = NULL)))
                                   ), 
                                   bsTooltip(id = "comparetree_minratio",
                                             title = "Double (Default:1/20). If min.ratio is not NULL,all edge length of a phylogenetic tree should be greater than min.ratio*the longest edge length.If not, the edge length will be reset as min.ratio*longest edge length.",
                                             placement = "top",
                                             trigger = "hover"),
                                   br(),
                                   textInput(inputId = "comparetree_commoncol",
                                             label = div(style = "font-size:1.5em; font-weight:600; ", 'Common color'),
                                             value = "red"),
                                   bsTooltip(id = "comparetree_commoncol",
                                             title = "Color of common branches.",
                                             placement = "right",
                                             trigger = "hover"),
                                   
                                   sliderInput(inputId='comparetree_width',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image width'),min = 700,max = 1400, value = 1100, width = 500),
                                   sliderInput(inputId='comparetree_height',label = div(style = "font-size:1.5em; font-weight:600; ", 'Image height'),min = 400,max = 1000, value = 700, width = 500),
                                   # 
                                   fluidRow(
                                       column(
                                           width = 9,
                                           div(
                                               tags$button(
                                                   id = "submit_comparetree", type = "button", class = "action-button bttn",
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
                                 uiOutput("treemutsig_patientid_ui"),
                                 selectInput("treemutsig_getphylotree_method", label = div(style = "font-size:1.5em; font-weight:600;  ", "Tree construct method"),
                                             choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                                             selected = "NJ"),
                                 bsTooltip(id = "treemutsig_getphylotree_method",
                                           title = "Approach to construct phylogenetic trees",
                                           placement = "top",
                                           trigger = "hover"),
                                 tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                             tags$td(width = "60%", textInput(inputId = "treemutsig_getphylotree_minvaf", value = 0.02, label = NULL)))
                                 ), 
                                 bsTooltip(id = "treemutsig_getphylotree_minvaf",
                                           title = "The minimum value of vaf",
                                           placement = "top",
                                           trigger = "hover"),
                                 br(),
                                 tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                             tags$td(width = "60%", textInput(inputId = "treemutsig_getphylotree_minccf", value = 0, label = NULL)))
                                 ), 
                                 bsTooltip(id = "treemutsig_getphylotree_minccf",
                                           title = "The minimum value of ccf",
                                           placement = "top",
                                           trigger = "hover"),
                                 br(),
                                 tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Boostrap repetitions")),
                                             tags$td(width = "30%", textInput(inputId = "treemutsig_getphylotree_bootstraprepnum", value = 100, label = NULL)))
                                 ), 
                                 bsTooltip(id = "treemutsig_getphylotree_bootstraprepnum",
                                           title = "Bootstrap iterations. Default 100.",
                                           placement = "top",
                                           trigger = "hover"),
                                 br(),
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
                                 tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "85%", div(style = "font-size:1.5em; font-weight:600; ", "Minimal mutation number")),
                                           tags$td(width = "15%", textInput(inputId = "treemutsig_minmutcount", value = 15, label = NULL)))
                                 ), 
                                 bsTooltip(id = "treemutsig_minmutcount",
                                           title = 'The threshold for the variants in a branch. Default 15.',
                                           placement = "top",
                                           trigger = "hover"),
                                 br(),
                                 tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Signature cutoff")),
                                           tags$td(width = "30%", textInput(inputId = "treemutsig_signaturecutoff", value = 0.1, label = NULL)))
                                 ), 
                                 bsTooltip(id = "treemutsig_signaturecutoff",
                                           title = 'Discard any signature contributions with a weight less than this amount.Default: 0.1.',
                                           placement = "top",
                                           trigger = "hover"),
                                 br(),
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
                                 checkboxInput('treemutsig_usetumorsamplelabel', label = div(style = "font-size:1.5em; font-weight:600; padding-left:15px", 'Use tumor sample label'),value = FALSE, width = 400),
                                 bsTooltip(id = "treemutsig_usetumorsamplelabel",
                                           title = "Logical (Default: FALSE). Rename the 'Tumor_Sample_Barcode' with 'Tumor_Label'.",
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
                                 uiOutput("muttrunkbranch_patientid_ui"),
                                 selectInput("muttrunkbranch_getphylotree_method", label = div(style = "font-size:1.5em; font-weight:600;  ", "Tree construct method"),
                                             choices = c("NJ","MP","ML","FASTME.ols","FASTME.bal"),
                                             selected = "NJ"),
                                 bsTooltip(id = "muttrunkbranch_getphylotree_method",
                                           title = "Approach to construct phylogenetic trees",
                                           placement = "top",
                                           trigger = "hover"),
                                 tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min vaf")),
                                             tags$td(width = "60%", textInput(inputId = "muttrunkbranch_getphylotree_minvaf", value = 0.02, label = NULL)))
                                 ), 
                                 bsTooltip(id = "muttrunkbranch_getphylotree_minvaf",
                                           title = "The minimum value of vaf",
                                           placement = "top",
                                           trigger = "hover"),
                                 br(),
                                 tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "Min ccf")),
                                             tags$td(width = "60%", textInput(inputId = "muttrunkbranch_getphylotree_minccf", value = 0, label = NULL)))
                                 ), 
                                 bsTooltip(id = "muttrunkbranch_getphylotree_minccf",
                                           title = "The minimum value of ccf",
                                           placement = "top",
                                           trigger = "hover"),
                                 br(),
                                 tags$table(
                                     tags$tr(id = "inline", 
                                             width = "100%",
                                             tags$td(width = "50%", div(style = "font-size:1.5em; font-weight:600; ", "Boostrap repetitions")),
                                             tags$td(width = "30%", textInput(inputId = "muttrunkbranch_getphylotree_bootstraprepnum", value = 100, label = NULL)))
                                 ), 
                                 bsTooltip(id = "muttrunkbranch_getphylotree_bootstraprepnum",
                                           title = "Bootstrap iterations. Default 100.",
                                           placement = "top",
                                           trigger = "hover"),
                                 checkboxInput('muttrunkbranch_ct',div(style = "font-size:1.5em; font-weight:600; padding-left:15px ", 'CT'),value = FALSE),
                                 bsTooltip(id = "muttrunkbranch_ct",
                                           title = 'Distinction between C>T at CpG and C>T at other sites, Default FALSE',
                                           placement = "top",
                                           trigger = "hover"),
                                 tags$table(
                                   tags$tr(id = "inline", 
                                           width = "100%",
                                           tags$td(width = "30%", div(style = "font-size:1.5em; font-weight:600; ", "P-value")),
                                           tags$td(width = "60%", textInput(inputId = "muttrunkbranch_pvalue", value = 0.05, label = NULL)))
                                 ), 
                                 br(),
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
                               div(strong("PhyloTree-based analysis"),style = "font-size:27px; font-weight:500;"),
                               p("Systematic understanding of evolutionary relationship among regions plays a fundamental role in MRS study, where phylogenetic tree the a primary tool for describing these associations and interpreting ITH. MesKit is capable of constructing and comparing phylogenetic trees based on different methods, visualizing the rooted phylogenetic trees with annotation, as well as charactering mutational patterns based on phylogenetic trees.",
                                 style = "font-size:20px; font-weight:500;line-height:40px;"),
                               tabBox(
                                 id = 'sgt',
                                 side = 'left',
                                 selected = 'S_plotphylotree',
                                 width = "100%",
                                 height = "100%",
                                 tabPanel(
                                     title = div(icon("tree"), "Plot phylotree", style = "font-size:1.5em; font-weight:600; "),
                                     value = 'S_plotphylotree',
                                     uiOutput("phylotree.patientlist"),
                                     div(plotOutput("phylotree_plot",height = "100%",width = "100%"),align = "left"),
                                     br(),
                                     uiOutput("phylotree_downloadbutton_ui")
                                 ),
                                 tabPanel(
                                     title = div(icon("tree"), "Compare tree", style = "font-size:1.5em; font-weight:600; "), 
                                     value = 'S_comparetree',
                                     verbatimTextOutput("comparetree_dist"),
                                     br(),
                                     div(plotOutput('comparetree_plot', height = "100%", width = "100%"),align = "left"),
                                     uiOutput('comparetree_db_ui')
                                     # br(),
                                     # uiOutput('treemutsig_table_ui')
                                 ),
                                 tabPanel(
                                     title = div(icon("image"), "Mutational signature", style = "font-size:1.5em; font-weight:600; "), 
                                     value = 'S_treemutsig',
                                     # uiOutput('warningMessage_treemutsig'),
                                     uiOutput("treemutsig.patientlist"),
                                     uiOutput("treemutsig.samplelist"),
                                     div(plotOutput('treemutsig_plot', height = "100%", width = "100%"),align = "left"),
                                     uiOutput("treemutsig_download_button_ui"),
                                     # br(),
                                     # uiOutput('treemutsig_table_ui')
                                 ),
                                 tabPanel(
                                     title = div(icon("image"), "MutTrunkBranch", style = "font-size:1.5em; font-weight:600; "),
                                     value = 'S_muttrunkbranch',
                                     uiOutput("muttrunkbranch.patientlist"),
                                     br(),
                                     uiOutput('muttrunkbranch_table_ui'),
                                     div(plotOutput('muttrunkbranch_plot', height = "100%", width = "100%"),align = "left"),
                                     uiOutput("muttrunkbranch_download_button_ui")
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
                                tags$li(p("Jianyu Chen, chenjy327@mail2.sysu.edu.cn")),
                                tags$li(p("Xin Wang, wangx555@mail2.sysu.edu.cn")),
                                tags$li(p("Chengwei Wang, wangchw8@outlook.com"))
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
        $("header").find("nav").append(\'<span class="textnvbar"> MesKit: a tool kit for dissecting cancer evolution from multi-region derived tumor biopsies via somatic alterations</span>\');
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
                         .main-sidebar { font-size: 20px; }


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
