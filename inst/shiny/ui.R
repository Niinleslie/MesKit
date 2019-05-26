options(spinner.type=4)

#required packages
suppressMessages(library(shiny))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyWidgets))
suppressMessages(library(shinycssloaders))

#sider bar----

sidebar <- dashboardSidebar(
  width = 300,
  sidebarMenu(id="sidername",selected='home',
    menuItem("Home", tabName = "home", icon = icon("home")),
    menuItem("MesKit overview", tabName = "MesKit", icon = icon("eye")),
    menuItem("ITH evaluation", tabName = "ITH", icon = icon("location-arrow")),
    menuItem("Clonal analysis", tabName = "clone", icon = icon("th",lib = "glyphicon")),
    menuItem("Functional analysis", tabName = "function", icon = icon("bar-chart")),
    menuItem("Other analysis", tabName = "Survival", icon = icon("line-chart"))
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

#Main function----
shinyUI(
  dashboardPage(skin = "blue",
    dashboardHeader(title = "Meskit: Analysis and visualize multi-sample whole-exome sequencing data",
                  titleWidth = 650),
    sidebar,
    dashboardBody(
      tags$head(
        tags$style(HTML(".shiny-output-error-validation {color: brown;}")),
        tags$link(rel = "stylesheet", type = "text/css", href = "css/main.css")
      ),
      tabItems(
        bodyHome
      )
    )
  )
)
