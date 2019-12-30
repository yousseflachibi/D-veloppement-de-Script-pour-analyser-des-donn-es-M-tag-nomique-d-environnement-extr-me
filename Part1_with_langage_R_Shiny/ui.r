## add sub menu items in shinydashboard sidebar

# load the required packages
library(shiny)

library(shinyWidgets)
#install.packages("shinydashboard")
library(shinydashboard)

shinyUI(
  dashboardPage(
    dashboardHeader(title = "Demo", titleWidth = 200),
    dashboardSidebar(
      sidebarMenu(id = 'sidebarmenu',
                  # first menu item
                  menuItem("Menu", tabName = "Dashboard", icon = icon("dashboard")),
                  
                  # second menu item with 2 sub menus
                  menuItem('Utils',
                           icon = icon('line-chart'),
                           menuSubItem('Search',
                                       tabName = 'Search_tap',
                                       icon = icon('line-chart')),
                           menuSubItem('Proteins Freq',
                                       tabName = 'Protein_Tap',
                                       icon = icon('line-chart')
                           )))),
    
    
    dashboardBody(
      tabItems(
        tabItem("Dashboard",h4("this is the Dashboard tab page")),
        tabItem("Search_tap",
                headerPanel("Search ..."),
                sidebarPanel(#shinythemes::themeSelector(),
          shinyjs::useShinyjs(),
          textInput(inputId = "input_1",
                    label = "Organisme : ",
                    value = "",
                    placeholder = "Organ..."),
          textInput(inputId = "input_2",
                    label = "Parametre : ",
                    value = "",
                    placeholder = "Param..."),
          fluidRow(
            column(6, align="center", offset = 3,
                   actionButton("go_1",label = "Search"),
                   tags$style(type='text/css'))),
          
          fluidRow(
            column(6, align="center", offset = 3,
                   shinyjs::hidden(p(id = "text_1", "Processing..."),
                                   tags$style(type='text/css')))
          )
        ),
        # Main panel for displaying outputs ----
        mainPanel()
        ),
        tabItem("Protein_Tap", tags$style(type = "text/css", " {height: 1000px important;}"),
                headerPanel("Protein ...."),
                sidebarPanel(#shinythemes::themeSelector(),
                  shinyjs::useShinyjs(),
                  fileInput(inputId = "input_3",
                            label = "Choise the file fna: ",
                            accept = c('.fna'),multiple = T),
                  hr(),
                  tags$head(
                    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
                  ),
                  fluidRow(
                    column(width = 12,
                           h5("Supprersion Proteins no commun :"),
                           progressBar(id = "pb1", value = 0,display_pct = TRUE)
                           ,hr())
                    ,
                    column(width = 12,
                           h5("Creer csv Proteins communs :"),
                           progressBar(id = "pb2", value = 0,display_pct = TRUE)
                           ,hr())
                  ),
                  fluidRow(
                    column(6, align="center", offset = 3,
                           actionButton("go_2",label = "Submit"),
                           tags$style(type='text/css'))),
                  fluidRow(
                    column(6, align="center", offset = 3,
                           shinyjs::hidden(p(id = "text_2", "Processing..."),
                                           tags$style(type='text/css')))
                  )
                ),
                
                # Main panel for displaying outputs ----
                mainPanel(
                  tabsetPanel(
                    # using iframe along with tags() within tab to display pdf with scroll, height and width could be adjusted
                    tabPanel("1..."
                             ,plotOutput(outputId = "distPlot")
                             ))
                  
                  
                  )
                )
      )
    )
  )
)
