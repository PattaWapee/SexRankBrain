library(shiny)
library(shinydashboard)

source('Healthy_reg.R')
source('Healthy_allreg2.R')
source('Summary.R')

WebApp <- function() {
  
  ########## UI ##############
  ui <- dashboardPage(
    
    dashboardHeader(title = "SexRankBrain"),
    
    dashboardSidebar(
      sidebarMenu(
        menuItem("Healthy (specific brain region)", tabName = "Healthy1"),
        menuItem("Healthy (across brain region)", tabName = "Healthy2"),
        menuItem("About", tabName = "About")
      )
    ),
    
    
    ### body ##
    dashboardBody(
      tabItems(
        tabItem(tabName = "Healthy1",
                h2("Sex-biased genes rank for healthy brain samples"),
                Healthy_UI('Healthy1')
        ),
        
        
        tabItem(tabName = "Healthy2",
                h2("Sex-biased genes rank across brain regions for healthy brain samples"),
                HealthyAllUI('Healthy2')
        ),
        
        
        tabItem(tabName = "About",
                SummaryUI('Summary')
                #h2("About"),
                
        )
      )
      
      
    )
  )
  
  
  
  ########## SERVER #######
  server <- function(input, output, session) {
    
    #RRA_data <- RankServer('Rank1')
    Healthy_Server('Healthy1')
    HealthyAllServer('Healthy2')
    
    
  }
  shinyApp(ui, server)  
}


WebApp()
