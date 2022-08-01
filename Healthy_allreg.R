library(shiny)
library(shinydashboard)
library(ggplot2)
library(reshape2)
library(DT)

source("global.R")
source("Fun_Rank_disease.R")
source("enrichr.R")

HealthyAllUI <- function(id,label="RankAll"){
  fluidPage(
    ###  First Row ###
    
    
    fluidRow(
      h3("Rank aggregation of Sex-biased genes of all brain regions"),
      h4("Please wait a few minutes for loading data"),
      column(width = 6,
             box(width = NULL, height = 300,imageOutput(NS(id, "image"))),
             box(width = NULL,
                 title = "Please choose cutoff for RRA rank aggregation",
                 #DTOutput(NS(id,"cutoff"))
             selectInput(NS(id,"DE_p_value"), "DE p-value cutoff:",
                         choices = list("0.05" = 0.05, "0.1" = 0.1
                                        ), selected = 0.05
                      ),
             selectInput(NS(id,"FC"),"FC cutoff", 
                      choices = list("1.2" = 1.2, "1.5" = 1.5, "2.0" = 2.0
                                     ), selected = 1.2
                      ),
             selectInput(NS(id,"RRA_p_value"),"RRA p-value cutoff", 
                         choices = list("0.05" = 0.05, "0.1" = 0.1
                                        ), selected = 0.05
                      ),
             
             )
      ),
      column(width = 6,
             box(width = NULL,
                 title = "The number of RRA sex-biased genes across brain regions",
                 #plotOutput(NS(id,"plt_RRAgene_num"))
                 #dataTableOutput(NS(id,"table"))
                 
             ),
             box(width = NULL,
                 title = "Correlation heatmap of sex-biased genes",
                 h4("Female"),
                 #plotOutput(NS(id,"F_corr_plot")),
                 h4("Male"),
                 #plotOutput(NS(id,"M_corr_plot"))
                 )
      )
      
      
    ),
    
    ### Second Row ###
    fluidRow(
      h3("Enrichment analysis of sex-biased genes"),
      tabBox(width = 12,
        title = "GO Enrichment",
             id = "GO",
             tabPanel("GO Biological Process",
                      h4("Female"),
                      #plotOutput(NS(id,"GO_BP_F")),
                      h4("Male"),
                      #plotOutput(NS(id,"GO_BP_M"))
             ),
             tabPanel("GO Molecular Function",
                      h4("Female"),
                      #plotOutput(NS(id,"GO_MO_F")),
                      h4("Male"),
                      #plotOutput(NS(id,"GO_MO_M"))
             ),
             tabPanel("GO Cellular Component",
                      h4("Female"),
                      #plotOutput(NS(id,"GO_CC_F")),
                      h4("Male"),
                      #plotOutput(NS(id,"GO_CC_M")),
             )
      )
      
     
    ),
    
    ## Third Row ###
    fluidRow(
      #h3("DisGenet"),
      box(width = 12,
          title = "DisGeNET Enrichment",
          h4("Female"),
          #plotOutput(NS(id,"Dis_F")),
          h4("Male"),
          #plotOutput(NS(id,"Dis_M"))
          )
      
    ),
    
    ## Forth Row ##
    fluidRow(
      box(width = 12,
          title = "KEGG Enrichment",
          h4("Female"),
          #plotOutput(NS(id,"KEGG_F")),
          h4("Male"),
          #plotOutput(NS(id,"KEGG_M"))
          )
      
      
    ),
    
    
  )
  
}

#SummaryServer <- function(id,data) {
#  stopifnot(is.reactive(data))
HealthyAllServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
    # Dataset summary plot image
    output$image <- renderImage({
      filename <- normalizePath(file.path(paste0('www/Healthy1.png')))
      list(
        src = filename, 
        height = 250
      )
    }, deleteFile = FALSE)
    
    # 1. Get input cutoff from user for RRA
    print("Get input cutoff")
    
    
    
    
    
    
  })
}
