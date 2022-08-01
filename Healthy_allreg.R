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
                 imageOutput(NS(id, "num_img"))
                 
             ),
             box(width = NULL,
                 title = "Correlation heatmap of sex-biased genes",
                 h4("Female"),
                 imageOutput(NS(id, "corrF_img")),
                 h4("Male"),
                 imageOutput(NS(id, "corrM_img")),
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
                      imageOutput(NS(id,"GO_BP_F")),
                      h4("Male"),
                      imageOutput(NS(id,"GO_BP_M"))
             ),
             tabPanel("GO Molecular Function",
                      h4("Female"),
                      imageOutput(NS(id,"GO_MO_F")),
                      h4("Male"),
                      imageOutput(NS(id,"GO_MO_M"))
             ),
             tabPanel("GO Cellular Component",
                      h4("Female"),
                      imageOutput(NS(id,"GO_CC_F")),
                      h4("Male"),
                      imageOutput(NS(id,"GO_CC_M"))
             )
      )
      
     
    ),
    
    ## Third Row ###
    fluidRow(
      h3("DisGenet"),
      box(width = 12,
          title = "DisGeNET Enrichment",
          h4("Female"),
          imageOutput(NS(id,"Dis_F")),
          h4("Male"),
          imageOutput(NS(id,"Dis_M"))
          )
      
    ),
    
    ## Forth Row ##
    fluidRow(
      box(width = 12,
          title = "KEGG Enrichment",
          h4("Female"),
          imageOutput(NS(id,"KEGG_F")),
          h4("Male"),
          imageOutput(NS(id,"KEGG_M"))
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
    
    # Image number of sex-biased genes
    output$num_img <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_num_img/num_allreg_DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'.png')))
      
      list(src = filename, height = 300, width = 400)
    }, deleteFile = FALSE)
    
    # Corelation heatmap
    output$corrF_img <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_corr_img/corr_DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_F.png')))
      
      list(src = filename, height = 300, width = 300)
    }, deleteFile = FALSE)

    output$corrM_img <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_corr_img/corr_DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_M.png')))
      
      list(src = filename, height = 300, width = 300)
    }, deleteFile = FALSE)
    
    # GO enrichment
    output$GO_BP_F <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_BP_female.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
    
    output$GO_BP_M <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_BP_male.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
    
    # GO MO
    output$GO_MO_F <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_MO_female.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
    
    output$GO_MO_M <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_MO_male.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
    
    # GO CC
    output$GO_CC_F <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_CC_female.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
    
    output$GO_CC_M <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_CC_male.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)

    # DisGenet
    output$Dis_F <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_female_DisGeNET_curated.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
   
    output$Dis_M <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_male_DisGeNET_curated.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
    
    # KEGG
    output$KEGG_F <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_female_KEGG.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)
    
    output$KEGG_M <- renderImage({
      filename <- normalizePath(file.path(paste0('data/pre_result_allreg/plt_enrich_img/DE_',
                                                 as.character(input$DE_p_value),'_FC_',
                                                 as.character(input$FC),'_RRA_',as.character(input$RRA_p_value),'_male_KEGG.png')))
      
      list(src = filename, height = 400, width = 400)
    }, deleteFile = FALSE)




  })
}
