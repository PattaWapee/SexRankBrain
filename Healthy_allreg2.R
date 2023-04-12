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
             selectInput(NS(id,"Cutoff"), "Please select a cutoff:",
                         choices = list(
                                        "DE FC cutoff:1.0, RRA p-adjust value:0.05" = 4,
                                        "DE p-value:0.05, DE FC:1.2, RRA p-value:0.05" = 1, 
                                        "DE p-value:0.05, DE FC:1.5, RRA p-value:0.05" = 2,
                                        "DE p-value:0.05, DE FC:2.0, RRA p-value:0.05" = 3
                                        ), selected = 4),
             
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
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_num_img/num_allreg_DE_0.05_FC_1.2_RRA_0.05.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_num_img/num_allreg_DE_0.05_FC_1.5_RRA_0.05.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_num_img/num_allreg_DE_0.05_FC_2_RRA_0.05.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_num_img/num_allreg_RRA_adj_0.05.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    
    
    # Corelation heatmap
    
    output$corrF_img <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_DE_0.05_FC_1.2_RRA_0.05_F.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_DE_0.05_FC_1.5_RRA_0.05_F.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_DE_0.05_FC_2_RRA_0.05_F.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_RRA_adj_0.05_F.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    output$corrM_img <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_DE_0.05_FC_1.2_RRA_0.05_M.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_DE_0.05_FC_1.5_RRA_0.05_M.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_DE_0.05_FC_2_RRA_0.05_M.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_corr_img/corr_RRA_adj_0.05_M.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    
    # GO enrichment
    output$GO_BP_F <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_BP_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_BP_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_BP_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_BP_female.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    output$GO_BP_M <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_BP_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_BP_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_BP_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_BP_male.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    # GO MO
    output$GO_MO_F <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_MF_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_MF_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_MF_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_MF_female.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    output$GO_MO_M <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_MF_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_MF_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_MF_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_MF_male.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    # GO CC
    output$GO_CC_F <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_CC_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_CC_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_CC_female.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_CC_female.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    output$GO_CC_M <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_CC_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_CC_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_CC_male.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_CC_male.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)

    # DisGenet
    output$Dis_F <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_female_DisGeNET_curated.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_female_DisGeNET_curated.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_female_DisGeNET_curated.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_Dis_F.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
   
    output$Dis_M <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_male_DisGeNET_curated.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_male_DisGeNET_curated.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_male_DisGeNET_curated.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_FC_1_RRA_adj_0.05_Dis_M.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    
    # KEGG
    output$KEGG_F <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_female_KEGG.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_female_KEGG.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_female_KEGG.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/No_enriched.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    

    output$KEGG_M <- renderImage({
      if (input$Cutoff == 1) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.2_RRA_0.05_male_KEGG.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 2) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_1.5_RRA_0.05_male_KEGG.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 3) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/DE_0.05_FC_2_RRA_0.05_male_KEGG.png'))
        list(src = filename, height = 300, width = 400)
      } else if (input$Cutoff == 4) {
        filename <- normalizePath(file.path('data/pre_result_allreg/plt_enrich_img/No_enriched.png'))
        list(src = filename, height = 300, width = 400)
      } else {
        list()
      }
    }, deleteFile = FALSE)
    



  })
}
