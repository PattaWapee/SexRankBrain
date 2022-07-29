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
      column(width = 6,
             box(width = NULL, height = 300,imageOutput(NS(id, "image"))),
             box(width = NULL,
                 title = "Editable cutoff for RRA rank aggregation",
                 DTOutput(NS(id,"cutoff"))
             )
      ),
      column(width = 6,
             box(width = NULL,
                 title = "The number of RRA sex-biased genes across brain regions",
                 plotOutput(NS(id,"plt_RRAgene_num"))
                 #dataTableOutput(NS(id,"table"))
                 
             ),
             box(width = NULL,
                 title = "Correlation heatmap of sex-biased genes",
                 h4("Female"),
                 plotOutput(NS(id,"F_corr_plot")),
                 h4("Male"),
                 plotOutput(NS(id,"M_corr_plot"))
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
                      plotOutput(NS(id,"GO_BP_F")),
                      h4("Male"),
                      plotOutput(NS(id,"GO_BP_M"))
             ),
             tabPanel("GO Molecular Function",
                      h4("Female"),
                      plotOutput(NS(id,"GO_MO_F")),
                      h4("Male"),
                      plotOutput(NS(id,"GO_MO_M"))
             ),
             tabPanel("GO Cellular Component",
                      h4("Female"),
                      plotOutput(NS(id,"GO_CC_F")),
                      h4("Male"),
                      plotOutput(NS(id,"GO_CC_M")),
             )
      )
      
     
    ),
    
    ## Third Row ###
    fluidRow(
      #h3("DisGenet"),
      box(width = 12,
          title = "DisGeNET Enrichment",
          h4("Female"),
          plotOutput(NS(id,"Dis_F")),
          h4("Male"),
          plotOutput(NS(id,"Dis_M"))
          )
      
    ),
    
    ## Forth Row ##
    fluidRow(
      box(width = 12,
          title = "KEGG Enrichment",
          h4("Female"),
          plotOutput(NS(id,"KEGG_F")),
          h4("Male"),
          plotOutput(NS(id,"KEGG_M"))
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
    
    cutoff <- reactiveValues(df=inputdata)
    output$cutoff <- renderDT(cutoff$df,selection = 'none', server = TRUE, editable ="cell")
    observeEvent(input$cutoff_cell_edit, {
      
      cutoff$df <- editData(cutoff$df, input$cutoff_cell_edit, 'cutoff')
      
    })
    
    
    # 2. Run RRA and filter from pvalue
    print("Run RRA and filter from pvalue")
    RRA_F <-reactive( {cal_filter_RRAreg(cutoff$df,"F","Healthy")} )
    RRA_M <-reactive( {cal_filter_RRAreg(cutoff$df,"M","Healthy")} )
    print("Calculate number of sex-biased genes")
    num_fil_df <- reactive( numfil_gene(RRA_F() ,RRA_M() ) )
    print("Plot number of sex-biased genes")
    plt_num_fil <- reactive({
      ggplot(num_fil_df(), aes(x=dataset, y=value, fill=gender)) +
        geom_bar(stat="identity",position='dodge') +
        scale_fill_manual(values=c("deeppink", "dodgerblue"))+
        theme(axis.text.x = element_text(angle = 45))    
    })
    
    
    # 3. Plot number of filtered genes for each dataset
    output$plt_RRAgene_num <- renderPlot(plt_num_fil())
    #output$table <- renderDataTable( num_fil_df()  )
    
    # 4. Plot correlation matrix heatmap
    library(RColorBrewer)
    m_coul <- colorRampPalette(brewer.pal(9, "YlGnBu"))(25)
    f_coul <- colorRampPalette(brewer.pal(9, "RdPu"))(25)
    output$F_corr_plot <-  renderPlot({
      F_corr <- correlation_reg( RRA_F() )
      heatmap(F_corr,col=f_coul,scale = "row")
    })
    
    output$M_corr_plot <-  renderPlot({
      M_corr <- correlation_reg( RRA_M() )
      heatmap(M_corr,col=m_coul,scale = "row")
    })
    
    # 4. Enrichment of GO, DisGeNET
    output$GO_BP_F <- renderPlot({
      run_plt_enrich(RRA_F(), "GO_Biological_Process_2021", "F")
    })
    
    output$GO_BP_M <- renderPlot({
      run_plt_enrich(RRA_M(), "GO_Biological_Process_2021", "M")
    })
    
    output$GO_MO_F <- renderPlot({
      run_plt_enrich(RRA_F(), "GO_Molecular_Function_2021", "F")
    })
    
    output$GO_MO_M <- renderPlot({
      run_plt_enrich(RRA_M(), "GO_Molecular_Function_2021", "M")
    })
    
    output$GO_CC_F <- renderPlot({
      run_plt_enrich(RRA_F(), "GO_Cellular_Component_2021", "F")
    })
    
    output$GO_CC_M <- renderPlot({
      run_plt_enrich(RRA_M(), "GO_Cellular_Component_2021", "M")
    })
    
    # DisGeNET
    output$Dis_F <- renderPlot({
      run_plt_enrich(RRA_F(), "DisGeNET", "F")
    })
    
    
    output$Dis_M <- renderPlot({
      run_plt_enrich(RRA_M(), "DisGeNET", "M")
    })
    
    # KEGG
    output$KEGG_F <- renderPlot({
      run_plt_enrich(RRA_F(), "KEGG_2021_Human", "F")
    })
    
    
    output$KEGG_M <- renderPlot({
      run_plt_enrich(RRA_M(), "KEGG_2021_Human", "M")
    })
    
  })
}
