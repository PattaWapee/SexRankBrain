library(shiny)
library(shinydashboard)
library(ggplot2)
library(reshape2)

source("Fun_Rank_disease.R")
source("global.R")
source("enrichr.R")


  
Healthy_UI <- function(id,label="Rank"){
  fluidPage(
    ## First Row: Filtering sex-biased genes from differential expression topTable.
    
    fluidRow(
      h3("Step1: Filtering sex-biased genes of each dataset before rank aggregation by p-value & logFC cutoff"),
      
      ##### INPUT BOX #######
      box(title = "Sex-biased gene differential gene rank",
          selectInput(NS(id,"region"),"brain region",
                      choices = list("AMY" = 1, "CBC" = 2, "CC" = 3,
                                     "FC" = 4, "HIP" = 5, "MED" = 6,
                                     "OC" = 7, "PL" = 8, "STR" = 9,
                                     "TC" = 10, "THA" = 11
                      ), selected = 1
          ),
          sliderInput(NS(id,"DE_p_value"), "p-value cutoff:", min=0, max = 1, value = 0.05),
          numericInput(NS(id,"FC"),"FC cutoff", value = 1.2)
          
          
      ),
      
      
      ######## OUTPUT BOX ######
      
      box(title='The number of sex-biased genes',
          plotOutput(NS(id,"plt_gene_num"))
          )
      
      
      
    ),
    ## Second Row: Filtering genes from RRA
    fluidRow(
      h3("Step2: Combind all dataset ranks by RRA and filtered sex-biased genes by RRA p-value"),
      
      box(title = 'Robust Rank aggregation cutoff',
          sliderInput(NS(id,"RRA_p_value"), "RRA p-value cutoff:", min=0, max = 1, value = 1)
          ),
      #box(title = 'Rank aggregation table of male-biased genes',
      #    dataTableOutput(NS(id,"RRAM")))
      tabBox(title = "RRA output",
             id= "tabRRA",
             tabPanel("Output","The number of RRA genes",
                      plotOutput(NS(id,"plt_RRAnum"))
                      #dataTableOutput(NS(id, "RRAF"))
             ),
             
             tabPanel("Female RRA table","Rank aggregation table of female-biased genes",        
                      dataTableOutput(NS(id, "RRAFtable"))
                      ),
             tabPanel("Male RRA table","Rank aggregation table of male-biased genes",
                      dataTableOutput(NS(id, "RRAMtable"))
                      #dataTableOutput(NS(id, "RRA"))
                      )
             )
    ),
    ## Third Row: GO Enrichment of Filtered sex-biased genes of RRA
    fluidRow(
      h3("Step3: Enrichment analysis of sex-biased genes of RRA"),
      
      tabBox(title = "GO Enrichment",
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
          ),
      
      box(title = "DisGeNET Enrichement",
          h4("Female"),
          plotOutput(NS(id,"Dis_F")),
          h4("Male"),
          plotOutput(NS(id,"Dis_M"))
          ),
      box(title = "KEGG enrichment",
          h4("Female"),
          plotOutput(NS(id,"KEGG_F")),
          h4("Male"),
          plotOutput(NS(id,"KEGG_M"))
          )
      
    
    )
    
  )
    
    
  
    
}


Healthy_Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    

    
    ########### MAIN RUN FOR EACH REGION #############
    # 1. filtered genes for each dataset
    
    sig_reg_F <- reactive( Filter_gene_reg_dis(input$region,input$DE_p_value,input$FC,'F','Healthy') )
    sig_reg_M <- reactive( Filter_gene_reg_dis(input$region,input$DE_p_value,input$FC,'M','Healthy') )
    
    # 2. Plot number of filtered genes for each dataset
    
    output$plt_gene_num <- renderPlot({
      num_fil_df <- numfil_gene(sig_reg_F() ,sig_reg_M() )
      
      ggplot(num_fil_df, aes(x=dataset, y=value, fill=gender)) +
        geom_bar(stat="identity",position='dodge') +
        scale_fill_manual(values=c("deeppink", "dodgerblue"))+
        theme(axis.text.x = element_text(angle = 45))      
    })
    
    # 3. Rank aggregation of filtered genes with RRA pvalue cutoff
    
    RRA_F <- reactive({
      Fsig <- lapply(sig_reg_F() , rownames)
      RRA_F <- RA.analysis(Fsig)
      sig_RRA_F <- RRA_F[ RRA_F['P.value'] < input$RRA_p_value,]
    })
    
    RRA_M <- reactive({
      Msig <- lapply(sig_reg_M() , rownames)
      RRA_M <- RA.analysis(Msig)
      sig_RRA_M <- RRA_M[ RRA_M['P.value'] < input$RRA_p_value,]
    })
    output$RRAFtable <- renderDataTable(RRA_F())
    output$RRAMtable <- renderDataTable(RRA_M())
    
    # 4. Plot number of filtered genes for each dataset and RRA
    output$plt_RRAnum <- renderPlot({
      num_RRA <- data.frame(Female=c(nrow(RRA_F() )),Male=c(nrow(RRA_M() )) )
      barplot(t(as.matrix(num_RRA)),beside=TRUE,col = c('deeppink','dodgerblue'))
    })
    
    # 5. Enrichment of GO, DisGeNET
    
    output$GO_BP_F <- renderPlot(plt_enriched_reg(RRA_F(), "GO_Biological_Process_2021"))
    output$GO_BP_M <- renderPlot(plt_enriched_reg(RRA_M(), "GO_Biological_Process_2021"))
    
    output$GO_MO_F <- renderPlot(plt_enriched_reg(RRA_F(), "GO_Molecular_Function_2021"))
    output$GO_MO_M <- renderPlot(plt_enriched_reg(RRA_M(), "GO_Molecular_Function_2021"))
    
    output$GO_CC_F <- renderPlot(plt_enriched_reg(RRA_F(), "GO_Cellular_Component_2021"))
    output$GO_CC_M <- renderPlot(plt_enriched_reg(RRA_M(), "GO_Cellular_Component_2021"))
    
    output$Dis_F <- renderPlot({
      plt_enriched_reg(RRA_F(), "DisGeNET")
    })
    output$Dis_M <- renderPlot({
      plt_enriched_reg(RRA_M(), "DisGeNET")
    })
    
    # KEGG
    output$KEGG_F <- renderPlot({
      plt_enriched_reg(RRA_F(), "KEGG_2021_Human")
    })
    output$KEGG_M <- renderPlot({
      plt_enriched_reg(RRA_M(), "KEGG_2021_Human")
    })
    
    
  
  })
}
