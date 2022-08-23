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
      h4("Please wait for a few minutes for calculation"),
      
      ##### INPUT BOX #######
      box(title = "Sex-biasedMne differential gene rank",
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
          plotOutput(NS(id,"plt_gene_num")),
          downloadButton(NS(id,'downloadNum'), 'Download plot' )
          )
      
      
      
    ),
    ## Second Row: Filtering genes from RRA
    fluidRow(
      h3("Step2: Combind all dataset ranks by RRA and filtered sex-biased genes by RRA p-value"),
      
      box(title = 'Robust Rank aggregation cutoff',
          sliderInput(NS(id,"RRA_p_value"), "RRA p-value cutoff:", min=0, max = 1, value = 0.05),
          
          ),
      #box(title = 'Rank aggregation table of male-biased genes',
      #    dataTableOutput(NS(id,"RRAM")))
      tabBox(title = "RRA output",
             id= "tabRRA",
             tabPanel("Output","The number of RRA genes",
                      plotOutput(NS(id,"plt_RRAnum")),
                      downloadButton(NS(id,'downloadRRANum'), 'Download plot' )
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
                      downloadButton(NS(id,'downloadGO_BP_F'), 'Download plot' ),
                      
                      h4("Male"),
                      plotOutput(NS(id,"GO_BP_M")),
                      downloadButton(NS(id,'downloadGO_BP_M'), 'Download plot' ),
                      ),
             tabPanel("GO Molecular Function",
                      h4("Female"),
                      plotOutput(NS(id,"GO_MO_F")),
                      downloadButton(NS(id,'downloadGO_MO_F'), 'Download plot' ),
                      h4("Male"),
                      plotOutput(NS(id,"GO_MO_M")),
                      downloadButton(NS(id,'downloadGO_MO_M'), 'Download plot' ),
                      ),
             tabPanel("GO Cellular Component",
                      h4("Female"),
                      plotOutput(NS(id,"GO_CC_F")),
                      downloadButton(NS(id,'downloadGO_CC_F'), 'Download plot' ),
                      h4("Male"),
                      plotOutput(NS(id,"GO_CC_M")),
                      downloadButton(NS(id,'downloadGO_CC_M'), 'Download plot' ),
                      )
          ),
      
      box(title = "DisGeNET Enrichement",
          h4("Female"),
          plotOutput(NS(id,"Dis_F")),
          downloadButton(NS(id,'downloadDis_F'), 'Download plot' ),
          h4("Male"),
          plotOutput(NS(id,"Dis_M")),
          downloadButton(NS(id,'downloadDis_M'), 'Download plot' ),
          ),
      box(title = "KEGG enrichment",
          h4("Female"),
          plotOutput(NS(id,"KEGG_F")),
          downloadButton(NS(id,'downloadKEGG_F'), 'Download plot' ),
          h4("Male"),
          plotOutput(NS(id,"KEGG_M")),
          downloadButton(NS(id,'downloadKEGG_M'), 'Download plot' ),
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

    # save plot
    num_df <- reactive( numfil_gene(sig_reg_F(), sig_reg_M() ) )
    num_plot <- reactive( ggplot(num_df(), aes(x=dataset, y=value, fill=gender)) +
                         geom_bar(stat="identity",position='dodge') +
                         scale_fill_manual(values=c("deeppink", "dodgerblue"))+
                         theme(axis.text.x = element_text(angle = 45))
                       )


    output$downloadNum <- downloadHandler(
    filename = "Sex_biased_genes.tiff",
    content = function(file) {
        tiff(file, res = 300, units = 'in', height =5, width = 5)
        plot(num_plot())
        dev.off()
    }
    )
    
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
    output$RRAFtable <- renderDataTable(RRA_F(),server = FALSE,extensions = 'Buttons',
                                        options = list(scrollX=TRUE, lengthMenu = c(5,10,15),
                                                       paging = TRUE, searching = TRUE,
                                                       fixedColumns = TRUE, autoWidth = TRUE,
                                                       ordering = TRUE, dom = 'tBp',
                                                       #ordering = TRUE, dom = 'Bfrtip',
                                                       buttons = c('copy', 'csv', 'excel','pdf'))
                                        )
    output$RRAMtable <- renderDataTable(RRA_M(),server = FALSE,extensions = 'Buttons',
                                        options = list(scrollX=TRUE, lengthMenu = c(5,10,15),
                                                       paging = TRUE, searching = TRUE,
                                                       fixedColumns = TRUE, autoWidth = TRUE,
                                                       ordering = TRUE, dom = 'tBp',
                                                       #ordering = TRUE, dom = 'Bfrtip',
                                                       buttons = c('copy', 'csv', 'excel','pdf')))
    
    # 4. Plot number of filtered genes for each dataset and RRA
    output$plt_RRAnum <- renderPlot({
      num_RRA <- data.frame(Female=c(nrow(RRA_F() )),Male=c(nrow(RRA_M() )) )
      barplot(t(as.matrix(num_RRA)),beside=TRUE,col = c('deeppink','dodgerblue'))
    })

    plt_RRA_num <- reactive({
      num_RRA <- data.frame(Number_of_genes=c(nrow(RRA_F() ), nrow(RRA_M() )),
                            Sex = c('Female', 'Male')
      )
      ggplot(num_RRA, aes(x= Sex, y = Number_of_genes, fill = Sex))+
        geom_bar(stat = 'identity')+
        scale_fill_manual(values=c("deeppink", "dodgerblue"))

    })

    # 5. User download RRA plot
    output$downloadRRANum <- downloadHandler(
    filename = "Number_of_sex_biased_RRA_filtered.tiff",
    content = function(file) {
        tiff(file, res = 300, units = 'in', height =5, width = 5)
        plot(plt_RRA_num())
        dev.off()
    }
    )



    # 6. Enrichment of GO, DisGeNET
    

    ############
    ## GO BP####
    ############
    output$GO_BP_F <- renderPlot(plt_enriched_reg(RRA_F(), "GO_Biological_Process_2021"))
    output$GO_BP_M <- renderPlot(plt_enriched_reg(RRA_M(), "GO_Biological_Process_2021"))
    pltGO_BP_F <- reactive( plt_enriched_reg(RRA_F(), "GO_Biological_Process_2021" ))
    pltGO_BP_M <- reactive( plt_enriched_reg(RRA_M(), "GO_Biological_Process_2021" ))
    # User download GO plot
    output$downloadGO_BP_F <- downloadHandler(
                                              filename = "GO_BP_F.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltGO_BP_F())
                                                dev.off()}
                              )
    output$downloadGO_BP_M <- downloadHandler(
                                              filename = "GO_BP_M.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltGO_BP_M())
                                                dev.off()}
                              )
    ###########
    ## GO MO ##
    ###########
    output$GO_MO_F <- renderPlot(plt_enriched_reg(RRA_F(), "GO_Molecular_Function_2021"))
    output$GO_MO_M <- renderPlot(plt_enriched_reg(RRA_M(), "GO_Molecular_Function_2021"))
    pltGO_MO_F <- reactive( plt_enriched_reg(RRA_F(), "GO_Molecular_Function_2021"))
    pltGO_MO_M <- reactive( plt_enriched_reg(RRA_M(), "GO_Molecular_Function_2021"))
    
    # User download GO plot
    output$downloadGO_MO_F <- downloadHandler(
                                              filename = "GO_MO_F.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltGO_MO_F())
                                                dev.off()}
                              )
    output$downloadGO_MO_M <- downloadHandler(
                                              filename = "GO_MO_M.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltGO_MO_M())
                                                dev.off()}
                              )
    output$GO_CC_F <- renderPlot(plt_enriched_reg(RRA_F(), "GO_Cellular_Component_2021"))
    output$GO_CC_M <- renderPlot(plt_enriched_reg(RRA_M(), "GO_Cellular_Component_2021"))
    pltGO_CC_F <- reactive( plt_enriched_reg(RRA_F(), "GO_Cellular_Component_2021"))
    pltGO_CC_M <- reactive( plt_enriched_reg(RRA_M(), "GO_Cellular_Component_2021"))

    output$downloadGO_CC_F <- downloadHandler(
                                              filename = "GO_CC_F.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltGO_CC_F())
                                                dev.off()}
                              )
    output$downloadGO_CC_M <- downloadHandler(
                                              filename = "GO_CC_M.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltGO_CC_M())
                                                dev.off()})
    ###############
    ## DisGeNET ###
    ##############
    output$Dis_F <- renderPlot({
      plt_enriched_reg(RRA_F(), "DisGeNET")
    })
    output$Dis_M <- renderPlot({
      plt_enriched_reg(RRA_M(), "DisGeNET")
    })
    pltDis_F <- reactive( plt_enriched_reg(RRA_F(), "DisGeNET"))
    pltDis_M <- reactive( plt_enriched_reg(RRA_M(), "DisGeNET"))
    
    output$downloadDis_F <- downloadHandler(
                                              filename = "DisGeNET_F.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltDis_F())
                                                dev.off()}
                              )
    output$downloadDis_M <- downloadHandler(
                                              filename = "DisGeNET_M.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltDis_M())
                                                dev.off()})
    ########
    # KEGG #
    ########
    output$KEGG_F <- renderPlot({
      plt_enriched_reg(RRA_F(), "KEGG_2021_Human")
    })
    output$KEGG_M <- renderPlot({
      plt_enriched_reg(RRA_M(), "KEGG_2021_Human")
    })
    pltKEGG_F <- reactive( plt_enriched_reg(RRA_F(), "KEGG_2021_Human"))
    pltKEGG_M <- reactive( plt_enriched_reg(RRA_M(), "KEGG_2021_Human"))
    output$downloadKEGG_F <- downloadHandler(
                                              filename = "KEGG_F.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltKEGG_F())
                                                dev.off()}
                              )
    output$downloadKEGG_M <- downloadHandler(
                                              filename = "KEGG_M.tiff",
                                              content = function(file) {
                                                tiff(file, res = 300, units = 'in', height =5, width = 5)
                                                plot(pltKEGG_M())
                                                dev.off()}
                              )
    
    
  
  })




}
