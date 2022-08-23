library(shiny)
library(shinydashboard)
library(ggplot2)
library(reshape2)


SummaryUI <- function(id, label = "summary"){
  fluidPage(
            # First Section
            fluidRow(
                     h2("About"),

                     h5("SexBrainRank application is build 
                        to help user explore the result data 
                        from the study named “Sex-Associated Differences 
                        in Human Brain”.Please check the abstract section below."),

                      h4("Abstract"),
                      h5("Background"),
                      h6("Sex dimorphism is highly prominent in mammals with
                         many systematic differences between male and female 
                         form of the species. Sex differences arise from genetic
                         and environmental factors. Accordingly, the fundamental 
                         social and cultural stratification factors for humans is 
                         sex. It distinguishes individuals most prominently on the 
                         reproductive traits, but also affects many of the other 
                         related traits. This can lead to different disease 
                         susceptibilities and treatment responses across sexes. 
                         Sex differences in brain have raised a lot of controversy 
                         due to small and sometimes contradictory sex-specific 
                         effects. Many studies have been published to identify 
                         sex-biased genes in one or several brain regions but 
                         the assessment of the robustness of these studies is 
                         missing. We therefore collected huge amount of publicly
                         available transcriptomic data to first estimate whether
                         consistent sex differences exists and further explore 
                         their likely origin and functional significance."),
                      h5("Results and Conclusion"),
                      h6("In order to systematically characterise sex specific 
                         differences across human brain regions, we collected 
                         transcription profiles for more than 16000 samples from 
                         over 46 datasets across eleven brain regions. By systematic 
                         integration of the data from multiple studies, we identified 
                         robust transcription level differences in human brain across eleven 
                         brain regions and classified male-biased and female-biased genes. 
                         Firstly, both male and female-biased genes were highly conserved across
                         primates and showed a high overlap with sex-biased genes in other species. 
                         Female-biased genes were enriched for neuron-associated processes while 
                         male-biased genes were enriched for membranes and nuclear structures. 
                         Male-biased genes were enriched on the Y chromosome while female-biased 
                         genes were enriched on the X chromosome, which included X chromosome inactivation 
                         escapees explaining the origins of some sex differences. We noted that age is 
                         a major co-variate of sex-differences. Finally, many more female-biased genes 
                         were affected by adverse drug reactions than male-biased genes."),
                      h6("In summary, by building a comprehensive resource of sex differences across 
                         human brain regions at gene expression level, we explored their likely origin and 
                         functional significance. We have also developed a web resource to make the entire
                         analysis available for the scientific community for further exploration."),


                      h4("Source code"),
                      h6("The open source code for this application is available on Github :"),
                      tags$a(href = "https://github.com/PattaWapee/SexRankBrain", "https://github.com/PattaWapee/SexRankBrain")
                      



                      




            )

  )
}
