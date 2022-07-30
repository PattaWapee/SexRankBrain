library(ggplot2)
library(reshape2)

source("Fun_rank.R")
#source("../../../enrichr.R")

plt_num_fil <- function(num_fil_df,output_file){
  plot_num <- ggplot(num_fil_df, aes(x=dataset, y=value, fill=gender)) +
    geom_bar(stat="identity",position='dodge') +
    scale_fill_manual(values=c("deeppink", "dodgerblue"))+
    theme(axis.text.x = element_text(angle = 45))    
  
  png(output_file, units="in",
      width=5, height=5, res = 300)
  print(plot_num)
  dev.off()
}

plt_corr <- function(RRA_df, output_file){
  dat_corr <- correlation_reg( RRA_df )
  png(output_file, units="in",
      width=5, height=5, res = 300)
  heatmap(dat_corr,scale = "row", col = hcl.colors(50))
  dev.off()

}

runplot <- function(FC, DE_pval, RRA_pval){
    #1 Run RRA and filter from pvalue
    RRA_F <-cal_filter_RRAreg(FC, DE_pval, RRA_pval,"F","Healthy")
    RRA_M <-cal_filter_RRAreg(FC, DE_pval, RRA_pval,"M","Healthy")
    #2 Calculate number of sex-biased genes"
    num_fil_df <- numfil_gene(RRA_F ,RRA_M ) 
    # 3 Plot number of sex-biased genes
    output_file  <-  paste0("plt_num_img/num_allreg_DE_",DE_pval,"_FC_",FC,"_RRA_",RRA_pval,".png")
    plt_num_fil(num_fil_df,output_file)
    # 4. Plot correlation matrix heatmap
    output_corr_M  <- paste0("plt_corr_img/corr_DE_",DE_pval,"_FC_",FC,"_RRA_",RRA_pval,"_M.png")
    plt_corr(RRA_M, output_corr_M)
    output_corr_F  <- paste0("plt_corr_img/corr_DE_", DE_pval,"_FC_",FC,"_RRA_",RRA_pval,"_F.png")
    plt_corr(RRA_F, output_corr_F)
}


#run_main <- function(region_num, FC, DE_p_value, output_path){

FC = c(1.2,1.5,2.0)
DE_p_value = c(0.05,0.1)
RRA_p_value = c(0.05,0.1)
combi  <-  expand.grid(FC, DE_p_value, RRA_p_value)

for (i in (1:nrow(combi) ) ){
  runplot(combi[i,1],combi[i,2], combi[i,3])
}



