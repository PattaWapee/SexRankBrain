library(ggplot2)
library(reshape2)

source("Fun_rank.R")
source("enrich_cluster.R")

plt_num_fil <- function(num_fil_df,output_file){
  plot_num <- ggplot(num_fil_df, aes(x=dataset, y=value, fill=gender)) +
    geom_bar(stat="identity",position='dodge') +
    scale_fill_manual(values=c("deeppink", "dodgerblue"))+
    theme(axis.text.x = element_text(angle = 45))    
  
  png(output_file, units="in",
      width=4, height=5, res = 300)
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

runRRA <- function(FC, DE_pval, RRA_pval){
    #1 Run RRA and filter from pvalue
    RRA_F <-cal_filter_RRAreg(FC, DE_pval, RRA_pval,"F","Healthy")
    RRA_M <-cal_filter_RRAreg(FC, DE_pval, RRA_pval,"M","Healthy")
     # 2. export RRA filtered table to Rdata object
    save(RRA_F, RRA_M, file = paste0("output_obj/RRA_DE_",DE_pval,
                                     "_FC_", FC,
                                     "_RRA_",RRA_pval,
                                     ".RData"))
    return(list(RRA_F, RRA_M))
}


runplot <- function(FC, DE_pval, RRA_pval, RRA_ls){
    RRA_F <- RRA_ls[[1]]
    RRA_M <- RRA_ls[[2]]
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

run_enrich <- function(FC, DE_pval, RRA_pval, RRA_ls){
  RRA_F <- RRA_ls[[1]]
  RRA_M <- RRA_ls[[2]]
  
  output_name  <-  paste0("DE_",DE_pval,"_FC_",FC,"_RRA_",RRA_pval)
  # GO enrich female
  F_gene_ls <- lapply(RRA_F, function(x) x$Gene )
  plot_GO(F_gene_ls,'BP','female',12, output_name)
  plot_GO(F_gene_ls,'CC','female',8, output_name)
  plot_GO(F_gene_ls,'MF','female',8, output_name)
  
  # GO enrich male
  M_gene_ls <- lapply(RRA_M, function(x) x$Gene )
  plot_GO(M_gene_ls,'BP','male',12, output_name)
  plot_GO(M_gene_ls,'CC','male',8, output_name)
  plot_GO(M_gene_ls,'MF','male',8, output_name)
  # KEGG
  plot_KEG(F_gene_ls,'female', output_name)
  plot_KEG(M_gene_ls,'male', output_name)

  # DisGeNET
  plot_Dis(RRA_F, 'female', output_name)
  plot_Dis(RRA_M, 'male', output_name)
}

#run_main <- function(region_num, FC, DE_p_value, output_path){

FC = c(1.2,1.5,2.0)
DE_p_value = c(0.05,0.1)
RRA_p_value = c(0.05,0.1)
combi  <-  expand.grid(FC, DE_p_value, RRA_p_value)

for (i in (1:nrow(combi) ) ){
  RRA_ls <- runRRA(combi[i,1], combi[i,2], combi[i,3])
  runplot(combi[i,1],combi[i,2], combi[i,3], RRA_ls)
  run_enrich(combi[i,1], combi[i,2], combi[i,3], RRA_ls)
}



