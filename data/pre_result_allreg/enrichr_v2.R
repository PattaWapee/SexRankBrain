library('ggplot2')
library('stringr')
library("enrichR")
setEnrichrSite("Enrichr") # Human genes

###################################
# Enrich function for all package #
###################################


#-------------#
####enrichr####
#-------------#

# list of database enrichr
# dbs <- listEnrichrDbs()
# ("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
# ("DisGeNET")
# ("KEGG_2021_Human")
enriched <- function(df,dbsx){
  gene_ls <- as.vector(df[,1])
  enriched <- enrichr(gene_ls, dbsx)
  # add gene_count column in df
  endf <- enriched[[1]]
  gene_count <-  unlist(lapply(endf$Overlap, function(x) { unlist(strsplit(x, "/"))[1] } ))
  endf$gene_count <- gene_count

  return(endf)
}

plt_enriched_reg <- function(df,dbsx){
  enrich_df <- enriched(df,dbsx)
  plotEnrich(enrich_df, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}


#--------------#
###disgenet2r###
#--------------#

#library(devtools)
#install_bitbucket("ibi_group/disgenet2r")
library(disgenet2r)
Sys.setenv(DISGENET_API_KEY='f327ba3156d80b7d2c9e4482cefffcbd6e688dc3' )

enrich_disgenet2r <- function(df){
  gene_ls <- as.vector(df[,1])
  res_enrich <-disease_enrichment( entities =gene_ls, 
                                  vocabulary = "HGNC", 
                                  database = "CURATED")
  
  res <- extract(res_enrich)
  return(res)
}


#########################################
# Function fo enrichment in all regions #
#########################################


enriched_allreg <- function(df_ls, package, dbs){
  if(package == 'enrichr'){
    enrich_ls <- lapply(df_ls, function(x) enriched(x,dbs))
  }
  else if(package == 'disgenet2r'){
    enrich_ls <- lapply(df_ls, function(x) enrich_disgenet2r(x))

  }

  # Add region column
  for (reg in names(enrich_ls)) {
    df <- enrich_ls[[reg]]
    if  ( !is.null(dim(df)) ){
      df$region <- rep(reg,nrow(df))
      enrich_ls[[reg]] <- df
    }
  }
  return(enrich_ls)
}

## choose only top5 enrich to combine in one df for plotting
top5_df <- function(df_ls){
  head_df <- lapply(df_ls, function(x) head(x,5))
  combind_df <- Reduce(rbind, head_df)
  # Add gene ratio column
  #combind_df$gene_ratio <- sapply(combind_df$Overlap, function(x) eval(parse(text = x)) )
  return(combind_df)
}

# Run & Plot function
run_plt_enrich <- function(df_ls,dbsx,package,type){
  enrich_df <- enriched_allreg(df_ls,package,dbsx)
  top5 <- top5_df(enrich_df)
  write.table(top5,paste0("output_obj/top5_",dbsx,"_",type,".txt"))
  # filter top df with only p.adjust < 0.05
  if(package == 'enrichr'){
    filtered_top <- top5[top5['Adjusted.P.value']< 0.05,]
    filtered_top$gene_count <- as.numeric(filtered_top$gene_count)
    filtered_top <- filtered_top[filtered_top['gene_count'] > 1,]
    size_var = 'gene_count'
    color_var = 'Adjusted.P.value'
    y_var = 'Term'
  }
  else if(package == 'disgenet2r'){
    filtered_top <- top5[top5['FDR']< 0.05,]
    size_var = 'Count'
    color_var = 'FDR'
    y_var = 'Description'
  }
  Title = (paste0( "Enrichment of ",type, "-biased genes in ",dbsx ))
  p <- ggplot(filtered_top, aes_string(x = 'region', y = y_var, size = size_var)) + 
    geom_point(aes_string(color = color_var)) + 
    scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=TRUE)) +
    ylab(NULL) + ggtitle(Title) + scale_size_continuous(range=c(3, 8)) + 
    guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2))+
    scale_y_discrete(labels=function(x) str_wrap(x,width=40))

  return(p)
}

run_plt_enrich2 <- function(df_ls,dbsx,type){
  enrich_df <- enriched_allreg(df_ls,dbsx)
  top5_plt <- top5_df(enrich_df)
  if(type=='hflm'){
    low_col <- "purple"
    high_col <- "deeppink"
  }
  else{
    low_col <- "darkblue"
    high_col <- "lightseagreen"
  }
  p <- ggplot(top5_plt,aes(region,Term))+
    geom_point(aes(size= gene_ratio,color = P.value))+
    scale_color_gradient(low=low_col,high =high_col) +
    geom_text(aes(region, Term, label = gene_count), colour = I(alpha("black", 0.85)), size = 2 ) +
    ggtitle(paste0( "Enrichment of ",type, "-variable genes in ",dbsx ))+
    labs(caption = "Each bubble point are labeled with gene count")

  return(p)
  
}
