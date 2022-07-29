library("enrichR")
setEnrichrSite("Enrichr") # Human genes
#dbs <- listEnrichrDbs()

# list of database enrichment
# ("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
# ("DisGeNET")
# ("KEGG_2021_Human")

enriched <- function(df,dbsx){
  gene_ls <- as.vector(df[,1])
  enriched <- enrichr(gene_ls, dbsx)
  
  return(enriched[[1]])
}

plt_enriched_reg <- function(df,dbsx){
  enrich_df <- enriched(df,dbsx)
  plotEnrich(enrich_df, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}

# Function fo enrichment in all regions
enriched_allreg <- function(df_ls,dbsx){
  enrich_ls <- lapply(df_ls, function(x) enriched(x,dbsx))
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
  combind_df$gene_ratio <- sapply(combind_df$Overlap, function(x) eval(parse(text = x)) )
  return(combind_df)
}

# Run & Plot function
run_plt_enrich <- function(df_ls,dbsx,gender){
  enrich_df <- enriched_allreg(df_ls,dbsx)
  top5_plt <- top5_df(enrich_df)
  if(gender=='F'){
    low_col <- "purple"
    high_col <- "deeppink"
  }
  else{
    low_col <- "darkblue"
    high_col <- "lightseagreen"
  }
  p <- ggplot(top5_plt,aes(region,Term))+geom_point(aes(size= gene_ratio,color = P.value))+scale_color_gradient(low=low_col,high =high_col)
  return(p)
  
}
