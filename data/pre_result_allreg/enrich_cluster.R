library("clusterProfiler")
source('enrichr_v2.R')

options(ggrepel.max.overlaps = Inf)


##########################
#GO.ont = 'CC', 'BP', 'MO'
##########################
compareGO <- function(input.ls, GO.ont){
  enrich <- compareCluster(geneCluster  =input.ls, 
		       fun          = "enrichGO",
		       keyType      = "SYMBOL",
		       ont          = GO.ont,
		       OrgDb        = "org.Hs.eg.db",
		       pAdjustMethod= "BH",
           minGSSize    = 1,
		       pvalueCutoff = 0.05,
		       qvalueCutoff = 0.05)
	return(enrich)
}

plot_GO <-  function(input_ls, GO.ont, gender, font_size, output_name){
  en.GO <- compareGO(input_ls,GO.ont)
  GO.plot <- dotplot(en.GO,by = 'count',
        title = paste0(GO.ont," of ",gender,"-biased genes"), font.size = font_size)
  tiff(paste0("plt_enrich_img/",output_name,"_",GO.ont,"_",gender,".tiff"), 
       units="in", width=7, height=7.5, res=300)
  print(GO.plot)
  dev.off()
  ego.df <- as.data.frame(en.GO)
  write.table(ego.df,
              paste0("enrich_table/",
                     output_name,"_",GO.ont,"_",gender,".txt"))
}


#########################
#####   KEGG    #########
#########################

compareKEGG <- function(input.ls){
	input.ls.keg <- lapply(input.ls, function(x) 
			      {gene.df <- GenetoENTREZ(x)
			      return(gene.df$ENTREZID)})
	length(input.ls.keg)
	names(input.ls.keg) <- names(input.ls)
  enrich <- compareCluster(geneCluster  =input.ls.keg, 
		       fun          = "enrichKEGG",
		       pAdjustMethod= "BH",
		       pvalueCutoff = 0.05,
		       qvalueCutoff = 0.05)
	return(enrich)
}

GenetoENTREZ <- function(symbol){
	library(org.Hs.eg.db)
	hs <- org.Hs.eg.db
	entrez.df <- AnnotationDbi::select(hs, 
					   keys = symbol,
					   columns = c("ENTREZID", "SYMBOL"),
					   keytype = "SYMBOL")
	return(entrez.df)
}

plot_KEG <-  function(input_ls, gender, output_name){
  en <- compareKEGG(input_ls)
  enplot <- dotplot(en,by = 'count',
        title = paste0("KEGG pathway of ",gender,"-biased genes"), font.size = 12)
  tiff(paste0("plt_enrich_img/",output_name,"_",gender,"_KEGG.tiff"), units="in", width=7, height=7.5, res=300)
  print(enplot)
  dev.off()
  en.df <- as.data.frame(en)
  write.table(en.df,paste0("enrich_table/",output_name,"_",gender,"_KEGG.txt"))
}

###############
## DisGeNET ###
###############
plot_Dis <-  function(RRA_df, gender, output_name){
  Dis <- run_plt_enrich(RRA_df, "DisGeNET (CURATED)", "disgenet2r",gender)
  tiff(paste0("plt_enrich_img/",
              output_name,"_",gender,
              "_DisGeNET_curated.tiff"), 
       units="in", width=7, height=7.5, res=300)
  print(Dis)
  dev.off()
}

