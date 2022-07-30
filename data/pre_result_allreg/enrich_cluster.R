library("clusterProfiler")

#GO.ont = 'CC', 'BP', 'MO'
compareGO <- function(input.ls, GO.ont){
  enrich <- compareCluster(geneCluster  =input.ls, 
		       fun          = "enrichGO",
		       keyType      = "SYMBOL",
		       ont          = GO.ont,
		       OrgDb        = "org.Hs.eg.db",
		       pAdjustMethod= "BH",
		       pvalueCutoff = 0.05,
		       qvalueCutoff = 0.05)
	return(enrich)
}

plot_GO <-  function(input_ls, GO.ont, gender, font_size, output_file){
  en.GO <- compareGO(input_ls,GO.ont)
  GO.plot <- dotplot(en.GO,by = 'count',
        title = paste0(GO.ont," of ",gender,"-biased genes"), font.size = font_size)
  #tiff(paste0("plt_enrich_img/",gender,"_GO_",GO.ont,".tiff"), units="in", width=7, height=7.5, res=300)
  tiff(output_file, units="in", width=7, height=7.5, res=300)
  print(GO.plot)
  dev.off()
  ego.df <- as.data.frame(en.GO)
  write.table(ego.df,paste0("output_obj/",gender,"_GO_",GO.ont,".txt"))
  cnet <- cnetplot(en.GO)
  tiff(paste0("output_plot/",gender,"_GO_",GO.ont,"_cnet.tiff"), units="in", width=7, height=7.5, res=300)
  print(cnet)
}
