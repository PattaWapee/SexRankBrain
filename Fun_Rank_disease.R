library(reshape2)
library(ggplot2)
source('global.R')


Filter_gene <- function( order.gene.df, P.val.cutoff,
                         FC.val.cutoff, gender ) {
  if(gender =='F'){
    logFC <- log2(1/as.numeric(FC.val.cutoff))
    gene.sig <- order.gene.df[  order.gene.df[["P.Value"]] <= P.val.cutoff
                                & order.gene.df[["logFC"]] <= logFC, ]
  }
  
  else{
    # Ordering df by logFC from positive to negative (male to female biased)
    M.order.gene.df <-  order.gene.df[order(-order.gene.df[["logFC"]]),]
    logFC <- log2(as.numeric(FC.val.cutoff))
    gene.sig <- M.order.gene.df[  M.order.gene.df[["P.Value"]] <= P.val.cutoff
                                & M.order.gene.df[["logFC"]] >= logFC, ]
  }
  
  #If there are sig genes add index number for each sig gene
  
  if(nrow(gene.sig) > 0) {
    gene.sig$index <- seq.int(nrow(gene.sig))
  }
  #print(nrow(gene.sig))
  
  return(gene.sig)
}

# This function is to import all DE topTable in each brain regions
import_topTable_reg_dis <- function(reg_num,disease){
  reg_list <- c('AMY','CBC', 'CC', 'FC', 'HIP','MED', 'OC', 'PL','STR','TC', 'THA')
  reg_str = reg_list[as.numeric(reg_num)]
  #print(reg_num)
  
  library("stringr")
  #print(reg_str)
  path_reg = paste0('data/',disease,'/',reg_str)
  files <- list.files(path = path_reg,pattern = ".txt$",full.names = TRUE)
  name_files <- list.files(path = path_reg,pattern = ".txt$",full.names = FALSE)
  name.files <- str_replace(name_files, ".txt","")
  #print(name.files)
  #import data table from file list
  topTable.dat <-  lapply(files, read.table)
  #ordering genes by logFC (negative to positive) (female to male biased)
  topTable.order <- lapply(topTable.dat, function(x) x[order(x$logFC),])
  
  #print(head(topTable.dat[[1]]))
  # add gene column
  topTable.order <- lapply(topTable.order, function(x){x$Gene <- rownames(x); return(x)})
  names(topTable.order) <- name.files
  return(topTable.order)
}

## Import & filter sig genes
Filter_gene_reg_dis <- function(reg_num, P.val.cutoff, FC.val.cutoff, gender,disease){
  
  topTable_reg <- import_topTable_reg_dis(reg_num,disease)
  sigTable_reg <- lapply(topTable_reg, function(x) Filter_gene(x,P.val.cutoff, FC.val.cutoff, gender))
  #print(head(sigTable_reg))
  return(sigTable_reg)
  #return(topTable_reg)
  
}

numfil_gene <- function(sigTable_reg_F, sigTable_reg_M){
  
  num_F <- lapply(sigTable_reg_F, nrow)
  print(num_F)
  num_M <- lapply(sigTable_reg_M, nrow)
  
  num_df <- data.frame(
    Female = unlist(num_F),
    Male = unlist(num_M) )
  
  num_df2 <- melt(as.matrix(num_df), id="dataset")
  colnames(num_df2) <- c('dataset','gender','value')
  
  return(num_df2)
  
} 


## RRA function #####

RA.analysis <- function(gene.list){
  print(str(gene.list))
  library("RobustRankAggreg")
  rank.RRA <- aggregateRanks(rmat = rankMatrix(gene.list),
                             method = "RRA")
  colnames(rank.RRA) <- c("Gene", "P.value")
  rank.RRA$P.adj <- p.adjust(rank.RRA$P.value, method = "bonferroni" )
  
  return(rank.RRA)
}

######################################################
######## Function for across brain region analysis ###
######################################################


RRA_allreg <- function(input_df,gender, disease){
  
  reg_list <- c('AMY','CBC', 'CC', 'FC', 'HIP','MED', 'OC', 'PL','STR','TC', 'THA')
  
  
  RRA_table_list <- list()
  regions <- input_df$Region
  
  for (i in 1:length(regions)){
    reg <- regions[i]
    
    DE_p_value <- input_df[input_df$Region==reg,]$DE_pval
    FC <- input_df[input_df$Region==reg,]$FC
    
    
    # 1. import and filtered genes for each dataset
    reg_num = match(reg,reg_list)
    sig_reg <- Filter_gene_reg_dis(reg_num,DE_p_value,FC,gender,disease)
    #print(sig_reg[1])
    
    # 2. Rank aggregation of filtered genes with RRA pvalue cutoff
    sig <- lapply(sig_reg , rownames)
    if (all(sapply(sig,function(x) identical(x,character(0))))) {
      print(paste0("There are no significant sex-biased genes in ",reg))
      # export emthy df if no sig genes
      RRA <- data.frame(Gene=character(),
                        P.value=numeric(), 
                        P.adj=numeric(), 
                        stringsAsFactors=FALSE) 
    }
    else {
      RRA <- RA.analysis(sig)
      }
    
    
    RRA_table_list <- append(RRA_table_list, list(RRA) )
    
  }
  
  return(RRA_table_list)
}

##### 2. Functions for filtering RRA gene with cutoff 

Filter_RRA_reg <- function(RRA_df, Pval){
  RRA.sig <- RRA_df[  RRA_df[["P.value"]] <= Pval, ]
  return(RRA.sig)
}

cal_filter_RRAreg <- function(inputcutoff, gender,disease){
  print("Calculate RRA of all regions for")
  print(gender)
  RRAreg <- RRA_allreg(inputcutoff,gender,disease)
  RRApval <- inputcutoff$RRA_pval
  filter_mat <- mapply(function(X,Y) { Filter_RRA_reg(X,Y)}, X=RRAreg, Y=RRApval)
  #making dataframe from matrix after filtering
  filter_df <- list()
  for (i in 1:dim(filter_mat)[2]){
    df <- data.frame(filter_mat["Gene",i],filter_mat["P.value",i], filter_mat["P.adj",i])
    #print(head(df))
    filter_df <-  append(filter_df, list(df) )
    
  }
  names(filter_df) <- inputcutoff$Region
  return(filter_df)
}


#### 4. Function for plotting correlation heatmap 


#reference: https://github.com/coriell-research/coriell/blob/master/R/list-to-matrix.R
#create binary matrix of intersection genes
list_to_matrix <- function(sets) {
  stopifnot("List of vectors must be supplied" = class(sets) == "list")
  union_all <- Reduce(union, sets)
  
  if (sum(is.na(union_all)) > 0) {
    message("NA values present in union of all sets. NA values will be dropped in final matrix")
    union_all <- union_all[!is.na(union_all)]
  }
  
  mat <- matrix(
    data = 0,
    nrow = length(union_all),
    ncol = length(sets)
  )
  colnames(mat) <- names(sets)
  rownames(mat) <- union_all
  
  for (i in seq_along(sets)) {
    mat[unique(sets[[i]][!is.na(sets[[i]])]), i] <- 1
  }
  mat
}


correlation_reg <- function(RRA_ls){
  #sig_gene_ls <- lapply(RRA_ls, function(x) x$Gene)
  sig_gene_ls <- lapply(RRA_ls, function(x) as.vector(x[,1]))
  # create intersect binary matrix
  sig_mat <- list_to_matrix(sig_gene_ls)
  # Remove columns with all zeros
  sig_mat2 <- sig_mat[, !sapply(colnames(sig_mat), function(col) {all(sig_mat[,col]==0) })]
  sig_cor <- cor(sig_mat2,method = "spearman")
  #print(sig_cor)
  return(sig_cor)
} 


