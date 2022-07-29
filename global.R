# Dataframe of input value for filtering genes for Rubust rank aggregation (RRA)
# cutoff input for each brain region
inputdata = data.frame(Region = c('AMY','CBC', 'CC', 'FC', 'HIP','MED', 'OC', 'PL','STR','TC', 'THA'),
                       FC     = rep( c(1.2), times =11),
                       DE_pval= rep( c(0.05), times=11),
                       RRA_pval=rep( c(0.05), times=11),
                       row.names=NULL, stringsAsFactors = FALSE
                       )

