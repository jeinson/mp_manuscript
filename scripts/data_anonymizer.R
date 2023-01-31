# Data anonymyzer. Read in data and get rid of individual labels, for sake
# of making this more difficult to track. 

rv_bytiss <- readRDS("data/fig2_data/rare_variants_exon_psi_zscore_bytissue.rds")
tissues <- names(rv_bytiss)

for(tiss in tissues){
  x <- rv_bytiss[[tiss]]
  
  x$indv <- "XXXXX"
  
  x$var <- "XXXXX"
  x$chr <- "XXXXX"
  x$pos <- "XXXXX"
  x$ref <- "XXXXX"
  x$alt <- "XXXXX"
  
  rv_bytiss[[tiss]] <- x
}

saveRDS(rv_bytiss, "data/fig2_data/rare_variants_exon_psi_zscore_bytissue_anon.rds")
