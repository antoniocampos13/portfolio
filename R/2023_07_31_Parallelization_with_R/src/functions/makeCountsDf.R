library(tidyverse)

makeCountsDf <- function(dataFrame, columnIndexes) {
  
  counts <- dataFrame %>%
    dplyr::select(1, all_of(columnIndexes)) %>%
    rename_with(~ str_replace(.x, "C9ALS|C9FTD", "case")) %>%
    rename_with(~ str_replace(.x, "Control", "control")) %>%
    as.data.frame() 
  
  row.names(counts) <- counts$gene_ID
  
  counts <- subset(counts, select = -c(gene_ID))
  
  return(counts)
}
