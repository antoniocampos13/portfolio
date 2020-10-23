eisenberg <- read.table(here("data", "HK_genes.txt"), header = TRUE, sep = "\t")

eisenberg_data <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = eisenberg$gene_name, # 3804 genes
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL"
)

eisenberg_data <- eisenberg_data %>% tidyr::drop_na(ENTREZID) # 3558 genes, dropped 246

housekeeping <- as.character(eisenberg_data$ENTREZID)
