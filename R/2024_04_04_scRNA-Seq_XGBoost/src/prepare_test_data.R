# Download data from GSE131135 project (Shu et al. 2020) (test dataset) ----
URL_TEST <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131135&format=file&file=GSE131135%5FBRD4JQ1%5Fsc10x%5Fraw%5Fcounts%2ERData%2Egz"

test_data_path <- here("input", "raw_count_matrix_test.RData")

if (file.exists(test_data_path)) {
  load(test_data_path)
} else {
  download.file(URL_TEST, destfile = glue("{test_data_path}.gz"))
  
  R.utils::gunzip(glue("{test_data_path}.gz"))
  
  load(test_data_path)
}

# Subset data to keep specific samples ----
test_counts <- CreateSeuratObject(counts = rawCounts,
                                  project = "GSE131135")

test_metadata_filtered <- tibble(
  umi = as.vector(names(test_counts$orig.ident)),
  original_label = rawCounts.grp
) %>%
  filter(original_label %in% c("SUM149DMSO", "SUM149RDMSO")) %>%
  mutate(label = ifelse(str_detect(original_label, "R"), 1, 0))

test_metadata_filtered %>% write_tsv(here("intermediate", "test_metadata_filtered.tsv"))

test_counts <- subset(test_counts, cells = test_metadata_filtered$umi)

## Transform, normalize, and reduce dimensionality ----
### https://satijalab.org/seurat/articles/sctransform_vignette
test_counts_transformed <- test_counts %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  SCTransform(vars.to.regress = "percent.mt") %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()

# Save scaled data for future use ----
test_counts_transformed_scaled <- test_counts_transformed[["SCT"]]$scale.data

test_gene_names <- rownames(test_counts_transformed_scaled)

test_gene_names_df <- as_tibble(rawCounts.ann) %>%
  rename("index" = "symbol")

test_counts_transformed_scaled %>% 
  as_tibble() %>% 
  mutate(id = test_gene_names) %>%
  left_join(test_gene_names_df, by = "id") %>%
  relocate(index) %>%
  select(-id) %>%
  write_tsv(here("intermediate", "test_counts_transformed_scaled.tsv.gz"))

rm(rawCounts.ann, rawCounts.grp, rawCounts)

save(test_counts_transformed, test_counts_transformed_scaled, test_gene_names_df, file = here("intermediate", "test_backup.RData"))
