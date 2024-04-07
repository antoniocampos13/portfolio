# Download data from GSE163836 project (Peluffo et al. 2023) (train dataset)
URL_TRAIN <-
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE163836&format=file&file=GSE163836%5FIBC%5Fsc10x%5Fraw%5Fcounts%2ERData%2Egz"

train_data_path <- here("input", "raw_count_matrix_train.RData.gz")

if (file.exists(train_data_path)) {
  load(train_data_path)
} else {
  download.file(URL_TRAIN, destfile = train_data_path)
  
  load(train_data_path)
}

# Prepare single-cell data labels (train dataset) ----
train_labels <-
  ifelse(str_detect(rawCounts.grp, "PR"),
         1,
         0)

## Import raw single-cell transcript counts ----
train_counts <-
  CreateSeuratObject(counts = rawCounts,
                     project = "GSE163836")

## Transform, normalize, and reduce dimensionality ----
### https://satijalab.org/seurat/articles/sctransform_vignette
train_counts_transformed <- train_counts %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  SCTransform(vars.to.regress = "percent.mt") %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()

# Save scaled data for future use ----
train_counts_transformed_scaled <- train_counts_transformed[["SCT"]]$scale.data

train_gene_names <- rownames(train_counts_transformed_scaled)

train_counts_transformed_scaled %>% 
  as_tibble() %>% 
  mutate(index = train_gene_names) %>% 
  relocate(index) %>%
  write_tsv(here("intermediate", "train_counts_transformed_scaled.tsv.gz"))

# Save metadata (data labels) ----
train_metadata = tibble(
  umi = names(as_tibble(train_counts_transformed_scaled)),
  original_label = rawCounts.grp,
  label = train_labels
)

train_metadata %>% write_tsv(here("intermediate", "train_metadata.tsv"))

rm(rawCounts, rawCounts.grp)

save(train_counts_transformed, train_counts_transformed_scaled, train_gene_names, file = here("intermediate", "train_backup.RData"))
