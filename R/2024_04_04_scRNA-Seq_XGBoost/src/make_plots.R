# Get list of cells according to identity ----
parental_cells <- train_metadata %>%
  filter(label == 0) %>%
  pull(umi)

ptx_resistant_cells <- train_metadata %>%
  filter(label == 1) %>%
  pull(umi)

# Set cells identities ----
train_cells <- Cells(train_counts_transformed)
Idents(object = train_counts_transformed, cells = parental_cells) <- "Parental"
Idents(object = train_counts_transformed, cells = ptx_resistant_cells) <- "Paclitaxel-resistant"

# Retrieve top 15 genes by importance metric (weight) ----
importance_df <- read_tsv(here("output", "xgboost", "feature_importance.tsv"))

## The TSV is ordered by descending weight values, take the first 15 rows
top15 <- importance_df %>% slice(1:15) %>% pull(gene)

# Make ridge and UMAP feature plots ----
RidgePlot(train_counts_transformed, features = top15, ncol = 3)
ggsave(filename = here("output", "plot", "ridge_plot.png"), width = 45, height = 20)

FeaturePlot(train_counts_transformed, features = top15)
ggsave(filename = here("output", "plot", "feature_plot.png"), width = 45, height = 20)

# Plot mitochondrial transcripts percent violin plot ----
## Save the cells identities as metadata column
train_counts_transformed[["identity"]] <- Idents(object = train_counts_transformed)

VlnPlot(train_counts_transformed, features = "percent.mt", split.by = "identity")
ggsave(filename = here("output", "plot", "percent_mt_violin_plot.png"), width = 10, height = 10)
