library(here)
library(tidyverse)
library(glue)
library(Seurat)
library(reticulate)

folders <-
  c("src", "input", "intermediate", "output", "output/xgboost", "output/plot")

lapply(folders, function(folder) {
  path <- do.call(here, as.list(folder))
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
})

# Scripts ----
## Prepare train data ----
source(here("src", "prepare_train_data.R"))

## Prepare test data ----
source(here("src", "prepare_test_data.R"))

## Gene (features) refinement ----
final_gene_list <-
  tibble(id = rownames(test_counts_transformed_scaled)) %>%
  inner_join(test_gene_names_df, by = "id") %>%
  inner_join(tibble(index = train_gene_names), by = "index")

final_gene_list %>% select(index) %>% write_tsv(here("intermediate", "final_gene_list.tsv"))

# Run XGBoost on Python ----
source_python(here("src", "run_xgboost.py"))

# Make plots ----
source(here("src", "make_plots.R"))
