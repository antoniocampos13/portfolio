# Install packages
install.packages(c("here", "tidyverse", "openxlsx", "BiocManager"))

## Install Bioconductor packages
BiocManager::install(c("edgeR", "AnnotationDbi", "annotate", "org.Hs.eg.db", "EnsDb.Hsapiens.v79","ensembldb"))

# Load packages
library(here)

# Shortcut: loading counts data frame from disk
source(here("data", "counts.RData"))

# Set up the output path string
out_path <- here("output", "prostate_cancer.xlsx")

# Run the custom function
## Will load all
edger_setup("prostate_cancer", counts, replicates = TRUE, filter = TRUE, gene_id = "ENSEMBL", out_path)
