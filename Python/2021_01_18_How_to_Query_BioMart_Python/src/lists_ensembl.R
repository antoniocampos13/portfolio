library(here)
library(biomaRt) # install with BiocManager::install("biomaRt")
library(openxlsx)

# Modify the output path as needed. Check here package documentation
output_path <- here("lists_ensembl.xlsx")

# Connecting to Ensembl
# The useast mirror apparently is more stable than the primary site
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "useast") 

# Retrieving Ensembl descriptors
datasets <- listDatasets(ensembl)
filters <- listFilters(ensembl)
Attributes <- listAttributes(ensembl)

# Saving data frames to a spreadsheet
wb <- createWorkbook()

addWorksheet(wb, sheetName = "Datasets", gridLines = TRUE)
writeData(wb, sheet = "Datasets", x = datasets)

addWorksheet(wb, sheetName = "Attributes", gridLines = TRUE)
writeData(wb, sheet = "Attributes", x = Attributes)

addWorksheet(wb, sheetName = "Filters", gridLines = TRUE)
writeData(wb, sheet = "Filters", x = filters)

saveWorkbook(wb, output_path, overwrite = TRUE)
