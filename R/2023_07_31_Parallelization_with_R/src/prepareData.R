# Load data ----
countDataOriginal <- read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE219nnn/GSE219278/suppl/GSE219278_allSamples_rsem_genes_results_counts_annotated.tsv.gz")

# Round up gene counts ----
countData <- countDataOriginal %>%
  mutate(across(where(is.numeric), ceiling))

# Identify samples and outcomes ----
sampleInfo <- tibble(Names = names(countData)[4:length(names(countData))]) %>%
  separate(Names, into = c("outcome", "patientId", "location", "cellType")) %>%
  mutate(columnIndex = row_number() + 3)

alsDfsList <- sampleInfo %>% 
  dplyr::filter(outcome == "C9ALS" | outcome == "Control") %>%
  group_by(location, cellType) %>%
  group_split()

ftdDfsList <- sampleInfo %>% 
  dplyr::filter(outcome == "C9FTD" | outcome == "Control") %>%
  group_by(location, cellType) %>%
  group_split()

allAnalysesList <- c(alsDfsList, ftdDfsList)

# Split count data into 12 dataframes for differential expression analysis with edgeR ----
countsDfsList <- lapply(seq_along(allAnalysesList), function(index) {
  
  df <- allAnalysesList[[index]]
  
  cIdx <- df$columnIndex
  
  makeCountsDf(dataFrame = countData, columnIndexes = cIdx)
  
})

allAnalysesNames <- sampleInfo %>%
  dplyr::filter(outcome != "Control") %>%
  dplyr::select(outcome, location, cellType) %>%
  distinct() %>%
  mutate(names = glue("{outcome}_{location}_{cellType}")) %>%
  pull(names)

names(countsDfsList) <- allAnalysesNames
