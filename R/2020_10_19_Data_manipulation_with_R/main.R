# Load packages ----
library(RPostgres)
library(here)
library(tidyverse)

# Set up connection to PostgreSQL  ----
# Add credentials to a .Renviron file and reload R session
con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "tcga",
  host = "localhost",
  port = 5432,
  user = Sys.getenv("userid"),
  password = Sys.getenv("pwd")
)

# Retrieve gene counts and cases table from the database ----
cases <- dbGetQuery(con, "SELECT * FROM gene_counts_cases")

# Pivot table to wide format ----
# Some samples had more than one gene expression quantification. I took the mean of those values
# NA values were substituted for 0
cases_pivoted <-
  cases %>% pivot_wider(
    names_from = case_id,
    values_from = gene_count,
    values_fill = 0,
    names_repair = "check_unique",
    values_fn = mean
  )

# Round eventual decimal numbers----
counts <-
  cases_pivoted %>% mutate(across(where(is.numeric), round, 0))


# Retrieve deduplicated follow_up table----
dbExecute(
  con,
  "CREATE TABLE followup_dedup AS SELECT case_id, STRING_AGG(followup_primarytherapyoutcomesuccess_1, ',') AS outcome FROM follow_up GROUP BY case_id"
)

dbExecute(
  con,
  "CREATE TABLE outcomes AS SELECT case_id, outcome FROM allcases INNER JOIN followup_dedup USING(case_id)"
)

outcomes <- dbGetQuery(con, "SELECT * FROM outcomes")

# Recode outcomes to a simpler case/control classification----
# Complete Remission/Response = control; otherwise = case
outcomes <- outcomes %>% mutate(class = case_when(
  str_detect(outcome, "Complete") ~ "control",
  str_detect(outcome, "Partial") ~ "case",
  str_detect(outcome, "Disease") ~ "case"
))

# Simplify case ID number and label them as case or control----
outcomes <-
  outcomes %>% mutate(new_names = paste0(str_sub(case_id, -12), "_", class))


# Apply new labels to the dataframe----
colnames(counts) <- dplyr::recode(colnames(counts), !!!setNames(as.character(outcomes$new_names), outcomes$case_id))

# Convert gene_ids into row names, then delete gene_id column----
counts <- as.data.frame(counts)

row.names(counts) <- counts$gene_id

counts <- subset(counts, select = -c(gene_id))

# Check data frame dimensions (shape) ----
dim(counts)

# Save counts data frame to use later
save(counts, file = "counts.RData")
