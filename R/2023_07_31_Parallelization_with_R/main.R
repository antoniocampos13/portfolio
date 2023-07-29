source("src/installPackages.R")

library(here)

# Paths ----
if (!dir.exists(here("outputs"))) { dir.create(here("outputs"))}

# Scripts ----
source(here("src", "functions.R"))

source(here("src", "prepareData.R"))

# Run DEA  ----
nWorkers <- c(1, 2, 4, 6)

lapply(nWorkers, runParallelDEA)
