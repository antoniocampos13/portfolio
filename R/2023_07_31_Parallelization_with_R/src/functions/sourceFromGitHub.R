
URLs <- c("https://raw.githubusercontent.com/antoniocampos13/portfolio/master/R/2020_10_22_DEA_with_edgeR/src/edgeR_setup.R")

lapply(URLs, devtools::source_url)
