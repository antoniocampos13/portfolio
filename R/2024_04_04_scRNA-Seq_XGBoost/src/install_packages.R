packages <-
  c(
    "here",
    "tidyverse",
    "glue",
    "reticulate",
    "Seurat"
  )

lapply(packages, function(pkg) {
  if (!pkg %in% rownames(installed.packages()))
    
    install.packages(pkg)
  
})