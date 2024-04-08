packages <-
  c(
    "here",
    "tidyverse",
    "glue",
    "reticulate",
    "R.utils",
    "Seurat"
  )

lapply(packages, function(pkg) {
  if (!pkg %in% rownames(installed.packages()))
    
    install.packages(pkg)
  
})