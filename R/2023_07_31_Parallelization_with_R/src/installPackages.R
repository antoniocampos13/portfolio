packages <-
  c(
    "here",
    "tidyverse",
    "glue",
    "openxlsx",
    "future",
    "future.callr",
    "devtools",
    "tictoc",
    "BiocManager"
  )

lapply(packages, function(pkg) {
  if (!pkg %in% rownames(installed.packages()))
    
    install.packages(pkg)
  
})

biocPackages <-
  c(
    "edgeR",
    "AnnotationDbi",
    "annotate",
    "org.Hs.eg.db",
    "EnsDb.Hsapiens.v79",
    "ensembldb"
    )

lapply(biocPackages, function(pkg) {
  if (!pkg %in% rownames(installed.packages()))
    
    BiocManager::install(pkg)
  
})
