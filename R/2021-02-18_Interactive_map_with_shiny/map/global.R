load("data/map_data.RData")

load("data/genes_variants.RData")

if(!require(leaflet)){
  install.packages(c("leaflet", "leaflet.extras"))
  library(leaflet)
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

