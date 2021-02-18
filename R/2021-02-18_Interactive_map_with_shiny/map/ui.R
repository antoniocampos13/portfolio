ui <- fluidPage(
  titlePanel("Toy Dataset: Variants in Brazil"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("gene", label = "Choose a gene", names(genes_variants)),
      selectInput("variant", label = "Choose a variant", genes_variants[[1]]),
    ),
    
    
    mainPanel(
      h4("The circles' diameters are proportional to sample size. Hover them with the cursor to see allelic frequencies. Click to see sample size (number of alleles)."),
      
      leafletOutput(outputId = "map")
    )
  )
)
