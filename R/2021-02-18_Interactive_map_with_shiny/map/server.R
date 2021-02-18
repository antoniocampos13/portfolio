server <- function(input, output, session) {
  
  observe({
    updateSelectInput(session, "variant", choices = genes_variants[[input$gene]])
  })

  data_subset <- reactive({
    df <- map_data %>% filter(gene == input$gene & variant == input$variant)

    return(df)
  })

  output$map <- renderLeaflet({
    leaflet(data_subset()) %>%
      setView(lat = -14.235004, lng = -51.92528, zoom = 4) %>%
      addTiles() %>%
      addCircles(
        lat = ~lat,
        lng = ~long,
        weight = 1,
        radius = ~ sqrt(alleles_total) * 5000,
        popup = ~ as.character(paste0("Alleles: ", alleles_total)),
        label = ~ as.character(paste0(input$variant, " Allele frequency: ", round(freq, digits = 2))),
        fillOpacity = 0.5
      )
  })
}
