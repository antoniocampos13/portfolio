library(dplyr)
library(here)
library(openxlsx)
library(ids)
library(maps)

GENES <- 3
VARIANTS <- 2
STATES <- 27
DATAPOINTS <- GENES * VARIANTS * STATES  

set.seed(123)
alleles_count <- sample(c(25:80), DATAPOINTS, replace = TRUE)
alleles_total <- 2 * sample(c(50:250), STATES, replace = TRUE)
genes <- random_id(GENES, 4, use_openssl = FALSE)
var_ids <- random_id(GENES * VARIANTS, 6, use_openssl = FALSE)

br <- read.xlsx("states_capitals.xlsx")

alleles_states <- bind_cols(state = br$state, alleles_total = alleles_total)

coords <- world.cities %>% filter(country.etc == "Brazil") %>% 
  filter(name %in% br$capital) %>%
  filter(lat != -26.48) %>%
  mutate(state = br$state) %>% 
  select(lat, long, state)

gene_var_comb <- bind_cols(gene = rep(genes, each = VARIANTS), variant = var_ids)

combinations <- expand.grid(state = br$state, variant = gene_var_comb$variant) %>% 
  inner_join(gene_var_comb, by = "variant") %>%
  inner_join(alleles_states, by = "state")

map_data <- bind_cols(combinations, 
                      alleles_count = alleles_count) %>%
  inner_join(coords, by = "state") %>%
  mutate(freq = alleles_count / alleles_total)

save(map_data, file = here("map","data","map_data.RData"))

gene_list <- map_data %>% select(gene, variant) %>%
  distinct(gene, variant) %>%
  group_by(gene) %>%
  mutate(varstring = paste0(variant, collapse = ",")) %>%
  select(-variant) %>%
  distinct(gene, varstring)

genes_variants <- lapply(seq_along(gene_list$gene), function(i)
  
  unlist(strsplit(gene_list$varstring[i], ","))
  
)

names(genes_variants) <- gene_list$gene

save(genes_variants, file = here("map","data","genes_variants.RData"))
