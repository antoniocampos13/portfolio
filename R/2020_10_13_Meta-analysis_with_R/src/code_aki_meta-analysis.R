# Install packages
# install.packages(c("meta", "here", "openxlsx"))

# Load packages
library(meta)
library(here)
library(openxlsx)

# Import spreadsheet
aki <- read.xlsx(here("data","aki_demo.xlsx"))

# Check column names
names(aki)

# Check shape (dimensions) of the data frame
dim(aki)

# Incidence meta-analysis 
aki_incidence_meta <- metaprop(event = case_aki,
                   n = sample_size,
                   studlab = paste0(authors,", ", location),
                   data = aki,
                   sm = "PLO",
                   predict = TRUE,
                   hakn = TRUE,
                   comb.fixed = TRUE,
                   comb.random = TRUE,
                   level.comb = 0.95,
                   method.tau = "ML",
                   method.bias = "linreg",
                   warn = TRUE
)

# Saving meta-analysis raw output to file
sink(here("output", "aki_incidence.txt"))
summary(aki_incidence_meta)
sink()

# Produce a Forest plot and save to file
tiff(filename = here("output", "aki_incidence.tiff"),
     compression = "lzw",
     res = 300,
     width = 3050,
     height = 1750)

meta::forest(aki_incidence_meta,
             xlim = c(0, 1),
             comb.fixed =  FALSE,
             leftlabs = c("Study, Location", "AKI cases", "Sample size"),
             rightlabs = c("AKI incidence", "95% CI", "Weight"),
             pooled.events = TRUE,
             col.predict = "black")

dev.off()

# Meta-regression
aki_incidence_metareg <- metareg(aki_incidence_meta, ~ design + setting + location + age + males_p + aki_criteria + hypertension_p + all_cardio_p + diabetes_p + copd_p + ckd_p + cancer_p)

summary(aki_incidence_metareg)
