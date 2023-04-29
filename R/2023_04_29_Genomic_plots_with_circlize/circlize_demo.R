library(tidyverse)
library(glue)
library(circlize)

# Magic numbers ----
intervalWidth <- 1e5

# Dataset links ----
variantSummaryPath <- "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"

genomicSuperDupsPath <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"

# Import ClinVar variant summary ----
variantSummary <- read_tsv(variantSummaryPath)

pathogenicSNVs <- variantSummary %>%
  filter(OriginSimple == "germline") %>%
  filter(Assembly == "GRCh38") %>%
  filter(ClinicalSignificance == "Pathogenic") %>%
  filter(Type == "single nucleotide variant") %>%
  filter(Chromosome != "MT") %>%
  mutate(chr = glue("chr{Chromosome}")) %>%
  group_by(chr) %>%
  mutate(posInterval = cut_width(Start, width = intervalWidth, boundary = 0)) %>%
  ungroup() %>%
  group_by(chr, posInterval) %>%
  summarise(nVariants = n()) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(newStart = as.integer(str_remove(
    str_split_1(as.character(posInterval), ",")[1], "\\[|\\("
  ))) %>%
  ungroup()

# Import segmental duplication dataset ----
genomicSuperDups <- read_tsv(
  genomicSuperDupsPath,
  col_names = c(
    "bin",
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "otherChr",
    "otherStart",
    "otherEnd",
    "otherSize",
    "uid",
    "posBasesHit",
    "testResult",
    "verdict",
    "chits",
    "ccov",
    "alignfile",
    "alignL",
    "indelN",
    "indelS",
    "alignB",
    "matchB",
    "mismatchB",
    "transitionsB",
    "transversionsB",
    "fracMatch",
    "fracMatchIndel",
    "jcK",
    "k2K"
  )
)

## Extract coordinates ----
genomicSuperDupsFiltered <- genomicSuperDups %>%
  filter(!str_detect(chrom, "random|chrUn")) %>%
  filter(!str_detect(otherChr, "random|chrUn"))

bed1 <- genomicSuperDups %>%
  select(chrom, chromStart, chromEnd) %>%
  setNames(c("chr", "start", "end"))

bed2 <- genomicSuperDups %>%
  select(otherChr, otherStart, otherEnd) %>%
  setNames(c("chr", "start", "end"))

# Create color function ----
minVariants <- log10(floor(min(pathogenicSNVs$nVariants)))

meanVariants <- log10(floor(mean(pathogenicSNVs$nVariants)))

maxVariants <- log10(ceiling(max(pathogenicSNVs$nVariants)))

colorFunction <-
  colorRamp2(c(minVariants, meanVariants, maxVariants),
             c("blue", "white", "red"))

# Make circlize plot ----
tiff(
  "circlize_demo.tiff",
  units = "cm",
  width = 17.35,
  height = 23.35,
  pointsize = 18,
  res = 300,
  compression = "lzw"
)
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg38")
circos.track(ylim = c(0, 100))
circos.trackLines(
  pathogenicSNVs$chr,
  x = pathogenicSNVs$newStart,
  y = rep(100, nrow(pathogenicSNVs)),
  type = "h",
  col = colorFunction(log10(pathogenicSNVs$nVariants))
)
circos.genomicLink(bed1, bed2, col = "coral", border = NA)
dev.off()
circos.clear()
