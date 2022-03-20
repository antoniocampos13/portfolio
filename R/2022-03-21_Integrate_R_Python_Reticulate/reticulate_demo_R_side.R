# Install packages
# install.packages(c("here", "reticulate", "tidyverse", "openxlsx", "BiocManager"))
# BiocManager::install(c("BSgenome", "BSgenome.Hsapiens.UCSC.hg38"))

library(here)
library(tidyverse)
library(reticulate)
library(glue)

# Biocondutctor packages
library(GenomicRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

# Functions
source(here("src", "getSequence.R"))

# Download at NCBI's FTP site: 
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
gffPath <- here("GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff")

genes <- c("HTT", "FMR1")

use_condaenv("bioenv")

pythonScript <- here("reticulate_demo_Python_side.py")

source_python(pythonScript)

exonsCoordinates <- py$df 

exonsCoordinatesClean <- exonsCoordinates %>% 
  rowwise() %>%
  mutate(sequence = getSequence(Chromosome, Start, End)) %>%
  ungroup() %>%
  select(Chromosome, Feature, Start, End, Score, Strand, Frame, ID, gene, sequence)
