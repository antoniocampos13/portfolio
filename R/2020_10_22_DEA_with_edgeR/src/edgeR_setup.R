library(tidyverse)
library(openxlsx)

# Install trough BiocManager
library(edgeR)
library(annotate)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v79)
library(ensembldb)
library(AnnotationDbi)

#' @title edger_setup
#' @description Perform differential expression analysis using a dataset of expression quantification from a experiment with case/control design.
#' @param name A string. An identifier for the experiment.
#' @param counts The data frame containing the transcript counts with G rows and S samples. Rows must be named with gene/transcript ids.
#' @param replicates A Boolean indicating if the samples are biological replicates. Defaults to TRUE.
#' @param filter A Boolean indicating if lowly expressed transcripts should be filter out. Defaults to `TRUE`.
#' @param gene_id A string indicating how transcripts are identified in the data frame. There are three options: NCBI (Entrez gene ID), ENSEMBL (ENSEMBL ids ENS#) or SYMBOL (Official HGNC gene symbol).
#' @param output_path A path and filename string where the results will be saved.
#' @return A spreadsheet with Excel file extension (.xlsx). It contains data on gene symbol and id, log fold-change (lfc), log counts per million transcripts (logcpm), quasi-likelihood F-test statistic (f) and corresponding p-values, both raw (pvalue) and false-discovery rate-adjusted (adjpvalue).

edger_setup <- function(name, counts, replicates = TRUE, filter = TRUE, gene_id = c("NCBI", "ENSEMBL", "SYMBOL"), output_path) {
  temp <- counts %>% dplyr::select(ends_with("control"), everything())

  group <- factor(ifelse(str_detect(names(temp), "control"), 1, 2))

  edger_list <- DGEList(counts = temp, group = group)

  if (filter) {
    keep <- filterByExpr(edger_list)

    edger_list <- edger_list[keep, , keep.lib.sizes = FALSE]
  }

  edger_list <- calcNormFactors(edger_list)

  design <- model.matrix(~group)

  if (replicates) {
    edger_list <- estimateDisp(edger_list, design)

    fit <- glmQLFit(edger_list, design)

    qlf <- glmQLFTest(fit, coef = 2)


    if (gene_id == "NCBI") {
      gene_names <- as.vector(getSYMBOL(rownames(qlf$table), "org.Hs.eg"))
    }

    else if (gene_id == "ENSEMBL") {
      
      keys <- gsub("\\..*","",rownames(qlf$table))
      
      gene_names <- ensembldb::select(EnsDb.Hsapiens.v79, keys = keys, keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
      
    }

    else {
      gene_names <- rownames(qlf$table)
    }

    output <- cbind(gene_names, qlf$table)
    output$adjpvalue <- p.adjust(output$PValue, method = "fdr")
  }

  else {
    source(here("src", "hk_genes.R"))

    edger_list1 <- edger_list

    edger_list1$samples$group <- 1

    index <- row.names(edger_list1)

    edger_list0 <- estimateDisp(edger_list1[(index %in% housekeeping), ],
      trend = "none",
      tagwise = FALSE
    )

    edger_list$common.dispersion <- edger_list0$common.dispersion

    fit <- glmFit(edger_list, design)

    lrt <- glmLRT(fit)

    if (gene_id == "NCBI") {
      gene_names <- as.vector(getSYMBOL(rownames(lrt$table), "org.Hs.eg"))
    }

    else if (gene_id == "ENSEMBL") {
      
      keys <- gsub("\\..*","",rownames(qlf$table))
      
      gene_names <- ensembldb::select(EnsDb.Hsapiens.v79, keys = keys, keytype = "GENEID", columns = c("SYMBOL", "GENEID"))
    }

    else {
      gene_names <- rownames(lrt$table)
    }

    output <- cbind(gene_names, lrt$table)
    output$adjpvalue <- p.adjust(output$PValue, method = "fdr")
  }

  output <- output %>% rename_all(tolower)

  if (dim(output)[1] > 0) {
    wb <- createWorkbook()

    addWorksheet(wb, sheetName = name, gridLines = TRUE)

    writeData(wb, sheet = name, x = output)

    saveWorkbook(wb, output_path, overwrite = FALSE)

    # write.xlsx(output, file = output_path, sheetName = name, append = TRUE, row.names = FALSE)

    paste0("DEGs detected. Saving to spreadsheet on ", output_path)
  }

  else {
    paste0("No DEGs detected.")
  }
}
