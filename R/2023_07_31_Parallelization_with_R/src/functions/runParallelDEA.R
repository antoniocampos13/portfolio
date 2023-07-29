library(glue)
library(furrr)
library(future.callr)
library(tictoc)

runParallelDEA <- function(nWorkers) {
  
  plan(callr, gc = TRUE, workers = nWorkers)
  
  opts <- furrr_options(scheduling = TRUE, seed = TRUE)
  
  parallelFileNames <- glue("{names(countsDfsList)}_parallel_{nWorkers}_workers.xlsx")
  
  logMessage <- ifelse(nWorkers == 1, "Sequential execution", glue("Parallel execution with {nWorkers} workers"))
  
  tic(logMessage)
  future_walk2(
    .x = parallelFileNames,
    .y = countsDfsList,
    ~ edger_setup(
      name = as.character(which(parallelFileNames == .x)),
      counts = .y,
      gene_id = "ENSEMBL",
      output_path = here("outputs", .x)
    ),
    .options = opts
  )
  toc(log = TRUE)
  parallelTime <- unlist(tic.log(format = TRUE))
  tic.clearlog()

  parallelTime %>% write_lines(here("outputs", "tictoc_log.txt"), append = file.exists(here("outputs", "tictoc_log.txt")))

}
