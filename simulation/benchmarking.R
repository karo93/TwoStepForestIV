#Setting

adj_method <- "holm"

Rule <- "BART_IV_cace"

heterogeneity <- "strong"
binary <- FALSE

dataset_names <- list.files(
  here::here(
    "simulation",
    "data",
    "N_2000"
  )
)


# Select datasets

filter_datasets <- function(
    dataset_names,
    heterogeneity = c("slight", "strong"),
    binary = TRUE
) {
  
  heterogeneity <- match.arg(heterogeneity)
  
  hetero_pattern <- heterogeneity
  type_pattern <- if (binary) "binary" else "continous"
  
  filtered <- dataset_names[
    stringr::str_detect(dataset_names, hetero_pattern) &
      stringr::str_detect(dataset_names, type_pattern)
  ]
  
  return(filtered)
}


dataset_names <- filter_datasets(
  dataset_names,
  heterogeneity = get("heterogeneity"),
  binary = get("binary")
)

dataset_names_selected <- dataset_names[c(1, 2, 3, 4)]


# Select Monte Carlo replications for benchmarking

set.seed(1)

mc_ids <- sort(
  sample(
    1:500,
    100,
    replace = FALSE
  )
)
# Run

for (j in 1:1) {
  
  dat1 <- readRDS(
    here::here(
      "simulation",
      "data",
      "N_2000",
      dataset_names_selected[j]
    )
  )
 out <- lapply(
    mc_ids,
    function(i) {
      
      cat(
        "sim:",
        i,
        "\n"
      )
      
      data_i <- dat1 %>%
        filter(sim_id == i)
      
      
      x_preds <- BART_IV_cace_tb(
        X = as.data.frame(data_i[, 6:15]),
        Y = data_i$Y,
        W = data_i$W,
        Z = data_i$Z,
        sim_id = i
      )
      
      
      results <- x_preds$results %>%
        extract_terminal_nodes() %>%
        mutate(
          true_CACE = sapply(
            node,
            function(node) {
              compute_subgroup_true_CACE(
                node,
                data_i
              )
            }
          ),
          id = i
        )
      
      
      list(
        results = results,
        discovery = x_preds$discovery,
        inference = x_preds$inference,
        pic = x_preds$pic,
        timings = x_preds$timings
      )
    }
  )
  
 saveRDS(
    out,
    file = here::here(
      "benchmark",
      "raw_results",
      Rule,
      "N_2000",
      dataset_names_selected[j]
    )
  )
  
 runtime_tbl <- dplyr::bind_rows(
    lapply(
      seq_along(out),
      function(idx) {
        
        tm <- out[[idx]]$timings
        
        tibble(
          dataset = dataset_names_selected[j],
          sim_id = mc_ids[idx],
          method = "BCF_IV_tb",
          component = names(tm),
          time_sec = as.numeric(tm)
        )
      }
    )
  )
  
  
  saveRDS(
    runtime_tbl,
    file = here::here(
      "benchmark",
      "runtime",
      Rule,
      "N_2000",
      paste0(
        "runtime_",
        dataset_names_selected[j]
      )
    )
  )
}