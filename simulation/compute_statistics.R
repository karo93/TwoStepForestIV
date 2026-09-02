# Settings

choose_method <- "BCF_IV"
heterogeneity <- "slight"
binary <- FALSE
dataset_names <- list.files(
  here::here(
    "simulation",
    "raw_results",
    choose_method,
    "N_10000"
  )
)

lapply(dataset_names, 
       function(k){
         choose_dataset <- k
         dataset <- list.files(
           here::here(
             "simulation",
             "raw_results",
             choose_method,
             "N_10000",
             choose_dataset
           )
         )
         results_list <- map(dataset, function(sim) {
           sim$results %>%
             mutate(true_CACE = as.numeric(true_CACE))  # <- force conversion here
         })         # Bind all extracted results into a single dataframe
         processed_data <- bind_rows(results_list) %>%
           mutate(id = rep(1:length(dataset), times = map_int(dataset, ~ nrow(.x$results)))) %>%
           group_by(id) %>%
           summarise(
             mse_CCACE = mean((as.numeric(CCACE) - true_CACE)^2, na.rm = TRUE),
             mae_CCACE = mean(abs(as.numeric(CCACE) - true_CACE), na.rm = TRUE),
             bias_CCACE = mean(as.numeric(CCACE) - true_CACE, na.rm = TRUE),
             .groups = "drop"
           ) %>%
         saveRDS(file = here::here( "simulation","aggregated_results",choose_method,"N_10000","MSE",choose_dataset))
       })

library(stringr)

filter_datasets <- function(dataset_names, heterogeneity = c("slight", "strong"), binary = TRUE) {
  # Validate inputs
  heterogeneity <- match.arg(heterogeneity)
  
  # Choose type strings
  hetero_pattern <- heterogeneity
  type_pattern <- if (binary) "binary" else "continous"
  
  # Filter the datasets
  filtered <- dataset_names[
    stringr::str_detect(dataset_names, hetero_pattern) &
      stringr::str_detect(dataset_names, type_pattern)
  ]
  
  return(filtered)
}
heterogeneity <- "slight"
binary <- FALSE
dataset_names <-filter_datasets(dataset_names, heterogeneity = get("heterogeneity"), binary = get("binary"))


choose_method_list <- c(choose_method)

# Function to process a dataset
process_dataset <- function(dataset_name, heterogeneity = get("heterogeneity"), binary = get("binary")) {
  print(paste("Processing:", dataset_name))
  
  # Load dataset
  dataset <- readRDS(
    here::here(
      "simulation",
      "raw_results",
      choose_method,
      "N_10000",
      dataset_name
    )
  )
  # Extract inference and results lists
  inference_list <- purrr::map(dataset, "inference")
  results_list <- purrr::map(dataset, "results")
  processed_inference_list <- purrr::map(inference_list, function(inf_data) {
    if (choose_method == "BCF_IV") {
      # X is already inside the inference data frame (e.g., for BART)
      return(as.data.frame(inf_data))
      
    } else (choose_method %in% c("GRF_IV", "DRRF_IV")) {
      # X is stored separately in the list
      return(bind_cols(
        tibble(Y = inf_data$Y, W = inf_data$W, Z = inf_data$Z, pihat = inf_data$pihat),
        as.data.frame(inf_data$X)
      ))
      
    } 
  })
  # Process terminal conditions
  processed_terminal_conditions <- map(results_list, function(res) {
    res[["node"]] %>% as.character() %>% trimws() %>% clean_conditions()
  })
  
  results_per_id <- list()
  performance_results <- list()
  
  for (i in seq_along(processed_inference_list)) {
    print(paste("Processing ID:", i))
    results_tbl<-as.data.frame(results_list[[i]])
    inference_unlisted <- processed_inference_list[[i]]
    if (binary == TRUE) {
      terminal_conditions <- convert_to_binary(processed_terminal_conditions[i])
    } else {
      terminal_conditions <- processed_terminal_conditions[[i]]
    }
    inference_unlisted <- inference_unlisted %>%
      rowwise() %>%
      mutate(
        true_node = assign_true_node_all(
          as.list(cur_data()), 
          heterogeneity = get("heterogeneity"), 
          is_binary = get("binary")
        ),
        pred_info = list(assign_node_all(
          as.list(cur_data()), 
          terminal_conditions,  
          heterogeneity = get("heterogeneity"), 
          is_binary = get("binary")
        ))
      ) %>%
      ungroup() %>%
      tidyr::unnest_wider(pred_info) %>%   # expands into node_label + node_index
      rename(pred_node = pred_nodes) %>%
      mutate(
        Adj_pvalue = ifelse(!is.na(node_index),
                            results_tbl$Adj_pvalue[node_index],
                            NA_real_)
      )
    
  
    # Compute performance metrics
    performance_metrics <- compute_heterogeneity_metrics(
      inference_unlisted,
      heterogeneity,
      true_col = "true_node",
      pred_col = "pred_node",
      pval_col = "Adj_pvalue",
      alpha = 0.05
    )

    results_per_id[[i]] <- inference_unlisted
    performance_results[[i]] <- performance_metrics
  }
  
  # Save processed results
  saveRDS(
    results_per_id,
    file = here::here(
      "simulation",
      "aggregated_results",
      choose_method,
      "N_10000",
      "Results_per_id",
      dataset_name
    )
  )
  saveRDS(
    performance_results,
    file = here::here(
      "simulation",
      "aggregated_results",
      choose_method,
      "N_10000",
      "TPR",
      dataset_name
    )
  )
  
  return(performance_results)
}

# Process all datasets
all_performance_results <- lapply(dataset_names, function(k) {
  process_dataset(k,heterogeneity=get("heterogeneity") , binary=get("binary"))
})

# Main loop to process all split rules
for (current_rule in choose_method_list) {
  
  assign("choose_method", current_rule, envir = .GlobalEnv)  # dynamically set global var
  
  message("====================================")
  message("Processing split rule: ", current_rule)
  message("====================================")
  
  # Process each dataset for the current rule
  all_performance_results <- lapply(dataset_names, function(k) {
    process_dataset(k, heterogeneity = get("heterogeneity"), binary = get("binary"))
  })
  assign(paste0("results_", current_rule), all_performance_results, envir = .GlobalEnv)
}




