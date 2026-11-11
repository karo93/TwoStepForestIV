source(here::here("simulation", "helper_functions.R"))

Method <- "BCF_IV"
n_repetitions <- 100
binary <- TRUE

# Load application data
data <- readRDS(
  here::here(
    "application",
    "data",
    "icu18_complete.rds"
  )
)

# Run application
app_repeats <- lapply(
  seq_len(n_repetitions),
  function(i) {
    
    print(i)
    
    set.seed(i)
    
    X <- data[, 4:65] %>%
      as.data.frame()
    
    x_preds <- BCF_IV_app(
      X = X,
      Y = data$dead28,
      W = data$icu_bed,
      Z = data$open_bin,
      binary = binary
    )
    
    results <- x_preds$results %>%
      extract_terminal_nodes() %>%
      mutate(id = i)
    
    application_results <- list(
      results = results,
      discovery = x_preds$discovery,
      inference = x_preds$inference
    )
    
    return(application_results)
  }
)

# Save application results
output_dir <- here::here(
  "application",
  "raw_results",
  Method
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  app_repeats,
  file = file.path(
    output_dir,
    "app_repeats_bcf_true_site_included.rds"
  )
)
