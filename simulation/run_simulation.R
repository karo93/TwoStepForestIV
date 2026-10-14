
Method <- "DRRF_IV"
heterogeneity <- "slight"
binary <- FALSE
dataset_names <- list.files(
  here::here(
    "simulation",
    "data",
    "N_2000"
  )
)
filter_datasets <- function(
    dataset_names,
    heterogeneity = c("slight", "strong"),
    binary = TRUE
) {
  
  heterogeneity <- match.arg(heterogeneity)
  
  hetero_pattern <- heterogeneity
  
  type_pattern <- if (binary) {
    "binary"
  } else {
    "continous"
  }
  
  filtered <- dataset_names[
    stringr::str_detect(
      dataset_names,
      hetero_pattern
    ) &
      stringr::str_detect(
        dataset_names,
        type_pattern
      )
  ]
  
  return(filtered)
}


dataset_names <- filter_datasets(
  dataset_names,
  heterogeneity = get("heterogeneity"),
  binary = get("binary")
)

for (j in 1:1) {
  dat1 <- readRDS(here::here("simulation","data","N_2000",dataset_names[j]))
  lapply(1:500,function(i) {
    print(i)
    data <- dat1 %>% filter(sim_id == i)
      x_preds <- DRR_IV
      X = data[, 6:15] %>% as.data.frame(), 
      Y = data$Y, 
      W = data$W, 
      Z = data$Z)
  
       results <- x_preds$results %>%
       extract_terminal_nodes() %>%
       mutate(true_CACE = sapply(node, function(node) {compute_subgroup_true_CACE(node, data) }),id = i)
      
      simulation_results <- list(
        results = results,
        discovery = x_preds$discovery,
        inference = x_preds$inference,
        pic = x_preds$pic
      )
      
      return(simulation_results)
    }
  ) %>%
    saveRDS(
      file = here::here(
        "simulation",
        "raw_results",
        Method,
        "N_2000",
        dataset_names[j]
      )
    )
}
