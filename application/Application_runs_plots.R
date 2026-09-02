source(file = here::here("simulation", "packages.R"))

Rule <- "BCF_IV"
data <- readRDS(
  here::here(
    "application",
    "data",
    "icu18_complete.rds"
  )
)
# Repeated subgroup discovery, for single run, change 100 to 1
app_repeats <- lapply(
  1:100,
  function(i) {
   X <- data[, 4:65] # here we included all site variables, see readme
   set.seed(i)
    
    x_preds <- BCF_IV(
      X = X %>% as.data.frame(),
      Y = data$dead28,
      W = data$icu_bed,
      Z = data$open_bin,
      binary = TRUE
    )
    
    results <- x_preds$results %>%
      extract_terminal_nodes() %>%
      mutate(id = i)
    
    simulation_results <- list(
      results = results,
      discovery = x_preds$discovery,
      inference = x_preds$inference
    )
    return(simulation_results)
  }
)


save_path <- here::here(
  "application",
  "raw_results",
  Rule,
  "repeats_bcf_iv.rds"
)

saveRDS(
  app_repeats,
  file = save_path
)

# Variable selection frequency

xtract_vars <- function(res) {
  
  vars <- stringr::str_extract_all(
    res$node,
    "[a-zA-Z0-9_]+(?=\\s*[<>=])"
  )
  
  unlist(vars)
}

# Collect variables from all repeated runs
all_vars <- unlist(
  lapply(
    app_repeats_grf,
    function(run) {
      extract_vars(run$results)
    }
  )
)

var_freq <- as.data.frame(
  table(all_vars)
) %>%
  arrange(
    desc(Freq)
  )

# Plot
png(
  here::here(
    "output",
    "application_plots",
    "var_impo_bcf.png"
  ),
  width = 1200,
  height = 800,
  res = 150
)

ggplot(
  var_freq,
  aes(
    x = reorder(all_vars, Freq),
    y = Freq
  )
) +
  geom_bar(
    stat = "identity",
    fill = "steelblue"
  ) +
  coord_flip() +
  labs(
    x = "Variable",
    y = "Frequency across runs (leaves)"
  )

dev.off()


# Plot single subgroup tree, run the code inside bcf_iv function to have everything in your environment
# Make a copy so the original fitted tree is unchanged
fit.tree.display <- fit.tree


# Replace fitted node values with CCACE estimates from the inference sample
ccace_vec <- as.numeric(
  ivResults$CCACE[
    match(
      rownames(fit.tree.display$frame),
      rownames(ivResults)
    )
  ]
)

fit.tree.display$frame$yval <- ccace_vec


# Construct node labels
all_node_ids <- as.numeric(
  rownames(fit.tree$frame)
)

all_labels <- sapply(
  all_node_ids,
  function(node_id) {
    
    row_idx <- which(
      rownames(ivResults) == node_id
    )
    
    if (length(row_idx) == 0) {
      return("")
    }
    
    ccace <- as.numeric(
      ivResults$CCACE[row_idx]
    )
    
    pi_obs <- as.numeric(
      ivResults$Pi_obs[row_idx]
    )
    
    paste0(
      "CCACE = ",
      round(ccace, 2),
      "\n",
      round(pi_obs * 100, 1),
      "%"
    )
  }
)


png(
  here::here(
    "output",
    "application_plots",
    "grf_iv_single_tree.png"
  ),
  width = 1200,
  height = 800,
  res = 150
)

rpart.plot::rpart.plot(
  fit.tree.display,
  extra = 0,
  under = TRUE,
  fallen.leaves = TRUE,
  node.fun = function(x, labs, digits, varlen) {
    
    node_ids <- as.numeric(
      rownames(x$frame)
    )
    
    labs <- all_labels[
      match(
        node_ids,
        all_node_ids
      )
    ]
    
    labs
  }
)

dev.off()