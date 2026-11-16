library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(here)
library(stringr)

# Simulation settings

N <- 2000
methods <- c("BCF_IV", "CF_IV", "DRRF_IV")

accuracy_results <- list()

for (rule in methods) {
  message("Processing split rule: ", rule)
  
  dataset_names <- list.files(here::here("simulation", "aggregated_results", rule, paste0("N_", N), "TPR"))  
  
  for (dataset in dataset_names) {
    result_path <- here::here("simulation", "aggregated_results", rule, paste0("N_", N), "TPR", dataset)    
    if (!file.exists(result_path)) {
      message("  Missing file: ", dataset)
      next
    }
    
    results <- readRDS(result_path)
    
    if (length(results) == 0 || !is.list(results)) {
      message("  Invalid or empty results: ", dataset)
      next
    }
    
    hetero <- unique(map_chr(results, "heterogeneity"))[1]
    
    # Holders
    det_matrix <- list()
    sig_matrix <- list()
    
    for (res in results) {
      
      #Detection metrics: per class ---
      if (!is.null(res$detection_metrics$per_class)) {
        det_metrics <- res$detection_metrics$per_class %>%
          select(Class, TrueDetectionRate, FalseDetectionRate)
        
        tdr_row <- det_metrics %>%
          select(Class, TrueDetectionRate) %>%
          pivot_wider(
            names_from = Class,
            values_from = TrueDetectionRate,
            names_prefix = "TDR_"
          )
        
        fdr_row <- det_metrics %>%
          select(Class, FalseDetectionRate) %>%
          pivot_wider(
            names_from = Class,
            values_from = FalseDetectionRate,
            names_prefix = "FDR_"
          )
        
        det_matrix[[length(det_matrix) + 1]] <- bind_cols(tdr_row, fdr_row)
      }
      
      ## --- Significance metrics: per class ---
      if (!is.null(res$significance_metrics$per_class)) {
        sig_metrics <- res$significance_metrics$per_class %>%
          select(Class, Precision, TPR, F1, FPR, Sig_Rate)
        
        sig_row <- sig_metrics %>%
          pivot_wider(
            names_from = Class,
            values_from = c(Precision, TPR, F1, FPR, Sig_Rate),
            names_sep = "_"
          )
        
        sig_matrix[[length(sig_matrix) + 1]] <- sig_row
      }
    }
    
    # class-specific detection metrics
    det_means <- if (length(det_matrix) > 0) {
      bind_rows(det_matrix) %>%
        summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
    } else {
      tibble()
    }
    
    # class-specific significance metrics
    sig_means <- if (length(sig_matrix) > 0) {
      bind_rows(sig_matrix) %>%
        summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
    } else {
      tibble()
    }
    
    #detection summary: use stored summary if available
    det_summary <- tibble()
    
    if (length(results) > 0) {
      det_summary_raw <- map(results, ~ .x$detection_metrics$summary) %>%
        compact() %>%
        bind_rows()
      
      if (nrow(det_summary_raw) > 0) {
        det_summary <- det_summary_raw %>%
          summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
      } else {
        det_all <- map(results, ~ .x$detection_metrics$per_class) %>%
          compact() %>%
          bind_rows(.id = "replication")
        
        if (nrow(det_all) > 0) {
          det_summary <- tibble(
            macro_TDR = mean(det_all$TrueDetectionRate, na.rm = TRUE),
            macro_FDR = mean(det_all$FalseDetectionRate, na.rm = TRUE),
            macro_TDR = mean(
              map_dbl(results, ~ .x$detection_metrics$summary$macro_TDR %||% NA_real_),
              na.rm = TRUE
            )
          )
        }
      }
    }
    
    #significance summary: ALWAYS recompute from per_class
    #this avoids inheriting an old macro_TDR that incorrectly includes L0
    sig_summary <- tibble()
    
    if (length(results) > 0) {
      sig_all <- map(results, ~ .x$significance_metrics$per_class) %>%
        compact() %>%
        bind_rows(.id = "replication")
      
      if (nrow(sig_all) > 0) {
        sig_all <- sig_all %>%
          mutate(Class = as.character(Class))
        
        base_summary <- tibble(
          macro_TPR       = mean(sig_all$TPR, na.rm = TRUE),
          macro_Precision = mean(sig_all$Precision, na.rm = TRUE),
          macro_F1        = mean(sig_all$F1, na.rm = TRUE),
          macro_FPR       = mean(sig_all$FPR, na.rm = TRUE)
        )
        
        if (hetero == "slight") {
          sig_summary <- base_summary %>%
            mutate(
              null_SigRate = NA_real_,
              macro_TDR = mean(sig_all$Sig_Rate[sig_all$Class != "L0"], na.rm = TRUE)
            )
        } else if (hetero == "strong") {
          sig_summary <- base_summary %>%
            mutate(
              null_SigRate = mean(sig_all$Sig_Rate[sig_all$Class == "L0"], na.rm = TRUE),
              macro_TDR = mean(sig_all$Sig_Rate[sig_all$Class != "L0"], na.rm = TRUE)
            )
        } else {
          sig_summary <- base_summary %>%
            mutate(
              null_SigRate = NA_real_,
              macro_TDR = mean(sig_all$Sig_Rate, na.rm = TRUE)
            )
        }
      }
    }
    
    ## --- Parse k from dataset name ---
    k_val <- dplyr::case_when(
      grepl("_0_5", dataset) ~ 0.5,
      grepl("_1_5", dataset) ~ 1.5,
      grepl("_1(_|\\.)", dataset) ~ 1.0,
      grepl("_2_5", dataset) ~ 2.5,
      grepl("_2(_|\\.)", dataset) ~ 2.0,
      grepl("_3_5", dataset) ~ 3.5,
      grepl("_3(_|\\.)", dataset) ~ 3.0,
      grepl("_4(_|\\.)", dataset) ~ 4.0,
      TRUE ~ NA_real_
    )
    
    # Store combined result
    accuracy_results[[length(accuracy_results) + 1]] <- tibble(
      dataset = dataset,
      split_rule = rule,
      heterogeneity = hetero,
      k = k_val
    ) %>%
      bind_cols(det_means, sig_means, det_summary, sig_summary)
  }
}

# Final tidy tibble
accuracy_df <- bind_rows(accuracy_results) %>%
  distinct()
  saveRDS(accuracy_df,
    here::here("output", "sim_tables", paste0("accuracy_df_N", N, ".rds")))

extract_k_vec <- function(filenames) {
  # Match the first k pattern: 0_2, 1, 4, etc.
  matches <- stringr::str_extract(filenames, "\\d+_\\d+|\\d+")
  ifelse(
    is.na(matches),
    NA_real_,
    as.numeric(gsub("_", ".", matches))
  )
}

accuracy_df<- accuracy_df%>%
  mutate(k = extract_k_vec(dataset))


extract_dgp_group <- function(filename) {
  # Remove the trailing _k.rds pattern (e.g., _0_2.rds, _1.rds, _4.rds)
  stringr::str_remove(filename, "_\\d+_\\d+\\.rds$|_\\d+\\.rds$")
}

accuracy_df <- accuracy_df %>%
  mutate(
    dgp_group = extract_dgp_group(dataset),
    split_rule = case_when(
      split_rule == "DRRF_IV" ~ "DRRF-IV",
      split_rule == "BCF_IV" ~ "BCF-IV",
      split_rule == "CF_IV" ~ "GRF-IV"
    )
  )

accuracy_df_clean <- accuracy_df %>%
  select(-dataset, -dgp_group)
#######################################

accuracy_df_strong <- accuracy_df %>%
  filter(str_detect(dataset, "strong")) %>%
  select(where(~ !all(is.na(.x))))

saveRDS(
  accuracy_df_strong,
  here::here("output", "sim_tables", paste0("accuracy_df_strong_N", N, ".rds"))
)

accuracy_df_slight <- accuracy_df %>%
  filter(str_detect(dataset, "slight")) %>%
  select(where(~ !all(is.na(.x))))

saveRDS(
  accuracy_df_slight,
  here::here("output", "sim_tables", paste0("accuracy_df_slight_N", N, ".rds"))
)

# graphs 


plot_df <- accuracy_df %>%
  filter(split_rule %in% c("BCF-IV", "DRRF-IV", "GRF-IV")) %>%
  select(split_rule, heterogeneity, k, macro_FPR, macro_FDR) %>%
  pivot_longer(
    cols = c(macro_FPR, macro_FDR),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    heterogeneity = factor(
      heterogeneity,
      levels = c("strong", "slight"),
      labels = c("Strong Heterogeneity", "Slight Heterogeneity")
    ),
    metric = factor(
      metric,
      levels = c("macro_FPR", "macro_FDR"),
      labels = c("False Positive Rate", "False Detection Rate")
    ),
    split_rule = factor(
      split_rule,
      levels = c("BCF-IV", "DRRF-IV", "GRF-IV")
    )
  )

ggplot(plot_df, aes(x = k, y = value, color = split_rule, group = split_rule)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.2) +
  facet_grid(heterogeneity ~ metric, scales = "free_y") +
  scale_color_manual(values = c(
    "BCF-IV" = "#b2182b",
    "DRRF-IV" = "#92c5de",
    "GRF-IV" = "#0571b0"
  )) +
  scale_x_continuous(breaks = c(1, 1.5, 2, 2.5, 3, 3.5, 4)) +
  labs(
    x = "Effect Size (k)",
    y = NULL,
    color = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#d1e5f0", color = "#d1e5f0"),
    strip.text = element_text(face = "bold"),
    panel.border = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_line(color = "grey92"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )
ggsave(filename = here::here("output", "sim_plots", paste0("FPR_FDR_N", N, ".png")),width = 9,height = 5.5, dpi = 300)

#table
accuracy_table <- accuracy_df %>%
  transmute(
    method = split_rule,
    heterogeneity,
    k,
    TPR = round(macro_TPR, 3),
    TDR = round(macro_TDR, 3),
    FPR = round(macro_FPR, 3),
    FDR = round(macro_FDR, 3),
    ODR = round(micro_ODR, 3),
    SigRate = round(macro_SigRate, 3),
    null_SigRate = round(null_SigRate, 3)
  ) %>%
  arrange(heterogeneity, method, k)

accuracy_df <- bind_rows(accuracy_results) %>% distinct()

  saveRDS(
  accuracy_table,
  here::here("output", "sim_tables", paste0("accuracy_table_N", N, ".rds"))
)
################## 

mse_results <- list()

for (rule in methods) {
  message("Processing split rule: ", rule)
  
  dataset_names <- list.files(here::here("simulation", "aggregated_results", rule, paste0("N_", N), "MSE"))  
  for (dataset in dataset_names) {
    result_path <- here::here( "simulation", "aggregated_results", rule, paste0("N_", N), "MSE", dataset)    
    if (!file.exists(result_path)) {
      message("  Missing file: ", dataset)
      next
    }
    
    results <- readRDS(result_path)
    
    if (!is.data.frame(results) || nrow(results) == 0) {
      message("  Invalid or empty results: ", dataset)
      next
    }
    
    hetero <- case_when(
      str_detect(dataset, "strong") ~ "strong",
      str_detect(dataset, "slight") ~ "slight",
      TRUE ~ NA_character_
    )
    
    k_val <- case_when(
      grepl("_0_5", dataset) ~ 0.5,
      grepl("_1_5", dataset) ~ 1.5,
      grepl("_1(_|\\.)", dataset) ~ 1.0,
      grepl("_2_5", dataset) ~ 2.5,
      grepl("_2(_|\\.)", dataset) ~ 2.0,
      grepl("_3_5", dataset) ~ 3.5,
      grepl("_3(_|\\.)", dataset) ~ 3.0,
      grepl("_4(_|\\.)", dataset) ~ 4.0,
      TRUE ~ NA_real_
    )
    
    mse_summary <- results %>%
      summarise(
        mean_MSE = mean(mse_CCACE, na.rm = TRUE),
        sd_MSE = sd(mse_CCACE, na.rm = TRUE),
        median_MSE = median(mse_CCACE, na.rm = TRUE),
        mean_Bias = mean(bias_CCACE, na.rm = TRUE),
        mean_MAE = mean(mae_CCACE, na.rm = TRUE),
        sd_Bias = sd(bias_CCACE, na.rm = TRUE),
        median_Bias = median(bias_CCACE, na.rm = TRUE)
      )
    
    mse_results[[length(mse_results) + 1]] <- tibble(
      dataset = dataset,
      split_rule = rule,
      heterogeneity = hetero,
      k = k_val
    ) %>%
      bind_cols(mse_summary)
  }
}


mse_df <- bind_rows(mse_results)  %>%
  distinct()
mse_table <- mse_df %>%
  transmute(
    Method = split_rule,
    heterogeneity,
    k,
    MSE = round(mean_MSE, 3),
    Bias = round(mean_Bias, 3),
    MAE = round(mean_MAE, 3)
  ) %>%
  arrange(heterogeneity, Method, k)

mse_table
saveRDS(  mse_table, here::here("output", "sim_tables", paste0("mse_table_N", N, ".rds")))