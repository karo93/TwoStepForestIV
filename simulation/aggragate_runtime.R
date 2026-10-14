# Read runtime files

extract_N_from_path <- function(path) {
  as.integer(
    stringr::str_extract(
      path,
      "(?<=N_)\\d+"
    )
  )
}

runtime_files <- list.files(
  here::here(
    "benchmark",
    "runtime"
  ),
  pattern = "^runtime_.*\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)


all_runtime <- purrr::map_dfr(
  runtime_files,
  function(f) {
    readRDS(f) %>%
      mutate(
        N = extract_N_from_path(f),
        outcome = "continuous"
      )
  }
)

all_runtime <- all_runtime %>%
  mutate(
    heterogeneity = case_when(
      str_detect(dataset, "strong") ~ "strong",
      str_detect(dataset, "slight") ~ "slight",
      TRUE ~ NA_character_
    ),
    k = case_when(
      str_detect(dataset, "_0_5\\.rds$") ~ 0.5,
      str_detect(dataset, "_4\\.rds$") ~ 4,
      TRUE ~ NA_real_
    )
  )


# Runtime comparison

runtime_plot <- all_runtime %>%
  filter(
    component == "total",
    method %in% c(
      "BCF_IV_tb",
      "DRRF_IV_tb",
      "GRF_IV_tb"
    ),
    k %in% c(0.5, 4)
  ) %>%
  mutate(
    
    method = dplyr::recode(
      method,
      "BCF_IV_tb" = "BCF-IV",
      "DRRF_IV_tb" = "DRRF-IV",
      "GRF_IV_tb" = "GRF-IV"
    ),
    
    k = factor(
      k,
      levels = c(0.5, 4),
      labels = c(
        "k = 0.5",
        "k = 4"
      )
    )
  ) %>%
  group_by(
    N,
    heterogeneity,
    k,
    method
  ) %>%
  summarise(
    mean_time = mean(
      time_sec,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  ggplot(
    aes(
      x = factor(N),
      y = mean_time,
      color = method,
      group = method
    )
  ) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(
    heterogeneity ~ k
  ) +
  scale_color_manual(
    values = c(
      "BCF-IV" = "#b2182b",
      "DRRF-IV" = "#92c5de",
      "GRF-IV" = "#0571b0"
    )
  ) +
  labs(
    x = "N",
    y = "Runtime (seconds)",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(
      fill = "#d1e5f0",
      color = "#d1e5f0"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    panel.border = element_blank(),
    panel.grid.major = element_line(
      color = "grey85"
    ),
    panel.grid.minor = element_line(
      color = "grey92"
    )
  )


ggsave(
  filename = here::here(
    "output",
    "benchmark_plots",
    "runtime_comparison.png"
  ),
  plot = runtime_plot,
  width = 8,
  height = 5,
  dpi = 600
)

# Runtime table

runtime_total_tbl <- all_runtime %>%
  filter(
    component == "total"
  ) %>%
  group_by(
    N,
    heterogeneity,
    k,
    method
  ) %>%
  summarise(
    mean_sec = mean(
      time_sec,
      na.rm = TRUE
    ),
    sd_sec = sd(
      time_sec,
      na.rm = TRUE
    ),
    median_sec = median(
      time_sec,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    N,
    heterogeneity,
    k,
    method
  )


runtime_total_tbl_clean <- runtime_total_tbl %>%
  filter(
    method != "BART_IV_cace_tb_slow"
  ) %>%
  mutate(
    method = dplyr::recode(
      method,
      "BCF_IV_tb" = "BCF-IV",
      "DRRF_IV_tb" = "DRRF-IV",
      "GRF_IV_tb" = "GRF-IV"
    )
  ) %>%
  arrange(
    N,
    heterogeneity,
    k,
    method
  )


print(
  runtime_total_tbl_clean
)
