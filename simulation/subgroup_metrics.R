# important functions
generate_true_conditions <- function(heterogeneity = "slight", is_binary = FALSE) {
  X_vars <- paste0("X", 1:10)
  if (is_binary) {
    # Binary Case
    if (heterogeneity == "strong") {
      conditions <- list(
        paste0(X_vars[1], " == 0 & ", X_vars[2], " == 0"),  # True Node L1 (CACE = k)
        paste0(X_vars[1], " == 1 & ", X_vars[2], " == 1"),  # True Node L2 (CACE = -k)
        "TRUE"  # Catch-all case: Zero Otherwise (L3)
      )
    } else if (heterogeneity == "slight") {
      conditions <- list(
        paste0(X_vars[1], " == 0 & ", X_vars[2], " == 0"),  # True Node L1 (CACE = k)
        paste0(X_vars[1], " == 1 & ", X_vars[2], " == 1"),  # True Node L2 (CACE = -k)
        paste0(X_vars[1], " == 1 & ", X_vars[2], " == 0"),  # True Node L3 (CACE = 0.5k)
        paste0(X_vars[1], " == 0 & ", X_vars[2], " == 1")   # True Node L4 (CACE = -0.5k)
      )
    } else {
      stop("Invalid heterogeneity scenario. Choose 'strong' or 'slight'.")
    }
  } else {
    # Continuous Case
    if (heterogeneity == "strong") {
      conditions <- list(
        paste0(X_vars[1], " < -0.5 & ", X_vars[2], " < -0.5"),  # True Node L1 (CACE = k)
        paste0(X_vars[1], " > 0.5 & ", X_vars[2], " > 0.5"),    # True Node L2 (CACE = -k)
        "TRUE"  # Catch-all case: Zero Otherwise (L3)
      )
    } else if (heterogeneity == "slight") {
      conditions <- list(
        paste0(X_vars[1], " < 0 & ", X_vars[2], " < 0"),  # True Node L1 (CACE = k)
        paste0(X_vars[1], " > 0 & ", X_vars[2], " > 0"),    # True Node L2 (CACE = -k)
        paste0(X_vars[1], " > 0 & ", X_vars[2], " < 0"),   # True Node L3 (CACE = 0.5k)
        paste0(X_vars[1], " < 0 & ", X_vars[2], " > 0")    # True Node L4 (CACE = -0.5k)
      )
    } else {
      stop("Invalid heterogeneity scenario. Choose 'strong' or 'slight'.")
    }
  }
  
  return(conditions)
}
assign_true_node_all <- function(row, heterogeneity = "strong", is_binary = FALSE) {
  row_list <- as.list(row)  # Convert row to a named list
  true_conditions <- generate_true_conditions(heterogeneity, is_binary)
  for (i in seq_along(true_conditions)) {
    condition <- tryCatch(parse(text = true_conditions[i]), error = function(e) return(NULL))
    if (!is.null(condition)) {
      result <- tryCatch(eval(condition, envir = row_list), error = function(e) return(FALSE))
      
      if (result) {
       # Continuous Cases: Differentiate between strong/slight
        if (heterogeneity == "slight") {
          return(paste("L", i, sep = ""))  # Assign L1-L4
        } else if (heterogeneity == "strong") {
          if (i == 1) return("L1")  
          if (i == 2) return("L2")  
        }
      }
    }
  }
  #  }
  return("L0")
}
assign_node_all <- function(row, terminal_conditions, heterogeneity = "strong", is_binary = FALSE) {
  row_list <- as.list(row)
  node_label <- "L0"
  node_index <- NA_integer_
  
  for (j in seq_along(terminal_conditions)) {
    condition <- tryCatch(parse(text = terminal_conditions[j]), error = function(e) NULL)
    if (!is.null(condition)) {
      result <- tryCatch(eval(condition, envir = row_list), error = function(e) FALSE)
      if (isTRUE(result)) {
        # Decide node_label
        if (is_binary) {
          if (heterogeneity == "strong") {
            if (grepl("X1 == 0", terminal_conditions[j]) & grepl("X2 == 0", terminal_conditions[j])) {
              node_label <- "L1"
            } else if (grepl("X1 == 1", terminal_conditions[j]) & grepl("X2 == 1", terminal_conditions[j])) {
              node_label <- "L2"
            }
          } else if (heterogeneity == "slight") {
            if (grepl("X1 == 0", terminal_conditions[j]) & grepl("X2 == 0", terminal_conditions[j])) {
              node_label <- "L1"
            } else if (grepl("X1 == 1", terminal_conditions[j]) & grepl("X2 == 1", terminal_conditions[j])) {
              node_label <- "L2"
            } else if (grepl("X1 == 1", terminal_conditions[j]) & grepl("X2 == 0", terminal_conditions[j])) {
              node_label <- "L3"
            } else if (grepl("X1 == 0", terminal_conditions[j]) & grepl("X2 == 1", terminal_conditions[j])) {
              node_label <- "L4"
            }
          }
        } else {
          # Continuous case
          if (heterogeneity == "strong") {
            if (grepl("X1 <", terminal_conditions[j]) & grepl("X2 <", terminal_conditions[j])) {
              node_label <- "L1"
            } else if (grepl("X1 >", terminal_conditions[j]) & grepl("X2 >", terminal_conditions[j])) {
              node_label <- "L2"
            } else {
              node_label <- "L0"
            }
          } else if (heterogeneity == "slight") {
            if (grepl("X1 <", terminal_conditions[j]) & grepl("X2 <", terminal_conditions[j])) {
              node_label <- "L1"
            } else if (grepl("X1 >", terminal_conditions[j]) & grepl("X2 >", terminal_conditions[j])) {
              node_label <- "L2"
            } else if (grepl("X1 >", terminal_conditions[j]) & grepl("X2 <", terminal_conditions[j])) {
              node_label <- "L3"
            } else if (grepl("X1 <", terminal_conditions[j]) & grepl("X2 >", terminal_conditions[j])) {
              node_label <- "L4"
            }
          }
        }
        
        if (is.null(node_label)) node_label <- paste("Node", j)
        node_index <- j
        return(tibble::tibble(pred_nodes = node_label, node_index = node_index))
      }
    }
  }
  tibble::tibble(pred_nodes = "L0", node_index = NA_integer_)
}
convert_to_binary <- function(conditions) {
  binary_conditions <- sapply(conditions, function(cond) {
    cond <- gsub("([A-Za-z0-9_]+) >= 0\\.5", "\\1 == 1", cond)  # Fix regex spacing
    cond <- gsub("([A-Za-z0-9_]+) < 0\\.5", "\\1 == 0", cond)   # Ensure decimal point match
    return(cond)
  }, USE.NAMES = FALSE)  # Avoid named vector
  
  return(binary_conditions)
}
clean_conditions <- function(cond_list) {
  sapply(cond_list, function(x) {
    x <- gsub("([<>]=?)", " \\1 ", x, perl = TRUE)  # Adds space before and after <, >, <=, >=
    x <- gsub("\\s+", " ", x)  # Remove extra spaces
    return(trimws(x))  # Trim leading/trailing spaces
  }, USE.NAMES = FALSE)
}

normalize_node <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)                 # collapse whitespace
  x <- gsub("\\s*&\\s*", "&", x)            # remove spaces around &
  x <- gsub("\\s*(>=|<=|==|!=|>|<)\\s*", "\\1", x)  # remove spaces around operators
  x
}

compute_heterogeneity_metrics <- function(df,
                                          heterogeneity = c("strong", "slight"),
                                          true_col = "true_node",
                                          pred_col = "pred_node",
                                          pval_col = "Adj_pvalue",
                                          alpha = 0.05) {
  
  heterogeneity <- match.arg(heterogeneity)
  
  stopifnot(true_col %in% names(df), pred_col %in% names(df), pval_col %in% names(df))
  
  true_labels <- df[[true_col]]
  pred_labels <- df[[pred_col]]
  pvals       <- df[[pval_col]]
  
  if (heterogeneity == "strong") {
    classes <- c("L0", "L1", "L2")
  } else {
    classes <- sort(unique(true_labels))   # no L0 in slight
  }
  
  ## Detection metrics
  per_class <- lapply(classes, function(class) {
    TP <- sum(true_labels == class & pred_labels == class, na.rm = TRUE)
    FN <- sum(true_labels == class & pred_labels != class, na.rm = TRUE)
    FP <- sum(true_labels != class & pred_labels == class, na.rm = TRUE)
    TN <- sum(true_labels != class & pred_labels != class, na.rm = TRUE)
    
    tibble::tibble(
      Class = class,
      TP, FN, FP, TN,
      TrueDetectionRate  = TP / (TP + FN + 1e-9),
      FalseDetectionRate = FP / (FP + TN + 1e-9)
    )
  }) %>% dplyr::bind_rows()
  
  detection_summary <- tibble::tibble(
    macro_TDR = mean(per_class$TrueDetectionRate, na.rm = TRUE),
    macro_FDR = mean(per_class$FalseDetectionRate, na.rm = TRUE),
    micro_ODR = mean(true_labels == pred_labels, na.rm = TRUE)
  )
  
  # Significance metrics
  detected <- pvals < alpha
  
  sig_per_class <- lapply(classes, function(class) {
    in_class <- true_labels == class
    
    if (class == "L0") {
      TP_c <- sum(in_class & pred_labels == "L0" & !detected, na.rm = TRUE)
      FN_c <- sum(in_class & detected, na.rm = TRUE)
      FP_c <- sum(!in_class & pred_labels == "L0" & !detected, na.rm = TRUE)
      TN_c <- sum(!in_class & !(pred_labels == "L0" & !detected), na.rm = TRUE)
    } else {
      TP_c <- sum(in_class & pred_labels == class & detected, na.rm = TRUE)
      FN_c <- sum(in_class & !(pred_labels == class & detected), na.rm = TRUE)
      FP_c <- sum(!in_class & pred_labels == class & detected, na.rm = TRUE)
      TN_c <- sum(!in_class & !(pred_labels == class & detected), na.rm = TRUE)
    }
    
    precision_c <- TP_c / (TP_c + FP_c + 1e-9)
    tpr_c       <- TP_c / (TP_c + FN_c + 1e-9)
    f1_c        <- 2 * precision_c * tpr_c / (precision_c + tpr_c + 1e-9)
    fpr_c       <- FP_c / (FP_c + TN_c + 1e-9)
    sig_rate    <- mean(detected[in_class], na.rm = TRUE)
    
    tibble::tibble(
      Class = class,
      TP = TP_c, FN = FN_c, FP = FP_c, TN = TN_c,
      Precision = precision_c,
      TPR = tpr_c,
      F1 = f1_c,
      FPR = fpr_c,
      Sig_Rate = sig_rate
    )
  }) %>% dplyr::bind_rows()
  
  # in slight we do NOT create null_SigRate
  significance_summary <- tibble::tibble(
    macro_TPR = mean(sig_per_class$TPR, na.rm = TRUE),
    macro_Precision = mean(sig_per_class$Precision, na.rm = TRUE),
    macro_F1 = mean(sig_per_class$F1, na.rm = TRUE),
    macro_FPR = mean(sig_per_class$FPR, na.rm = TRUE),
    mean_signal_SigRate = mean(sig_per_class$Sig_Rate, na.rm = TRUE)
  )
  
  if (heterogeneity == "strong") {
    null_rate <- sig_per_class %>%
      dplyr::filter(Class == "L0") %>%
      dplyr::pull(Sig_Rate)
    
    significance_summary <- dplyr::mutate(
      significance_summary,
      null_SigRate = mean(null_rate, na.rm = TRUE),
      .before = mean_signal_SigRate
    )
  }
  return(list(
    heterogeneity = heterogeneity,
    detection_metrics = list(
      per_class = per_class,
      summary = detection_summary
    ),
    significance_metrics = list(
      per_class = sig_per_class,
      summary = significance_summary
    )
  ))
}


