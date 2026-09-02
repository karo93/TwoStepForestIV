adj_method="holm"
extract_terminal_nodes <- function(results) {
  results <- results %>% filter(!is.na(node))
  all_nodes <- results$node
  # Identify parent nodes
  parent_nodes <- unique(unlist(lapply(all_nodes, function(x) {
    node_parts <- strsplit(x, " & ")[[1]]
    parent_nodes <- sapply(seq_along(node_parts) - 1, function(i) paste(node_parts[1:i], collapse = " & "))
    return(parent_nodes)
  })))
  terminal_nodes <- results %>% filter(!(node %in% parent_nodes))
  return(terminal_nodes)
}
compute_subgroup_true_CACE <- function(node_condition, data) {
  subset_data <- tryCatch(
    {
      data %>% filter(eval(parse(text = node_condition)))
    },
    error = function(e) {
      message("Error processing node: ", node_condition)
      return(data.frame())
    }
  )
  if (nrow(subset_data) == 0) {
    return(NA)
  }
  return(mean(subset_data$CACE, na.rm = TRUE))
}
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}
DRRF_IV <- function(Y, W, Z, X) {
  #Propensity Scores
  p.score <- glm(Z ~ ., family = binomial, data = X)
  pihat <- predict(p.score, type = "response")
  # Identify indices for which propensity scores are within [0.1, 0.9]
  valid_idx <- which(pihat > 0.1 & pihat < 0.9)
  subset_data <- list(
    Y = Y[valid_idx],
    W = W[valid_idx],
    Z = Z[valid_idx],
    X = X[valid_idx, ],
    pihat = pihat[valid_idx]
  )
  set.seed(42)  
  # Number of valid observations
  n_valid <- length(subset_data$Y)
  discovery_id <- sample(n_valid, size = floor(0.5 * n_valid))  # 50% for training
  inference_id <- setdiff(1:n_valid, discovery_id)  #  50% for validation
  # Create train/validation splits
  discovery <- list(
    Y = subset_data$Y[discovery_id],
    W = subset_data$W[discovery_id],
    Z = subset_data$Z[discovery_id],
    X = subset_data$X[discovery_id, ],
    pihat = subset_data$pihat[discovery_id]
  )
  inference <- list(
    Y = subset_data$Y[inference_id],
    W = subset_data$W[inference_id],
    Z = subset_data$Z[inference_id],
    X = subset_data$X[inference_id, ],
    pihat = subset_data$pihat[inference_id]
  )
  inference_combined <- data.frame(
    Y = inference$Y,
    W = inference$W,
    Z = inference$Z,
    inference$X,  # This expands the columns of X
    pihat = inference$pihat
  )
  
  m_model <- randomForest::randomForest(x =as.data.frame(discovery$X), y = discovery$Y)
  m_hat <- predict(m_model)
  p.score_itt <- glm(discovery$Z ~ ., family = binomial, data = discovery$X)
  pihat_itt <- predict(p.score_itt, type = "response")
  transformed_Y <- (discovery$Z - pihat_itt) * (discovery$Y - m_hat) / (pihat_itt * (1 - pihat_itt))
  # ITT and PIC estimation
  itt_model_rf<-grfRemote::regression_forest(Y=transformed_Y,X=as.data.frame(discovery$X),num.trees=2000)
  itt<- itt_model_rf$predictions
  pic_model<- grf::causal_forest(Y=discovery$W , W=discovery$Z, X=discovery$X  ,num.trees=2000)
  
  pic <- pic_model$predictions
  epsilon <- 1e-6
  tauhat <- itt / (pic+epsilon)

  exp_data <- as.data.frame(cbind(tauhat, as.data.frame(discovery$X)))
  
  #Build a CART with ITT
  fit.tree <- rpart(tauhat ~ .,
                    data = exp_data,
                    cp=0.001)
  
  
  rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
  
  # Initialize Outputs
 ivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol = 10))
  names(ivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test",
                    "Pi_obs", "ITT", "Pi_compliers",
                    "SE", "CI_lower", "CI_upper")
  lvs <- leaves <- numeric(length(rules)) 
  lvs[unique(fit.tree$where)] <- 1
  leaves[rules[lvs==1]] <- 1
  
  #Run an IV Regression 
  iv.root <- AER::ivreg(Y ~ W | Z,  
                        data = inference)
  summary <- summary(iv.root, diagnostics = TRUE)
  iv.effect.root <-  summary$coef[2,1]
  p.value.root <- summary$coef[2,4]
  p.value.weak.iv.root <- summary$diagnostics[1,4]
  proportion.root <- 1
  compliers.root <- length(which(inference$Z==inference$W))/nrow(inference$X)
  itt.root <- iv.effect.root*compliers.root
  
  se_ccace.root <- summary$coef[2,2]   # Standard error at root
  
  # Compute 95% Confidence Interval
  ci_lower_root <- iv.effect.root - 1.96 * se_ccace.root
  ci_upper_root <- iv.effect.root + 1.96 * se_ccace.root
  
  ivMat[1,] <- c(NA, round(iv.effect.root, 4),
                 round(p.value.root, 4),
                 round(p.value.weak.iv.root, 4),
                 round(proportion.root, 4),
                 round(itt.root, 4),
                 round(compliers.root, 4),
                 round(se_ccace.root, 4),
                 round(ci_lower_root, 4),
                 round(ci_upper_root, 4))
  names(inference) <- paste(names(inference), sep="")
  
  # Run a loop to get the rules (sub-populations)
  for (i in rules[-1]) {
    sub <- as.data.frame(matrix(NA, nrow = 1,
                                ncol = nrow(as.data.frame(path.rpart(fit.tree, node = i, print.it = FALSE))) - 1))
    quiet(capture.output(for (j in 1:ncol(sub)) {
      sub[, j] <- as.character(print(as.data.frame(path.rpart(fit.tree, node = i, print.it = FALSE))[j + 1, 1]))
      sub_pop <- noquote(paste(sub, collapse = " & "))
    }))
    
    subset <- with(inference_combined, inference_combined[which(eval(parse(text = sub_pop))), ])
    
    if (nrow(subset) > 3 &&
        length(unique(subset$W)) > 1 &&
        length(unique(subset$Z)) > 1 &&
        all(!is.na(subset$Y)) && all(!is.na(subset$W)) && all(!is.na(subset$Z))) {
      
      X_model <- tryCatch(model.matrix(~ W | Z, data = subset), error = function(e) return(NULL))
      if (is.null(X_model) || qr(X_model)$rank < ncol(X_model)) {
        message("Skipping node: ", sub_pop, " due to rank-deficient model matrix.")
        next
      }
      
      iv.reg <- tryCatch({
        AER::ivreg(Y ~ W | Z, data = subset)
      }, error = function(e) {
        message("IV regression failed for node: ", sub_pop, " - Skipping...")
        return(NULL)
      })
      
      if (!is.null(iv.reg)) {
        summary <- tryCatch({
          summary(iv.reg, diagnostics = TRUE)
        }, error = function(e) {
          message("Summary failed for node: ", sub_pop, " - Skipping...")
          return(NULL)
        })
        
        if (!is.null(summary) && summary$df[2] > 0) {
          iv.effect <- summary$coef[2, 1]
          se_ccace <- summary$coef[2, 2]
          ci_lower <- iv.effect - 1.96 * se_ccace
          ci_upper <- iv.effect + 1.96 * se_ccace
          
          p.value <- summary$coef[2, 4]
          p.value.weak.iv <- summary$diagnostics[1, 4]
          compliers <- length(which(subset$Z == subset$W)) / nrow(subset)
          itt <- iv.effect * compliers
          proportion.node <- nrow(subset) / nrow(inference$X)
          
          ivMat[i, ] <- c(sub_pop,
                          round(iv.effect, 4),
                          round(p.value, 4),
                          round(p.value.weak.iv, 4),
                          round(proportion.node, 4),
                          round(itt, 4),
                          round(compliers, 4),
                          round(se_ccace, 4),
                          round(ci_lower, 4),
                          round(ci_upper, 4))
        } else {
          message("Skipping node: ", sub_pop, " due to zero residual degrees of freedom.")
        }
      }
    } else {
      message("Skipping node: ", sub_pop, " due to missing values, small sample size, or lack of variation.")
    }
  }    
  # Adjust P.values 
  iv_correction <- cbind(as.data.frame(ivMat), leaves)
  adj <- round(p.adjust( as.numeric(iv_correction$pvalue[which(iv_correction$leaves==1)]) ,  paste(adj_method)), 5)
  Adj_pvalue <- rep(NA, length(rules)) 
  Adj_pvalue[which(iv_correction$leaves==1)] <- adj
  
  # Store Results
  ivResults <- cbind(as.data.frame(ivMat), Adj_pvalue)
  return(list(results = ivResults, discovery = discovery, inference = inference, pic=pic))
  
}
BCF_IV <- function(Y, W, Z, X,  binary = FALSE) {
  ######################################################
  ####         Step 1: Propensity Scores            ####
  ######################################################
  set.seed(42)  # For reproducibility
  index <- sample(nrow(X), nrow(X)*0.5, replace=FALSE) #inference Ratio = 0.5
  
  # Initialize total dataset
  iv.data <- as.data.frame(cbind(Y, W, Z, X))
  names(iv.data)[1:3] <- c("Y", "W", "Z")
  iv.data <- iv.data %>%
    mutate(across(c(W, Z, starts_with("X")), as.numeric))
  # Estimate propensity scores
  x_vars <- setdiff(names(iv.data), c("Y", "W", "Z"))
  p.score <- glm(
    Z ~ .,
    family = binomial,
    data = data.frame(Z = iv.data$Z, iv.data[, x_vars])
  )
  pihat_all <- predict(p.score, type = "response")
  valid_idx <- which( pihat_all > 0.1 & pihat_all < 0.9)
  iv.data <- iv.data[valid_idx, ]
  pihat_all <- pihat_all[valid_idx]
  index <- sample(nrow(iv.data), nrow(iv.data) * 0.5, replace = FALSE)
 
  # Discovery and Inference Samples
  discovery <- iv.data[-index,]
  inference <- iv.data[index,]
  discovery_X<-discovery[,4:13]
 
  # Perform the Bayesian Causal Forest for the Proportion of Compliers (pic)
  pic_bcf <- quiet(bartCause::bartc(discovery$W, discovery$Z, discovery_X, n.samples = 500, n.burn = 500, n.chains = 2L))
  tau_pic <- bartCause::extract(pic_bcf, type = "ite")
  pic <- apply(tau_pic, 2, mean)
  if (binary == FALSE){
    
    bcf_fit <- SparseBCF::SparseBCF(
      y = discovery$Y,
      z = as.numeric(discovery$Z),
      x_control = as.matrix(discovery_X),
      x_moderate = as.matrix(discovery_X),
      pihat = pihat,
      sparse = FALSE,  # <<< THIS IS THE KEY
      nburn = 500,
      nsim = 500,
      save_trees_mu_dir = NULL,
      save_trees_tau_dir = NULL,
      OOB = FALSE
    )
    
    # Compute Posterior
    tau_post <- bcf_fit$tau
    tauhat <- colMeans(tau_post)
    exp <- as.data.frame(cbind(tauhat, discovery_X))
    fit.tree <- rpart(tauhat ~ .,
                      data = exp,
                      cp=0.001,
                      maxdepth = 2)
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol = 10))
    names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test",
                         "Pi_obs", "ITT", "Pi_compliers",
                         "SE", "CI_lower", "CI_upper")
    # Generate Leaves Indicator
    lvs <- leaves <- numeric(length(rules)) 
    lvs[unique(fit.tree$where)] <- 1
    leaves[rules[lvs==1]] <- 1
    
    #Run an IV Regression on each Node
    
    
    # Run an IV Regression on the Root
    iv.root <- AER::ivreg(Y ~ W | Z,  
                          data = inference)
    summary <- summary(iv.root, diagnostics = TRUE)
    iv.effect.root <-  summary$coef[2,1]
    p.value.root <- summary$coef[2,4]
    p.value.weak.iv.root <- summary$diagnostics[1,4]
    proportion.root <- 1
    compliers.root <- length(which(inference$Z==inference$W))/nrow(inference)
    itt.root <- iv.effect.root*compliers.root
    se_ccace.root <- summary$coef[2,2]   # Standard error at root
    
    # Compute 95% Confidence Interval
    ci_lower_root <- iv.effect.root - 1.96 * se_ccace.root
    ci_upper_root <- iv.effect.root + 1.96 * se_ccace.root
    bcfivMat[1,] <- c(NA, round(iv.effect.root, 4),
                      round(p.value.root, 4),
                      round(p.value.weak.iv.root, 4),
                      round(proportion.root, 4),
                      round(itt.root, 4),
                      round(compliers.root, 4),
                      round(se_ccace.root, 4),
                      round(ci_lower_root, 4),
                      round(ci_upper_root, 4))
    names(inference) <- paste(names(inference), sep="")
    # Run a loop to get the rules (sub-populations)
    
    for (i in rules[-1]){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree, node=i, print.it = FALSE)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i,print.it=FALSE))[j+1,1]))
        sub_pop <- noquote(paste(sub , collapse = " & "))
      }))
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      if (length(unique(subset$W))!= 1 | length(unique(subset$Z))!= 1){
        iv.reg <- ivreg(Y ~ W | Z,  
                        data = subset)
        summary <- summary(iv.reg, diagnostics = TRUE)
        iv.effect <- summary$coef[2, 1]
        se_ccace <- summary$coef[2, 2]
        ci_lower <- iv.effect - 1.96 * se_ccace
        ci_upper <- iv.effect + 1.96 * se_ccace
        
        p.value <- summary$coef[2, 4]
        p.value.weak.iv <- summary$diagnostics[1, 4]
        compliers <- length(which(subset$Z == subset$W)) / nrow(subset)
        itt <- iv.effect * compliers
        proportion.node <- nrow(subset)/nrow(inference)
        
        bcfivMat[i, ] <- c(sub_pop,
                           round(iv.effect, 4),
                           round(p.value, 4),
                           round(p.value.weak.iv, 4),
                           round(proportion.node, 4),
                           round(itt, 4),
                           round(compliers, 4),
                           round(se_ccace, 4),
                           round(ci_lower, 4),
                           round(ci_upper, 4))
      }
      # Delete data
      rm(subset)  
    } 
    # Adjust P.values 
    bcfiv_correction <- cbind(as.data.frame(bcfivMat), leaves)
    adj <- round(p.adjust( as.numeric(bcfiv_correction$pvalue[which(bcfiv_correction$leaves==1)]) ,  paste(adj_method)), 5)
    Adj_pvalue <- rep(NA, length(rules)) 
    Adj_pvalue[which(bcfiv_correction$leaves==1)] <- adj
    
    # Store Results
    bcfivResults <- cbind(as.data.frame(bcfivMat), Adj_pvalue)
  }    
  

  #Binary Outcomes 

  if (binary == TRUE){
    # Perform the binary Bayesian Causal Forest for the ITT
    fit_itt <- quiet(bartCause::bartc(discovery$Y, as.numeric(discovery$Z), as.matrix(discovery_X), n.samples = 300, n.burn = 300, n.chains = 2L))
    
    # Get posterior of treatment effects
    ites <- bartCause::extract(fit_itt, type = "ite")
    itt <- apply(ites, 2, mean)
    
    # Get posterior of treatment effects
    epsilon <- 1e-6
    tauhat <- itt / (pic+epsilon)
    exp <- as.data.frame(cbind(tauhat, discovery_X))
    
    #Build a CART
    fit.tree <- rpart(tauhat ~ .,
                      data = exp,
                      maxdepth = 2,
                      cp=0.001)
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Initialize Outputs
    bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol=7))
    names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test", "Pi_obs", "ITT", "Pi_compliers")
    
    # Generate Leaves Indicator
    lvs <- leaves <- numeric(length(rules)) 
    lvs[unique(fit.tree$where)] <- 1
    leaves[rules[lvs==1]] <- 1
    
    #Run an IV Regression on each Node
    iv.root <- AER::ivreg(Y ~ W | Z,  
                          data = inference)
    summary <- summary(iv.root, diagnostics = TRUE)
    iv.effect.root <-  summary$coef[2,1]
    p.value.root <- summary$coef[2,4]
    p.value.weak.iv.root <- summary$diagnostics[1,4]
    proportion.root <- 1
    compliers.root <- length(which(inference$Z==inference$W))/nrow(inference)
    itt.root <- iv.effect.root*compliers.root
    
    # Store Results for Root
    bcfivMat[1,] <- c( NA , round(iv.effect.root, 4), round(p.value.root, 4), round(p.value.weak.iv.root, 4), round(proportion.root, 4), round(itt.root, 4), round(compliers.root, 4))
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules[-1]){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree, node=i, print.it = FALSE)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i,print.it=FALSE))[j+1,1]))
        sub_pop <- noquote(paste(sub , collapse = " & "))
      }))
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      # Run the IV Regression
      if (length(unique(subset$W))!= 1 | length(unique(subset$Z))!= 1){
        iv.reg <- AER::ivreg(Y ~ W | Z,  
                             data = subset)
        summary <- summary(iv.reg, diagnostics = TRUE)
        iv.effect <-  summary$coef[2,1]
        p.value <- summary$coef[2,4]
        p.value.weak.iv <- summary$diagnostics[1,4]
        compliers <- length(which(subset$Z==subset$W))/nrow(subset)
        itt <- iv.effect*compliers
        
        # Proportion of observations in the node
        proportion.node <- nrow(subset)/nrow(inference)
        bcfivMat[i,] <- c(sub_pop, round(iv.effect, 4), round(p.value, 4), round(p.value.weak.iv, 4), round(proportion.node, 4), round(itt, 4), round(compliers, 4))
      }
      
      # Delete data
      rm(subset)
    }
    
    # Adjust P.values 
    bcfiv_correction <- cbind(as.data.frame(bcfivMat), leaves)
    adj <- round(p.adjust( as.numeric(bcfiv_correction$pvalue[which(bcfiv_correction$leaves==1)]) ,  paste(adj_method)), 5)
    Adj_pvalue <- rep(NA, length(rules)) 
    Adj_pvalue[which(bcfiv_correction$leaves==1)] <- adj
    # Store Results
    bcfivResults <- cbind(as.data.frame(bcfivMat), Adj_pvalue)
  }
  # Return Results
  return(list(results = bcfivResults, discovery = discovery, inference = inference,pic=pic))
  
}
GRF_IV <- function(Y, W, Z, X) {
  p.score <- glm(Z ~ ., family = binomial, data = X)
  pihat <- predict(p.score, type = "response")
  valid_idx <- which(pihat > 0.1 & pihat < 0.9)
  subset_data <- list(
    Y = Y[valid_idx],
    W = W[valid_idx],
    Z = Z[valid_idx],
    X = X[valid_idx, ],
    pihat = pihat[valid_idx]
  )
  
  set.seed(42)  # For reproducibility
  
  # Number of valid observations
  n_valid <- length(subset_data$Y)
  
  # Create random indices for training 
  discovery_id <- sample(n_valid, size = floor(0.5 * n_valid))  # 50% for training
  inference_id <- setdiff(1:n_valid, discovery_id)  #  50% for validation
  
  # Create train/validation splits
  discovery <- list(
    Y = subset_data$Y[discovery_id],
    W = subset_data$W[discovery_id],
    Z = subset_data$Z[discovery_id],
    X = subset_data$X[discovery_id, ],
    pihat = subset_data$pihat[discovery_id]
  )
  
  inference <- list(
    Y = subset_data$Y[inference_id],
    W = subset_data$W[inference_id],
    Z = subset_data$Z[inference_id],
    X = subset_data$X[inference_id, ],
    pihat = subset_data$pihat[inference_id]
  )
  inference_combined <- data.frame(
    Y = inference$Y,
    W = inference$W,
    Z = inference$Z,
    inference$X,  # This expands the columns of X
    pihat = inference$pihat
  )

  iv_forest <- grf::instrumental_forest(X=as.data.frame(discovery$X), Y=discovery$Y, W=discovery$W, Z=discovery$Z, num.trees = 2000, min.node.size=20)
  tauhat <- predict(iv_forest)$predictions
 
  exp <- as.data.frame(cbind(tauhat, as.data.frame(discovery$X)))
  fit.tree <- rpart(tauhat ~ .,
                    data = exp,
                    maxdepth = 2,
                    cp=0.001)
  rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
  
  # Initialize Outputs
  ivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol = 10))
  names(ivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test",
                    "Pi_obs", "ITT", "Pi_compliers",
                    "SE", "CI_lower", "CI_upper")
  lvs <- leaves <- numeric(length(rules)) 
  lvs[unique(fit.tree$where)] <- 1
  leaves[rules[lvs==1]] <- 1
  # Run an IV Regression on the Root
  iv.root <- AER::ivreg(Y ~ W | Z,  
                        data = inference)
  summary <- summary(iv.root, diagnostics = TRUE)
  iv.effect.root <-  summary$coef[2,1]
  p.value.root <- summary$coef[2,4]
  p.value.weak.iv.root <- summary$diagnostics[1,4]
  proportion.root <- 1
  compliers.root <- length(which(inference$Z==inference$W))/nrow(inference$X)
  itt.root <- iv.effect.root*compliers.root
  
  se_ccace.root <- summary$coef[2,2]   # Standard error at root
  
  # Compute 95% Confidence Interval
  ci_lower_root <- iv.effect.root - 1.96 * se_ccace.root
  ci_upper_root <- iv.effect.root + 1.96 * se_ccace.root
  
  ivMat[1,] <- c(NA, round(iv.effect.root, 4),
                 round(p.value.root, 4),
                 round(p.value.weak.iv.root, 4),
                 round(proportion.root, 4),
                 round(itt.root, 4),
                 round(compliers.root, 4),
                 round(se_ccace.root, 4),
                 round(ci_lower_root, 4),
                 round(ci_upper_root, 4))
  names(inference) <- paste(names(inference), sep="")
  
  # Run a loop to get the rules (sub-populations)
  for (i in rules[-1]) {
    sub <- as.data.frame(matrix(NA, nrow = 1,
                                ncol = nrow(as.data.frame(path.rpart(fit.tree, node = i, print.it = FALSE))) - 1))
    quiet(capture.output(for (j in 1:ncol(sub)) {
      sub[, j] <- as.character(print(as.data.frame(path.rpart(fit.tree, node = i, print.it = FALSE))[j + 1, 1]))
      sub_pop <- noquote(paste(sub, collapse = " & "))
    }))
    
    subset <- with(inference_combined, inference_combined[which(eval(parse(text = sub_pop))), ])
    
    if (nrow(subset) > 3 &&
        length(unique(subset$W)) > 1 &&
        length(unique(subset$Z)) > 1 &&
        all(!is.na(subset$Y)) && all(!is.na(subset$W)) && all(!is.na(subset$Z))) {
      
      X_model <- tryCatch(model.matrix(~ W | Z, data = subset), error = function(e) return(NULL))
      if (is.null(X_model) || qr(X_model)$rank < ncol(X_model)) {
        message("Skipping node: ", sub_pop, " due to rank-deficient model matrix.")
        next
      }
      
      iv.reg <- tryCatch({
        AER::ivreg(Y ~ W | Z, data = subset)
      }, error = function(e) {
        message("IV regression failed for node: ", sub_pop, " - Skipping...")
        return(NULL)
      })
      
      if (!is.null(iv.reg)) {
        summary <- tryCatch({
          summary(iv.reg, diagnostics = TRUE)
        }, error = function(e) {
          message("Summary failed for node: ", sub_pop, " - Skipping...")
          return(NULL)
        })
        
        if (!is.null(summary) && summary$df[2] > 0) {
          iv.effect <- summary$coef[2, 1]
          se_ccace <- summary$coef[2, 2]
          ci_lower <- iv.effect - 1.96 * se_ccace
          ci_upper <- iv.effect + 1.96 * se_ccace
          
          p.value <- summary$coef[2, 4]
          p.value.weak.iv <- summary$diagnostics[1, 4]
          compliers <- length(which(subset$Z == subset$W)) / nrow(subset)
          itt <- iv.effect * compliers
          proportion.node <- nrow(subset) / nrow(inference$X)
          
          ivMat[i, ] <- c(sub_pop,
                          round(iv.effect, 4),
                          round(p.value, 4),
                          round(p.value.weak.iv, 4),
                          round(proportion.node, 4),
                          round(itt, 4),
                          round(compliers, 4),
                          round(se_ccace, 4),
                          round(ci_lower, 4),
                          round(ci_upper, 4))
        } else {
          message("Skipping node: ", sub_pop, " due to zero residual degrees of freedom.")
        }
      }
    } else {
      message("Skipping node: ", sub_pop, " due to missing values, small sample size, or lack of variation.")
    }
  }    
  
  # Adjust P.values 
  iv_correction <- cbind(as.data.frame(ivMat), leaves)
  adj <- round(p.adjust( as.numeric(iv_correction$pvalue[which(iv_correction$leaves==1)]) ,  paste(adj_method)), 5)
  Adj_pvalue <- rep(NA, length(rules)) 
  Adj_pvalue[which(iv_correction$leaves==1)] <- adj
  ivResults <- cbind(as.data.frame(ivMat), Adj_pvalue)
  # Return Results
  return(list(results = ivResults, discovery = discovery, inference = inference))
  
}
