generate_dgp_continous_sim <- function(sim=1, n, k, heterogeneity = "strong") {
  set.seed(123)  # For reproducibility
  # Initialize an empty list to store results from each simulation
  all_data <- vector("list", sim)
  
  for (s in 1:sim) {
    X <- matrix(rnorm(n * 10, mean = 0, sd = 1), ncol = 10)
    colnames(X) <- paste0("X", 1:10)
    Z <- rbinom(n, size = 1, prob = 0.5)
    compliance_rate<-0.75
    # Generate potential treatments
    W_0 <- rep(0, n)  # one-sided noncompliance
    W_1 <- rbinom(n, size = 1, prob = compliance_rate) 
    Y_0 <- rnorm(n, mean = 0, sd = 1)  
    # Define tau
    if (heterogeneity == "strong") {
      tau_cace <- ifelse(X[, 1] < -0.5 & X[, 2] < -0.5,  k,          # L1
                         ifelse(X[, 1] >  0.5 & X[, 2] >  0.5, -k, 0))      # L2, else L0
      
    } else if (heterogeneity == "slight") {
      tau_cace <- ifelse(X[, 1] < 0 & X[, 2] < 0,        k,          # L1
                         ifelse(X[, 1] > 0 & X[, 2] > 0,       -k,          # L2
                                ifelse(X[, 1] > 0 & X[, 2] < 0,   0.5 * k,         # L3
                                       -0.5 * k)))     # L4
      
    } else {
      stop("Invalid heterogeneity scenario. Choose 'strong' or 'slight'.")
    }
    Y_1 <- Y_0 + W_1 * tau_cace
    
    W <- ifelse(Z == 1, W_1, W_0)  
    Y <- ifelse(Z == 1, Y_1, Y_0)  
    im_data <- data.frame(sim_id = s, Y = Y, W = W, CACE = tau_cace, Z = Z, X)
    all_data[[s]] <- sim_data
  }
  # Combine all simulated datasets into a single data frame
  final_data <- do.call(rbind, all_data)
  return(final_data)
}




dgp_slight_1 <- generate_dgp_continous_sim(sim = 500,n = 10000,k = 1,heterogeneity = "slight")

saveRDS(
  dgp_slight_1,
  file = here::here(
    "data",
    "N_10000",
    "dgp_slight_1.rds"
  )
)