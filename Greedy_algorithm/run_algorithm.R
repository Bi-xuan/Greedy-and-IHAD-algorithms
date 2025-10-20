# Tests on identifiable graphs
# For n small
# -------- Library Loading --------
suppressPackageStartupMessages({
  library('fBasics')
  library('Matrix')
  library("matrixcalc")
  library(future.apply)
  # library(progressr)
})

# -------- Parallel setup --------
# This sets up parallel processing to speed up the computations.
plan(multisession, workers = parallel::detectCores() - 1)

# Enable progress bar globally
# handlers(global = TRUE)
# handlers("txtprogressbar") # simple text bar works in background jobs

# -------- Data Loading --------
# Load the ground truth Lambda matrices
source("Greedy_algorithm/algorithm.R")
load("Lambda_star/trees/n_3_trees_sig_0_2.RData")

n <- ncol(Lambda_star[[1]])
omega <- 1 # The true omega
omega_est <- 1 # Changeable parameter, using the default

# -------- Parallel setup --------
plan(multisession, workers = parallel::detectCores() - 1)

n <- ncol(Lambda_star[[1]])

# -------- Main Processing Loop --------
n_inits <- length(Lambda_star)
# n_inits <- 10  # For testing
cat("Starting analysis for", n_inits, "Lambda matrices...\n")

# A progressor allows the parallel loop to have a progress bar
# p <- progressor(steps = n_inits)

results_list <- future_lapply(1:n_inits, function(i) {
  # Get the true Lambda from the loaded data
  Lambda_star_curr <- Lambda_star[[i]]
  
  # 1. Calculate Sigma based on the model equation
  vec_Sigma <- solve(diag(n**2)-(t(Lambda_star_curr) %x% t(Lambda_star_curr))) %*% vec(omega * diag(n))
  Sigma <- matrix(vec_Sigma, nrow = n)
  
  # Initialize an empty result structure
  result <- list(
    lambda_index = i,
    converged = FALSE,
    objective_value = NA,
    ground_truth = Lambda_star_curr,
    estimated_lambda = matrix(NA, n, n),
    runtime_secs = NA,
    error_message = ""
  )
  
  # Proceed only if Sigma could be calculated and is positive definite
  if (!is.null(Sigma) && is.positive.definite(Sigma)) {
    # 2. Run the greedy algorithm
    start_time <- Sys.time()
    tryCatch({
      greedy_result <- algo_greedy(Sigma, omega = omega_est)
      end_time <- Sys.time()
      
      # 3. Store the results
      result$converged <- greedy_result[[1]]
      result$estimated_lambda <- greedy_result[[2]]
      result$objective_value <- greedy_result[[3]]
      result$runtime_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
    }, error = function(e) {
      result$error_message <- e$message
    })
  } else {
    result$error_message <- "Generated Sigma was NULL or not positive definite."
  }
  
  # Save partial result
  partial_file <- sprintf("Greedy_algorithm/results/partial_result/result_%03d.RData", i)
  save(result, file = partial_file)
  
  # p() # This updates the progress bar
  return(result)
}, future.seed = TRUE) # future.seed = TRUE makes parallel processing reproducible

cat("\nAnalysis complete. Consolidating results...\n")

# -------- Result Consolidation --------
# Convert the list of results into a dataframe
results_df <- do.call(rbind, lapply(results_list, function(res) {
  data.frame(
    lambda_index = res$lambda_index,
    converged = res$converged,
    objective_value = res$objective_value,
    runtime_secs = res$runtime_secs,
    error_message = res$error_message
  )
}))

# Storing matrices in a dataframe is best done using list-columns
results_df$true_lambda <- lapply(results_list, `[[`, "ground_truth")
results_df$estimated_lambda <- lapply(results_list, `[[`, "estimated_lambda")

# Add a metric to see if the sparsity pattern (support) was correctly identified
results_df$support_correct <- mapply(function(true, est) {
  # The check needs to handle the case where est is NA
  if(any(is.na(est))) return(FALSE)
  # Compare the boolean matrices representing where non-zero elements are
  all((true != 0) == (est != 0))
}, results_df$true_lambda, results_df$estimated_lambda)

# -------- Display result --------
cat("\n--- Results Summary (first 10 entries) ---\n")
print(results_df[1:10, c("lambda_index", "converged", "support_correct", "objective_value", "runtime_secs", "error_message")])

# -------- Store result --------
save(results_df, file = "Greedy_algorithm/results/trees/n_3_ADMM_VdC_10inits_sig_0_2.RData")



