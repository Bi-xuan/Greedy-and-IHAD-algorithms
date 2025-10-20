# Quasi monte carlo with Van der Corput sequence.
# --- Van der Corput (1D base sequence) ---
van_der_corput <- function(n, base = 2) {
  seq <- numeric(n)
  for (i in 1:n) {
    x <- i - 1
    f <- 1
    r <- 0
    while (x > 0) {
      f <- f / base
      r <- r + f * (x %% base)
      x <- floor(x / base)
    }
    seq[i] <- r
  }
  return(seq)
}

# --- Halton sequence: n points in 'dim' dimensions ---
halton_sequence <- function(n, dim) {
  # Generate the first 'dim' prime numbers
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
              53, 59, 61, 67, 71, 73, 79, 83, 89, 97)
  if (dim > length(primes)) stop("Increase number of available primes")
  
  seq <- matrix(0, n, dim)
  for (d in 1:dim) {
    seq[, d] <- van_der_corput(n, base = primes[d])
  }
  return(seq)
}

# ADMM with given support
ADMM_supp <- function(Sigma,omega,Lambda_t,beta=0.5) {
  n <- nrow(Sigma)
  for (j in 1:10) {
    # Initialization
    halton_mat <- halton_sequence(n = j, dim = n^2)  # get first j Halton points
    halton_vec <- as.numeric(halton_mat[j, ])        # use the j-th quasi-random point
    halton_vec <- 2 * (halton_vec - 0.5)             # scale to [-1, 1]
    halton_mat <- matrix(halton_vec, n, n)
    
    Lambda_1_prev <- halton_mat
    Lambda_2_prev <- halton_mat
    alpha_prev <- matrix(0, n, n)
    
    # Pre-calculate constant part
    Sigma_minus_omega <- Sigma - omega * diag(n)
    
    
    dis_1_2 <- list(1)
    obj <- list(1)
    
    for (i in 1:1000) {
      # Update Step
      A <- 2*(Sigma %*% Lambda_2_prev %*% t(Lambda_2_prev) %*% Sigma) + beta * diag(n)
      Lambda_1_curr <- chol2inv(chol(A)) %*% (2 * (Sigma %*% Lambda_2_prev %*% Sigma_minus_omega) - alpha_prev + beta * Lambda_2_prev)   # Use the more stable and more efficient chol2inv() function instead of solve()
      B <- 2*Sigma %*% Lambda_1_curr %*% t(Lambda_1_curr) %*% Sigma + beta * diag(n)
      Lambda_2_curr <- chol2inv(chol(B)) %*% (2 * Sigma %*% Lambda_1_curr %*% Sigma_minus_omega + alpha_prev + beta * Lambda_1_curr)
      
      # Apply sparsity
      Lambda_1_curr[Lambda_t == 0] <- 0
      Lambda_2_curr[Lambda_t == 0] <- 0
      
      # Update alpha
      alpha_curr <- alpha_prev + beta * (Lambda_1_curr - Lambda_2_curr)
      
      # Check convergence
      dis_1_2 <- norm(Lambda_1_curr - Lambda_2_curr,type="F")
      obj <- norm(Sigma - t(Lambda_1_curr) %*% Sigma %*% Lambda_2_curr - omega * diag(n),type="F")**2
      
      
      if(dis_1_2 < 1e-6 && obj < 1e-6) {
        # if(dis_1_2 < 1e-6) {
        return(list(TRUE,Lambda_1_curr,obj)) # Return the converged \Lambda and the objective function value
      }
      
      # Prepare for next iteration
      Lambda_1_prev <- Lambda_1_curr
      Lambda_2_prev <- Lambda_2_curr
      alpha_prev <- alpha_curr
    }
  }
  return(list(FALSE,Lambda_1_curr,obj))
}