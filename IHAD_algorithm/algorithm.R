#---------- Quasi monte carlo with Van der Corput sequence ---------------------#
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

#-------------------------- ADMM with given support ----------------------------#
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

#------------------------------- IHT -------------------------------------------#

ft_IHT_update <- function(Lambda_0,Sigma,omega,gma,eta) {
  n <- nrow(Sigma)
  
  grad <- -2 * ((Sigma %*% Lambda_0) %*% (Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n))) / frobenius.norm(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n))
  grad_step <- Lambda_0 - eta * grad
  Lambda_updt <- grad_step
  
  mask_diag <- diag(n) == 1        # diagonal entries
  mask_off <- !mask_diag
  
  
  mask_keep <- abs(grad_step) > sqrt(2 * eta * gma)
  # mask_keep[2,3] <- TRUE
  Lambda_updt[mask_off & !mask_keep] <- 0
  # Lambda_updt[c(2,3,4)] <- 0
  
  return (list(Lambda_updt,grad))
}

ft_IHT_bktr <- function(Lambda_0,Sigma,omega,gma,eta_0, eta_inc) {
  n <- nrow(Sigma)
  
  # eta_0 is the initial stepsize (before backtracking)
  eta <- list(eta_0)
  
  Lambda <- list(Lambda_0)
  obj <- list(frobenius.norm(Sigma - t(Lambda_0)%*%Sigma%*%Lambda_0 - omega * diag(n)) + gma * (sum(abs(c(Lambda_0)) >0)-n))
  obj_gap <- list(obj[[1]])
  
  for (i in 1:1000) {
    eta_t <- eta[[i]]
    l <- 1
    repeat {
      IHT_updt <- ft_IHT_update(Lambda[[i]], Sigma, omega, gma, eta_t)
      Lambda_updt <- IHT_updt[[1]]
      grad_updt <- IHT_updt[[2]]
      
      lhs <- frobenius.norm(Sigma - t(Lambda_updt) %*% Sigma %*% Lambda_updt - omega * diag(n))
      rhs <- frobenius.norm(Sigma - t(Lambda[[i]]) %*% Sigma %*% Lambda[[i]] - omega * diag(n)) + sum((Lambda_updt - Lambda[[i]]) * grad_updt) + (frobenius.norm(Lambda_updt - Lambda[[i]])^2) / (2*eta_t)
      
      if (lhs <= rhs) break
      eta_t <- eta_t / eta_inc
      
      l <- l+1
      if(l > 1000) break
    }
    eta[[i+1]] <- eta_t
    
    IHT_updt <- ft_IHT_update(Lambda[[i]],Sigma,omega,gma,eta[[i+1]])
    Lambda[[i+1]] <- IHT_updt[[1]]
    obj[[i+1]] <- frobenius.norm(Sigma - t(Lambda[[i+1]])%*%Sigma%*%Lambda[[i+1]] - omega * diag(n)) + gma * (sum(abs(c(Lambda[[i+1]]) )>0)-n)
    obj_gap[[i+1]] <- abs(obj[[i+1]] - obj[[i]])
    lmd_gap <- frobenius.norm(Lambda[[i+1]] - Lambda[[i]])
    
    if(obj_gap[[i+1]] < 1e-4 & lmd_gap < 1e-4) {
      obj_real <- frobenius.norm(Sigma - t(Lambda[[i+1]])%*%Sigma%*%Lambda[[i+1]] - omega * diag(n))
      return(list(Lambda[[i+1]],obj_real,obj[[i+1]],TRUE))
    }
  }
  print("IHT did not converge within 1000 iterations. Objective gap:")
  print(lmd_gap)
  obj_real <- frobenius.norm(Sigma - t(Lambda[[i+1]])%*%Sigma%*%Lambda[[i+1]] - omega * diag(n))
  return(list(Lambda[[i+1]],obj_real,obj[[i+1]],FALSE))
}

#---------------- Construct a list of possible graphs --------------------------#

run_IHT_chains <- function(l=10, Lambda_init, Sigma, omega) {
  res_list <- vector("list", l)
  
  for (i in seq_len(l)) {
    cat("Running IHT_chain for Lambda_init[[", i, "]]...\n")
    res_IHT <- IHT_chain(Lambda_init[[i]], Sigma, omega)
    res_list[[i]] <- res_IHT
  }
  
  return(res_list)
}

merge_IHT_supports <- function(res_list) {
  # Flatten all Lambda matrices from the chains
  Lambda_all <- unlist(res_list, recursive = FALSE)
  
  # Step 1: Extract support masks (TRUE where nonzero)
  support_list <- lapply(Lambda_all, function(L) (L != 0))
  
  # Step 2: Convert each support into a unique string key
  support_keys <- sapply(support_list, function(M) paste(as.integer(M), collapse = ""))
  
  # Step 3: Identify unique supports
  unique_support_keys <- unique(support_keys)
  unique_indices <- match(unique_support_keys, support_keys)
  
  # Step 4: Keep only the unique supports
  support_unique <- support_list[unique_indices]
  
  # Step 5: Sort by number of nonzero elements (increasing)
  nonzero_counts <- sapply(support_unique, sum)
  order_idx <- order(nonzero_counts)
  
  support_sorted <- support_unique[order_idx]
  
  # Step 6: Convert TRUE/FALSE → 1/0
  merged_support_list_num <- lapply(support_sorted, function(M) {
    M_num <- matrix(as.numeric(M), nrow = nrow(M))
    return(M_num)
  })
  
  cat("Merged", length(merged_support_list_num), "unique supports (sorted by sparsity).\n")
  
  return(merged_support_list_num)
}

#---------------------- Run ADMM with the constructed support list -------------#

run_ADMM_on_supports <- function(Sigma, omega, merged_support_list_num) {
  res_ADMM <- NULL
  
  for (i in seq_along(merged_support_list_num)) {
    cat("Running ADMM_supp for support", i, "of", length(merged_support_list_num), "...\n")
    
    res_ADMM <- ADMM_supp(Sigma, omega, merged_support_list_num[[i]])
    
    # Check if the first element (flag) is TRUE
    if (isTRUE(res_ADMM[[1]])) {
      cat("✓ ADMM succeeded at support index", i, "\n")
      print(res_ADMM)
      return(res_ADMM)
    }
  }
  
  cat("No valid ADMM result found.\n")
  return(res_ADMM)  # Will return last (possibly failed) result or NULL
}




