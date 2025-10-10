# ADMM with given support
ADMM_supp <- function(Sigma,omega,Lambda_t,beta=0.5) {
  for (j in 1:10) {
    # Initialization
    Lambda_1_prev <- diag(n)
    Lambda_2_prev <- Lambda_1_prev
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

# Greedy algorithm
algo_greedy <- function(Sigma,omega=1) {
  n <- ncol(Sigma)
  
  # List of indices of possible edges
  edge_ind <- setdiff(1:(n**2),seq(1,n**2,n+1))
  
  # Explore all possible support with num_edge - 1 and num_edge edges
  for (i in 1:length(edge_ind)) {
    edge_cand <- combn(edge_ind,i)
    for (j in seq(1,ncol(edge_cand),1)) {
      Lambda_t <- diag(n)
      Lambda_t[edge_cand[,j]] <- 1
      res_ADMM_supp <- ADMM_supp(Sigma,omega,Lambda_t)
      if(res_ADMM_supp[[1]]) {
        return (res_ADMM_supp)
      }
    }
  }
  print("The greedy algorithm never finds an optimal support.")
  return (res_ADMM_supp)
}