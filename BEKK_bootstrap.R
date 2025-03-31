source("Helpers.R")

bootstrap_bekk <- function(data, # Data, centered, scale-free
                           bekk_model, # Input object from bekk_fit
                           bekk_spec_model, # Input object from bekk_spec
                           bootsamp = 999){
  # Centering data
  if (abs(mean(data[,1]))>1e-15 || abs(mean(data[,2]))>1e-15) {
    data <- scale(data, center = TRUE, scale = FALSE)
    print("Data was centered")
  }
  
  # Parameter
  A <- bekk_model$A
  G <- bekk_model$G
  C <- bekk_model$C0
  B <- bekk_model$B
  if (is.null(B)){asym = FALSE}else{asym = TRUE}

  H_t <- matrix(0,nrow(data),4)
  xi <- matrix(0,nrow(data),2)

  H_t[1,] <- t(data) %*% data / nrow(data)
  for (i in 2:nrow(data)){
    eta <- calculate_eta(data[i-1,])
    H_t[i,] <- as.vector(C %*% t(C) + t(A) %*% data[i-1,] %*% t(data[i-1,]) %*% A + t(G) %*% matrix(H_t[i-1,],2,2) %*% G)
    if (asym == TRUE){H_t[i,] <- H_t[i,] + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
    xi[i,] <- solve(matroot(matrix(H_t[i,],2,2))) %*% data[i,]
  }
  
  if (asym == FALSE){
    simulated_params <- array(0, dim = c(6,2,bootsamp))}else{
      simulated_params <- array(0, dim = c(8,2,bootsamp))
    }
  
  xi1_for_sim <- xi
  xi2_for_sim <- xi
  
  # Bootstrap BEKK 
  for (k in 1:bootsamp) {
    # Simulation BEKK returns
    set.seed(123 + k)
    xi_sim <- cbind(matrix(sample(xi1_for_sim,nrow(xi), replace = TRUE),nrow(xi),1),matrix(sample(xi2_for_sim, nrow(xi), replace = TRUE),nrow(xi),1))
    
    H_t_sim <- matrix(0,nrow(data),4)
    H_t_sim[1,] <- t(data) %*% data / nrow(data)
    ret_sim <- matrix(0,nrow(data),ncol(data))
    ret_sim[1,] <- data[1,]
    for (i in 2:nrow(data)){
      eta <- calculate_eta(ret_sim[i-1,])
      H_t_sim[i,] <- as.vector(C %*% t(C) + t(A) %*% ret_sim[i-1,] %*% t(ret_sim[i-1,]) %*% A + t(G) %*% matrix(H_t_sim[i-1,],2,2) %*% G)
      if (asym == TRUE){H_t_sim[i,] <- H_t_sim[i,] + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
      ret_sim[i,] <- matroot(matrix(H_t_sim[i,],2,2)) %*% xi_sim[i,]
    }
    ret_sim <- scale(ret_sim, center = TRUE, scale = FALSE)
    
    # Re-estimate BEKK
    bekk_sim <- bekk_fit(bekk_spec_model, ret_sim)
    simulated_params[1:2,,k] <- bekk_sim$C0
    simulated_params[3:4,,k] <- bekk_sim$A
    simulated_params[5:6,,k] <- bekk_sim$G
    if (asym == TRUE){
      simulated_params[7:8,,k] <- bekk_sim$B
    }
    
    if (k %% 10 == 0){
      cat("Iteration", k, "done \n")
    }
  }
  
  result <- list()
  result$C0 <- simulated_params[1:2,,]
  result$A <- simulated_params[3:4,,]
  result$G <- simulated_params[5:6,,]
  if (asym == TRUE){
  result$B <- simulated_params[7:8,,]
  }
  
  return(result)
}

bootstrap_bekk_result <- bootstrap_bekk(data,
                                        bekk_model = bekk_m2,
                                        bekk_spec_model = bekk_s2,
                                        bootsamp = 999)

# save(bootstrap_bekk_result, file = "bootstrap_bekk_result.RData")
# load("bootstrap_bekk_result.RData")
