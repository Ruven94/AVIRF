source("Helpers.R")

AVIRF <- function(data, # returns
                  C, # BEKK parameter C
                  A, # BEKK parameter A
                  G, # BEKK parameter G
                  B = NULL, # BEKK parameter B (set NULL if sym model)
                  timeforVIRF, # time period for VIRF
                  n.ahead, # h-steps ahead forecast
                  asym = FALSE, # using a sym oder asym model?
                  bootsamp = NULL) # Number of bootstrap samples, set NULL if it should be chosen automatically
                  {
  
  # Check conditions
  if (is.null(B) & asym == TRUE){
    stop("Asym model requires B!")
  }
  
  if (is.null(bootsamp)){
  bootsamp = 2000 * n.ahead
  }
  
  # Centering data
  if (abs(mean(data[,1]))>1e-15 || abs(mean(data[,2]))>1e-15) {
    data <- scale(data, center = TRUE, scale = FALSE)
    print("Data was centered")
  }
 
  # Reconstruction BEKK
  H_t <- matrix(0,nrow(data),4)
  xi <- matrix(0,nrow(data),2)
  
  H_t[1,] <- t(data) %*% data / nrow(data)
  for (i in 2:nrow(data)){
    eta <- calculate_eta(data[i-1,])
    H_t[i,] <- as.vector(C %*% t(C) + t(A) %*% data[i-1,] %*% t(data[i-1,]) %*% A + t(G) %*% matrix(H_t[i-1,],2,2) %*% G)
    if (asym == TRUE){H_t[i,] <- H_t[i,] + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
    xi[i,] <- solve(matroot(matrix(H_t[i,],2,2))) %*% data[i,]
  }
  
  # VIRF
  VIRF <- array(0,dim = c(n.ahead,3, bootsamp))
  CIRF <- matrix(0, n.ahead, bootsamp)
  set.seed(bootsamp)
  
  for (j in 1:bootsamp) {
    xisamp0 <- xi[sample(nrow(xi), n.ahead), , drop = FALSE]
    xisamp1 <- xi[sample(nrow(xi), n.ahead), , drop = FALSE]
    
    H_0_init <- matrix(H_t[timeforVIRF,],2,2)
    H_1_init <- matrix(H_t[timeforVIRF + 1,],2,2)
    
    # n.ahead = 1
    e_0 <- matroot(H_0_init) %*% xisamp0[1,]
    eta <- calculate_eta(e_0)
    H_0 <- C %*% t(C) + t(A) %*% e_0 %*% t(e_0) %*% A + t(G) %*% H_0_init %*% G
    if (asym == TRUE){H_0 <- H_0 + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
    
    H_1 <- H_1_init
    VIRF[1,,j] <-  vech(H_1 - H_0)
    CIRF[1,j] <- (H_1[1,2] / sqrt(H_1[1,1] * H_1[2,2])) - (H_0[1,2] / sqrt(H_0[1,1] * H_0[2,2]))
    
    # n.ahead > 1
    for (i in 2:n.ahead) {
      e_0 <- matroot(H_0) %*% xisamp0[i,]
      eta <- calculate_eta(e_0)
      H_0 <- C %*% t(C) + t(A) %*% e_0 %*% t(e_0) %*% A + t(G) %*% H_0 %*% G
      if (asym == TRUE){H_0 <- H_0 + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
      
      e_1 <- matroot(H_1) %*% xisamp1[i,]
      eta <- calculate_eta(e_1)
      H_1 <- C %*% t(C) + t(A) %*% e_1 %*% t(e_1) %*% A + t(G) %*% H_1 %*% G
      if (asym == TRUE){H_1 <- H_1 + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
      
      VIRF[i,,j] <- vech(H_1 - H_0)
      CIRF[i,j] <- (H_1[1,2] / sqrt(H_1[1,1] * H_1[2,2])) - (H_0[1,2] / sqrt(H_0[1,1] * H_0[2,2]))
    }
  }
  
  VIRF <- apply(VIRF, c(1, 2), sum)
  VIRF <- VIRF / bootsamp
  
  CIRF <- apply(CIRF, 1, sum)
  CIRF <- CIRF / bootsamp
  
  result <- list()
  result$VIRF <- VIRF
  result$CIRF <- CIRF
  
  return(result)
}

out <- AVIRF(data = data,
      C = bekk_m1$C0,
      G = bekk_m1$G,
      A = bekk_m1$A,
      B = bekk_m1$B,
      timeforVIRF = 444,
      n.ahead = 250,
      asym = TRUE)
out
