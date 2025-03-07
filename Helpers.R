################
### Packages ###
################
library(BEKKs)
library(matrixcalc)

########################
### Helper Functions ###
########################
matroot <- function(mat) {
  ev <- eigen(mat, symmetric = TRUE)
  
  sqrt_vals <- sqrt(ev$values)
  Lambda <- diag(sqrt_vals)
  
  root <- ev$vectors %*% Lambda %*% t(ev$vectors)
  return(root)
}

calculate_eta <- function(e_t) {
  eta_t <- rep(0, length(e_t))
  if (all(e_t < 0)){
    eta_t <- e_t
  }
  return(eta_t)
}
