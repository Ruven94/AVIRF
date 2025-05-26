#######################################################
### Point-wise bootstrap based confidence intervals ###
#######################################################

# This R-script gives an example of how to calculate bootstrap based CI
# BEKK_model_and_data: Load your BEKK model and data
  ## Data should be centred and scale free
  ## For BEKK model, we recommend the "BEKKs" package
  ## Use either symmetric or asymmetric model
# bootstrap_bekk_result: Load the bootstrap_bekk_result
  ## Script: BEKK_bootstrap.R

##################
### Estimation ###
##################

### Packages ###
library(parallel)
source("aVIRF.R")
num_cores <- min(detectCores() - 1, 80)

### Server ###
load("BEKK_model_and_data.RData")
load("bootstrap_bekk_result.RData")

### Parameters ###
n.ahead <- 250
timeforVIRF <- 444

R <- dim(bootstrap_bekk_result[[1]])[3]

### AVIRF ###
VIRF_result <- AVIRF(data = data,
                      C = bekk_m2$C0,
                      G = bekk_m2$G,
                      A = bekk_m2$A,
                      timeforVIRF = timeforVIRF,
                      n.ahead = n.ahead,
                      asym = FALSE,
                      simsamp = 50000)

# save(VIRF_result, file = "VIRF_result.RData")
load("Run 2/VIRF_result.RData")

### Bootstrap CI ###
result_list <- mclapply(1:R, function(k) {
  if (k %% 10 == 0) cat("Iteration", k, "von", R, "abgeschlossen\n")
  AVIRF(data = data,
        C = bekk_m2$C0,
        G = bekk_m2$G,
        A = bekk_m2$A,
        timeforVIRF = timeforVIRF,
        n.ahead = n.ahead,
        asym = FALSE,
        simsamp = 50000,
        cores = 1,
        C_boot = bootstrap_bekk_result$C0[,,k],
        A_boot = bootstrap_bekk_result$A[,,k],
        G_boot = bootstrap_bekk_result$G[,,k])
}, mc.cores = num_cores)

all_virfs <- array(0, dim = c(n.ahead, 3, R))
all_cirfs <- matrix(0, n.ahead, R)
for (k in 1:R) {
  all_virfs[,,k] <- result_list[[k]]$VIRF
  all_cirfs[,k] <- result_list[[k]]$CIRF
}

# save(all_virfs, all_cirfs, file = "VIRF_sim_result.RData")
load("Run 3/VIRF_sim_result_bootstrap.RData")
# load("Run 2/VIRF_sim_result_withoutBEKK.RData")

##################
### Evaluation ###
##################
virf_CI_lower <- matrix(0, dim(all_virfs)[1], dim(all_virfs)[2])
virf_CI_upper <- matrix(0, dim(all_virfs)[1], dim(all_virfs)[2])
virf_mean <- matrix(0, dim(all_virfs)[1], dim(all_virfs)[2])
alpha <- 0.1

for (j in 1:dim(all_virfs)[2]) {
  for (i in 1:dim(all_virfs)[1]) {
    virf_CI_lower[i,j] <- quantile(all_virfs[i,j,],alpha)
    virf_CI_upper[i,j] <- quantile(all_virfs[i,j,],1-alpha)
    virf_mean[i,j] <- mean(all_virfs[i,j,])
  }
}

# library(ggplot2)
# library(dplyr)

df_list <- list()
for (i in 1:3) {
  df_temp <- data.frame(
    index           = 1:nrow(VIRF_result$VIRF),
    response        = paste0("Response_", i),
    est             = virf_mean[, i],
    lower           = virf_CI_lower[, i],
    upper           = virf_CI_upper[, i],
    theoretical_virf = VIRF_result$VIRF[, i]   # e.g. your “true” VIRF
  )
  df_list[[i]] <- df_temp
}
df_all <- dplyr::bind_rows(df_list)

ggplot(df_all, aes(x = index, y = est)) +
  geom_line(size = 1) +  
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +  
  facet_wrap(~ response, scales = "free_y") +  
  theme_minimal() +
  labs(x = "Index", y = "Wert", 
       title = "VIRF mit Konfidenzbändern") +
  theme(plot.title = element_text(hjust = 0.5))
