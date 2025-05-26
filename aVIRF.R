### Packages ###
source("Helpers.R")

### Server ###
load("BEKK_model_and_data.RData")
#load("VIRFs/bootstrap_bekk_result.RData")

### Parameters ###
n.ahead <- 250
timeforVIRF <- 444

# --------------------------------------
#  AVIRF  –  Asymmetric / Symmetric VIRF
# --------------------------------------
AVIRF_parallel <- function(data,                 # returns
                           C, A, G,              # BEKK-Parameter
                           B = NULL,             # only for asymmetric model
                           timeforVIRF,          # start index t
                           n.ahead,              # h-ahead forecast
                           asym = FALSE,         # Using sym / asym?
                           simsamp = NULL,      # simulation sample size
                           C_boot = NULL, A_boot = NULL,
                           G_boot = NULL, B_boot = NULL,
                           cores = detect.cores(),
                           seed = 123)
{
  ## -- 1  Checks --------------------------------------------------------------
  if (is.null(B) && asym){stop("Asym-Modell benötigt B!")}
  if (is.null(simsamp)){simsamp <- 2000 * n.ahead}
  if (abs(mean(data[, 1])) > 1e-15 || abs(mean(data[, 2])) > 1e-15) {
    data <- scale(data, center = TRUE, scale = FALSE)
    message("Data was centered")
  }
  
  ## -- 2  BEKK-Reconstruction--------------------------------------------------
  N  <- nrow(data)
  H_t <- matrix(0, N, 4)
  xi  <- matrix(0, N, 2)
  
  H_t[1, ] <- t(data) %*% data / N
  for (i in 2:N) {
    eta <- calculate_eta(data[i - 1, ])
    H_t[i, ] <- as.vector(C %*% t(C) +
                            t(A) %*% data[i - 1, ] %*% t(data[i - 1, ]) %*% A +
                            t(G) %*% matrix(H_t[i - 1, ], 2, 2) %*% G)
    if (asym){H_t[i, ] <- H_t[i, ] + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
    xi[i, ] <- solve(matroot(matrix(H_t[i, ], 2, 2))) %*% data[i, ]
  }
  
  ## -- 3  Optional: Only relevant for bootstrap CI ----------------------------
  if (!is.null(C_boot) | !is.null(A_boot) | !is.null(G_boot) | !is.null(B_boot)) {
    C <- C_boot; A <- A_boot; G <- G_boot; B <- B_boot
  }
  
  ## -- 4  Preparation simulation ----------------------------------------------
  library(parallel)
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)                               
  seed_vec <- sample.int(.Machine$integer.max, simsamp)
  
  # Ergebnis-Container
  VIRF <- array(0, dim = c(n.ahead, 3, simsamp))
  CIRF <- matrix(0, n.ahead, simsamp)
  
  ## -- 5  sim with mcapply ----------------------------------------------------
  res_list <- mclapply(seq_len(simsamp), function(j) {
    if (j %% 10000 == 0) {
      cat("Iteration", j, "done\n")
    }
    set.seed(seed_vec[j])                          # eindeutiger Worker-Seed
    
    # a) Stichproben der standardisierten Residuen
    xisamp0 <- xi[sample(N, n.ahead, replace = TRUE), , drop = FALSE]
    xisamp1 <- xi[sample(N, n.ahead, replace = TRUE), , drop = FALSE]
    
    # b) Initiale Kovarianzmatrizen
    H_0 <- matrix(H_t[timeforVIRF,     ], 2, 2)
    H_1 <- matrix(H_t[timeforVIRF + 1, ], 2, 2)
    
    VIRF_j <- matrix(0, n.ahead, 3)
    CIRF_j <- numeric(n.ahead)
    
    # ----- Schritt 1 ----------------------------------------------------------
    e_0  <- matroot(H_0) %*% xisamp0[1, ]
    eta  <- calculate_eta(e_0)
    H_0p <- C %*% t(C) + t(A) %*% e_0 %*% t(e_0) %*% A + t(G) %*% H_0 %*% G
    if (asym){H_0p <- H_0p + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
    
    VIRF_j[1, ] <- vech(H_1 - H_0p)
    CIRF_j[1]   <- (H_1[1, 2] / sqrt(H_1[1, 1] * H_1[2, 2])) -
      (H_0p[1, 2] / sqrt(H_0p[1, 1] * H_0p[2, 2]))
    
    # ----- Schritte 2 … h ------------------------------------------------------
    for (i in 2:n.ahead) {
      # „0-Pfad“
      e_0  <- matroot(H_0p) %*% xisamp0[i, ]
      eta  <- calculate_eta(e_0)
      H_0p <- C %*% t(C) + t(A) %*% e_0 %*% t(e_0) %*% A + t(G) %*% H_0p %*% G
      if (asym){H_0p <- H_0p + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
      
      # „1-Pfad“
      e_1  <- matroot(H_1) %*% xisamp1[i, ]
      eta  <- calculate_eta(e_1)
      H_1  <- C %*% t(C) + t(A) %*% e_1 %*% t(e_1) %*% A + t(G) %*% H_1 %*% G
      if (asym){H_1 <- H_1 + as.vector(t(B) %*% eta %*% t(eta) %*% B)}
      
      VIRF_j[i, ] <- vech(H_1 - H_0p)
      CIRF_j[i]   <- (H_1[1, 2] / sqrt(H_1[1, 1] * H_1[2, 2])) -
        (H_0p[1, 2] / sqrt(H_0p[1, 1] * H_0p[2, 2]))
    }
    
    list(VIRF = VIRF_j, CIRF = CIRF_j)             # Rückgabe für diesen Bootstrap
  }, mc.cores = cores)
  
  ## -- 6  Results -------------------------------------------
  for (j in seq_along(res_list)) {
    VIRF[ , , j] <- res_list[[j]]$VIRF
    CIRF[ ,   j] <- res_list[[j]]$CIRF
  }
  
  ## -- 7  
  VIRF_mean <- apply(VIRF, c(1, 2), mean)
  # VIRF_q_l <- apply(VIRF, c(1,2), quantile, probs=c(0.05), na.rm=TRUE)
  # VIRF_q_u <- apply(VIRF, c(1,2), quantile, probs=c(0.95), na.rm=TRUE)
  # VIRF_sd <- apply(VIRF, c(1,2), sd)
  CIRF_mean <- apply(CIRF, 1, mean)
  
  list(VIRF_mean = VIRF_mean,
       CIRF_mean = CIRF_mean,
       # VIRF_q_l = VIRF_q_l,
       # VIRF_q_u = VIRF_q_u,
       # VIRF_q_sd = VIRF_sd,
       VIRF      = VIRF,
       CIRF      = CIRF)
}

start <- Sys.time()
VIRF_result <- AVIRF(data = data,
                     C = bekk_m2$C0,
                     G = bekk_m2$G,
                     A = bekk_m2$A,
                     timeforVIRF = timeforVIRF,
                     n.ahead = n.ahead,
                     asym = FALSE,
                     simsamp = 100000,
                     cores = detectCores(),
                     seed = 123)

cat("Time: ", Sys.time()- start)
