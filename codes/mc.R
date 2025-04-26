###################################################
trueVar_MC_syspps <- function(M, pdata, n){
  muGR  <- numeric(M)
  y <- pdata$y
  x <- pdata$x
  z <- pdata$z
  N <- nrow(pdata)
  mux <- mean(x)
  
  for (m in seq_len(M)) {
    sam   <- syspps(z, n)
    ys    <- y[sam]
    xs    <- x[sam]
    pis   <- n * z[sam] / sum(z)      # π_i for PPS sys sample
    di    <- 1 / pis
    muYHT <- sum(di * ys) / N
    muXHT <- sum(di * xs) / N
    muGR[m] <- (muYHT / muXHT) * mux
    if (m %% 100 == 0)
      message(" iterations completed: ", m, "/", M)
  }
  
  ##--- Monte-Carlo estimate of the true design variance -----
  trueVar_MC <- var(muGR)           # empirical variance
  return(trueVar_MC)
}

trueVar_MC_syspps(1000, pdata, 500)

##############################################################################
## 1.  Poisson PPS  ––  Monte-Carlo design variance of  μ̂_y,GR
##############################################################################
trueVar_MC_poisson <- function(M, pdata, n) {
  
  y   <- pdata$y
  x   <- pdata$x
  z   <- pdata$z                 # size variable
  N   <- length(y)
  mux <- mean(x)
  
  muGR <- numeric(M)
  
  ## first–order inclusion probabilities for Poisson PPS
  piN <- n * z / sum(z)
  
  for (m in seq_len(M)) {
    
    ## ---- draw Poisson sample ---------------------------------------------
    sam <- which(runif(N) <= piN)     # selected unit labels
    if (length(sam) == 0L) next       # (extremely rare) empty sample
    
    ys  <- y[sam]
    xs  <- x[sam]
    pis <- piN[sam]
    di  <- 1 / pis                    # Horvitz-Thompson weights
    
    ## ---- HT estimates of μ_y and μ_x -------------------------------------
    muYHT <- sum(di * ys) / N
    muXHT <- sum(di * xs) / N
    
    ## ---- Generalised ratio estimator -------------------------------------
    muGR[m] <- (muYHT / muXHT) * mux
    
    if (m %% 100 == 0)
      message("Poisson PPS – iterations completed: ", m, "/", M)
  }
  
  return(var(muGR, na.rm = TRUE))            # Monte-Carlo design variance
}

trueVar_MC_poisson(10000, pdata, 500)
##############################################################################
## 2.  PPS *with replacement* (Hansen–Hurwitz) –– Monte-Carlo design variance
##############################################################################
trueVar_MC_wrpps <- function(M, pdata, n) {
  
  y   <- pdata$y
  x   <- pdata$x
  z   <- pdata$z / sum(pdata$z)   # normalised selection probs
  N   <- length(y)
  mux <- mean(x)
  
  muGR <- numeric(M)
  
  for (m in seq_len(M)) {
    
    ## ---- draw n i.i.d. selections with replacement -----------------------
    sam <- sample.int(N, size = n, replace = TRUE, prob = z)
    
    ## accumulate Hansen–Hurwitz weights (one per *draw*)
    ys  <- y[sam]
    xs  <- x[sam]
    zs  <- z[sam]                 # selection probs for each draw
    di  <- 1 / (n * zs)
    
    ## HH estimators of μ_y and μ_x
    muYHH <- sum(di * ys) / N
    muXHH <- sum(di * xs) / N
    
    ## Generalised ratio estimator
    muGR[m] <- (muYHH / muXHH) * mux
    
    if (m %% 100 == 0)
      message("WR-PPS – iterations completed: ", m, "/", M)
  }
  
  return(var(muGR))
}

trueVar_MC_wrpps(10000, pdata, 500)
