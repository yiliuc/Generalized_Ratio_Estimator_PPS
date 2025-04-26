# Read the generated data
pdata <- read.csv("data/pdata.csv")
sdata_syspps <- read.csv("data/sdata_syspps.csv")
sdata_poisson <- read.csv("data/sdata_poisson.csv")
sdata_wrpps <- read.csv("data/sdata_wrpps.csv")

################################################################################
# 1. Using HH to estimate HH for RSPPS and PPS sampling with replacement
linVarGR_HH <- function(sdata, mux, N, conf = 0.95){
  y <- sdata$y
  x <- sdata$x
  z <- sdata$zs
  n <- length(y)
  
  # HH weights & point estimate
  w  <- 1 / (n * z)         # HH weight per draw
  Ty <- sum(w * y)
  Tx <- sum(w * x)
  mu_GR <- (Ty / Tx) * mux
  Rhat  <- Ty / Tx
  
  # Linearised residuals
  e <- y - Rhat * x
  
  # HH variance for residual total
  Te_hat <- sum(w * e) # HH total of residuals
  v_lin <- (1 / N^2) *
    (1 / (n * (n - 1))) *
    sum((e / z - Te_hat)^2)
  
  # Confidence intervals
  z  <- qnorm(1 - (1 - conf)/2)
  se <- sqrt(v_lin)
  ci <- c(mu_GR - z * se,
          mu_GR + z * se)
  
  return(list(point = mu_GR,
              var = v_lin,
              se = se,
              ci = ci,
              level = conf))
}

################################################################################
# 2. Estimating the variance and confidence interval for Poisson Sampling under linearization
linVarGR_poisson <- function(sdata, mux, N, conf = 0.95) {
  y  <- sdata$y
  x  <- sdata$x
  pi <- sdata$pis # first-order inclusion probs
  d  <- 1 / pi    # HT weights
  n  <- length(y)
  
  # Point estimate
  muY_HT  <- sum(d * y) / N
  muX_HT  <- sum(d * x) / N
  mu_GR   <- (muY_HT / muX_HT) * mux
  Rhat    <- muY_HT / muX_HT
  
  # Linearised residuals
  e_hat <- y - Rhat * x
  
  # Poisson HT variance
  v_lin <- (1 / N^2) * sum((1 - pi) / pi^2 * e_hat^2)
  
  # Confidence intervals
  z   <- qnorm(1 - (1 - conf)/2)
  se  <- sqrt(v_lin)
  ci  <- c(mu_GR - z * se,
           mu_GR + z * se)
  
  return(list(point = mu_GR,
              var = v_lin,
              se = se,
              ci = ci,
              level = conf))
}

################################################################################
# Use the functions to estimate the variance for each PPS sampling method
mux <- mean(pdata$x)

# 1. For radomized systematic PPS sampling
vlinear_syspps  <- linVarGR_HH(sdata_syspps,
                              mux = mux,
                              N   = 5000)
vlinear_syspps
# 2. For Poisson Sampling
vlinear_poisson <- linVarGR_poisson(sdata_poisson,
                             mux = mux,
                             N   = 5000)
vlinear_poisson
# 3. For PPS sampling with replacement
vlinear_wrpps  <- linVarGR_HH(sdata_wrpps,
                          mux = mux,
                          N   = 5000)
vlinear_wrpps
