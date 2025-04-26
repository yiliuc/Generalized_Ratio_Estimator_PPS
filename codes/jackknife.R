# Read the simulated data
pdata <- read.csv("data/pdata.csv")
sdata_syspps <- read.csv("data/sdata_syspps.csv")
sdata_poisson <- read.csv("data/sdata_poisson.csv")
sdata_wrpps <- read.csv("data/sdata_wrpps.csv")

################################################################################
# 1. The Delete-1 Jackknife variance estimation for H-T estimators
varJackGR <- function(sdata, mux, N, level = 0.95) {
  ys  <- sdata$y
  xs  <- sdata$x
  pis <- sdata$pis
  
  n    <- length(ys)
  di   <- 1 / pis
  Nhat <- sum(di)
  
  # Point estimate
  muY_HT <- sum(di * ys) / N
  muX_HT <- sum(di * xs) / N
  muGR   <- (muY_HT / muX_HT) * mux
  
  muGR_mi <- numeric(n)
  for (i in seq_len(n)) {
    gamma   <- Nhat / (Nhat - di[i])
    di_mi   <- di * gamma
    di_mi[i] <- 0
    muY_HT_mi <- sum(di_mi * ys) / sum(di_mi)
    muX_HT_mi <- sum(di_mi * xs) / sum(di_mi)
    muGR_mi[i] <- (muY_HT_mi / muX_HT_mi) * mux
  }
  
  # Jackknife variance
  mu_dot <- mean(muGR_mi)
  vJack  <- (n - 1) / n * sum((muGR_mi - mu_dot)^2)
  seJack <- sqrt(vJack)
  
  # Confidence interval
  alpha <- 1 - level
  z     <- qnorm(1 - alpha/2)
  ci    <- c(lower = muGR - z * seJack,
             upper = muGR + z * seJack)
  
  return(list(point = muGR,
              var = vJack,
              ci = ci,
              estimates = muGR_mi))
}

################################################################################
# 2. Delete-1 Jackknife estimator using HH estimators
varJackGR_wrpps <- function(sdata, mux, level = 0.95) {
  ys <- sdata$y
  xs <- sdata$x
  zs <- sdata$zs
  n  <- length(ys)
  
  # HH weights
  z_i <- zs / sum(zs)
  wi  <- 1 / (n * z_i)
  
  # Point estimate
  muY_HH <- sum(wi * ys) / sum(wi)
  muX_HH <- sum(wi * xs) / sum(wi)
  muY_GR <- (muY_HH / muX_HH) * mux
  
  muGR_mi <- numeric(n)
  for (i in seq_len(n)) {
    # Replicate weights
    wi <- wi
    wi_mi <- (n / (n - 1)) * wi
    wi_mi[i] <- 0
    
    muY_HH_mi <- sum(wi_mi * ys) / sum(wi_mi)
    muX_HH_mi <- sum(wi_mi * xs) / sum(wi_mi)
    muGR_mi[i] <- (muY_HH_mi / muX_HH_mi) * mux
  }
  
  # Jackknife variance
  mu_dot <- mean(muGR_mi)
  vJack  <- (n - 1) / n * sum((muGR_mi - mu_dot)^2)
  seJack <- sqrt(vJack)
  
  alpha  <- 1 - level
  zcrit  <- qnorm(1 - alpha/2)
  ci     <- c(lower = muY_GR - zcrit * seJack,
              upper = muY_GR + zcrit * seJack)
  
  return(list(point = muY_GR,
              variance = vJack,
              ci = ci,
              estimates = muGR_mi))
}

################################################################################
mux <- mean(pdata$x)
N   <- 5000

vJack_syspps <- varJackGR(sdata_syspps, mux, N)
vJack_syspps
vJack_poisson <- varJackGR(sdata_poisson, mux, N)
vJack_poisson

vJack_wrpps <- varJackGR_wrpps(sdata_wrpps, mux)
vJack_wrpps
