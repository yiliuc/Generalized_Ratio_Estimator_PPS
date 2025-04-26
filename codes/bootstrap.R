# Read the generated data
pdata <- read.csv("data/pdata.csv")
sdata_syspps <- read.csv("data/sdata_syspps.csv")
sdata_poisson <- read.csv("data/sdata_poisson.csv")
sdata_wrpps <- read.csv("data/sdata_wrpps.csv")

################################################################################
# 1. Bootstrap function for randomized systematic PPS and Poisson sampling (use pi)
varBootGR <- function(sdata, mux, N, B = 1000, conf = 0.95){
  set.seed(3)
  ys <- sdata$y
  xs <- sdata$x
  pis <- sdata$pis
  n <- length(ys)
  
  # Hortitz-Thompson weights
  di  <- 1 / pis
  
  # Point estimate
  muY_HT <- sum(di * ys) / N
  muX_HT <- sum(di * xs) / N
  muY_GR <- (muY_HT / muX_HT) * mux
  
  muGR_star <- numeric(B)
  for (b in seq_len(B)) {
    bsam <- sample.int(n, n, replace = TRUE)
    yb <- ys[bsam]
    xb <- xs[bsam]
    db <- di[bsam]
    
    muY_HT_star <- sum(db * yb) / N
    muX_HT_star <- sum(db * xb) / N
    muGR_star[b] <- (muY_HT_star / muX_HT_star) * mux
  }
  
  # Compute the variance of generalized ratio estimator
  vB <- (B - 1) * var(muGR_star) / B
  
  # Confidence interval
  alpha <- 1 - conf
  ci    <- quantile(muGR_star,
                    probs = c(alpha/2, 1 - alpha/2),
                    names = FALSE)
  
  return(list(var = vB,
       ci  = ci,
       level = conf,
       point = muY_GR,
       estimates = muGR_star))
}

################################################################################
# 2. Bootstrap function for PPS sampling with replacement (use z)
varBootGR_wrpps <- function(sdata, mux, N, B = 1000, conf = 0.95){
  set.seed(854)
  ys <- sdata$y
  xs <- sdata$x
  zs <- sdata$zs
  n  <- length(ys)
  
  # Hansen–Hurwitz weights
  dHH <- 1 / (n * zs)
  
  # Point estimate
  Ty_HH <- sum(dHH * ys)
  Tx_HH <- sum(dHH * xs)
  muY_GR <- (Ty_HH / Tx_HH) * mux
  
  muGR_star <- numeric(B)
  for (b in seq_len(B)) {
    bsam <- sample.int(n, n, replace = TRUE)
    
    yb <- ys[bsam]
    xb <- xs[bsam]
    zb <- zs[bsam]
    dB <- 1 / (n * zb) 
    
    Ty_star <- sum(dB * yb)
    Tx_star <- sum(dB * xb)
    muGR_star[b] <- (Ty_star / Tx_star) * mux
  }
  
  # Compute the variance
  vB <- (B - 1) * var(muGR_star) / B
  
  # Confidence interval of mu_y
  alpha <- 1 - conf
  ci <- quantile(muGR_star,
                 probs = c(alpha/2, 1 - alpha/2),
                 names = FALSE)
  
  return(list(var = vB,
              ci = ci,
              level = conf,
              point = muY_GR,
              estimates = muGR_star))
}

################################################################################
mux <- mean(pdata$x)
N <- 5000
B <- 1000

vBoot_syspps <- varBootGR(sdata_syspps, mux, N, B)
vBoot_syspps
vBoot_poisson <- varBootGR(sdata_poisson, mux, N, B)
vBoot_poisson
vBoot_wrpps <- varBootGR_wrpps(sdata_wrpps, mux, N, B = 1000)
vBoot_wrpps
