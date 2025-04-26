# Generate the population data
set.seed(854)
N <- 5000
n <- 100

x <- runif(N)
z <- 0.5 + rexp(N)

y <- 1 + 2 * x + 2.5 * z + rnorm(N)
pdata <- cbind(y, x, z)

################################################################################
# 1. Randomized Systematic PPS Sampling Method
syspps <- function(x, n){
  N = length(x) # The population size
  U = sample(N,N) # Sample the index without replacement
  xx = x[U] # Find the corresponding size values
  z = rep(0,N) # Initial the size variable
  for(i in 1:N) z[i] = n * sum(xx[1:i]) / sum(x) # The grid of b
  r = runif(1)
  s = numeric()
  for(i in 1:N){
    if(z[i] >= r){
      s = c(s, U[i])
      r = r + 1
    }
  }
  return(s[order(s)])
}
set.seed(854)
sam = syspps(z, n) # The n units selected
ys = y[sam] # The y values in S, same as ys=y[sam]
piN = n * z/sum(z) # The inclusion probabilities
sum(piN) # Check: sum(pi) = n
pis = piN[sam] # pi_i for i in S

# We normalize the z here for the approximate of HT using HH later
zs <- z/sum(z)
zs <- zs[sam]

sdata_syspps = pdata[sam,] # The sample data matrix
sdata_syspps <- cbind(sdata_syspps, pis, zs)

################################################################################
# 2. PPS sampling with replacement
set.seed(854)
z <- z/sum(z) # normalized size variable
sam <- sample(N, n, replace = T, prob = z)
ys <- y[sam] # Values of y in the sample
zs <- z[sam] # Values of the size variable

sdata_wrpps <- pdata[sam,]
sdata_wrpps <- cbind(sdata_wrpps, zs)
sdata_wrpps <- sdata_wrpps[, !(colnames(sdata_wrpps) == "z")]

################################################################################
# 3. Poisson sampling method
set.seed(854)
piN <- n * z / sum(z) # Inclusion probabilities
sam <- numeric()
for(i in 1:N){
  r = runif(1)
  if(r <= piN[i]) sam = c(sam,i)
}
ys = y[sam]
pis = piN[sam]

sdata_poisson = pdata[sam,]
sdata_poisson <- cbind(sdata_poisson, pis)
################################################################################
# Optional: Export the population and the sample data
# write.csv(pdata, "./data/pdata.csv", row.names = FALSE)
# write.csv(sdata_syspps, "./data/sdata_syspps.csv", row.names = FALSE)
# write.csv(sdata_wrpps, "./data/sdata_wrpps.csv", row.names = FALSE)
# write.csv(sdata_poisson, "./data/sdata_poisson.csv", row.names = FALSE)

