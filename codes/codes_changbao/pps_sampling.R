#############################################################
# STAT 854/454 Final Project: Three PPS Sampling Procedures #
#                                                           #
# 1. Randomized Systematic PPS Sampling Method              #
# 2. PPS Sampling With Replacement                          #
# 3. Poisson Sampling                                       #
#                                                           #
# Sample code by Changbao Wu, April 5, 2020                 #
#############################################################

#------------------------------
N=1000
n=50
#------------------------------
# Create the population
#------------------------------
x1 = rbinom(N, 1, 0.5) #The Bernoulli RV
x2 = runif(N) #The Uniform[0,1]
z = 0.5 + rexp(N) #The standard exponential distribution + 0.5
#-----------------------------
# The z variable will be used as the size variable
#-----------------------------
# A linear model for y
#-----------------------------
y=1+x1+2*x2+2.5*z+rnorm(N)
#-----------------------------
# The population data matrix: Nx4
#-----------------------------
pdata=cbind(y,x1,x2,z)
#==========================================
# Randomized Systematic PPS Sampling Method
# x: the size variable (i.e., the z variable)
# n: sample size
#==========================================
syspps=function(x,n){
  ##
  ##Population is first randomized!
  ##
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
#==================================================
# 1) Take a sample by the randomized systematic PPS
#==================================================
sam = syspps(z,n) #The n units selected
sdata = pdata[sam,] #The sample data matrix: nx4
ys = sdata[,1] #The y values in S, same as ys=y[sam]
#--------------------------------------------------
piN = n * z/sum(z) # The inclusion probabilities
sum(piN) # Check: sum(pi) = n
pis = piN[sam] # pi_i for i in S
#--------------------------------------------------
# HT estimators for T_y and mu_y
#--------------------------------------------------
TyHT = sum(ys/pis)
muyHT = sum(ys/pis)/N
#==================================================
# 2) PPS sampling with replacement
#==================================================
z=z/sum(z) #normalized size variable
sam=sample(N,n,replace=T,prob=z)
ys=y[sam] #Values of y in the sample
zs=z[sam] #Values of the size variable
TyHH=mean(ys/zs) #The HH estimator for Ty
#==================================================
# 3) Poisson sampling method
#==================================================
piN=n*z/sum(z) #Inclusion probabilities
sam=numeric()
for(i in 1:N){
  r=runif(1)
  if(r<=piN[i]) sam=c(sam,i)
  
}
#-------------------
ys=y[sam]
pis=piN
TyHT=sum(ys/pis)
muyHT=sum(ys/pis)/N
###################################################

