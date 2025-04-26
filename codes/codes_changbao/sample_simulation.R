#########################################################
# STAT 854/454 Final Project: An Example for Simulation #
#                                                      #
# 1. Bias and MSE of Two Point Estimators              #
# 2. Length and Coverage Probabilities of CIs          #
# 3. Samples are taken by SRSWOR                       #
#                                                      #
# By Changbao Wu, February 20, 2009                    #
#########################################################

#--------------------------
set.seed(7654321,kind=NULL)
#--------------------------

nsim=1000 # Number of simulation runs
N=10000 # Popu size
n=160 # Sample size
bias=c(0,0)
mse=c(0,0)
cp=matrix(0,2,3)
len=c(0,0)

#--------------------------------------
# The Finite Population: y=4+x+e
#--------------------------------------
x=1+rexp(N)
e=rnorm(N)
y=4+x+e

#----------------------
# Population Parameters
#----------------------
muy=mean(y) # popu mean
p=mean(y>=5) # popu proportion

#======================================
# Repeated simulation starts from here:
#======================================
for(m in 1:nsim){
  sam=sample(N,n) # Take a SRSWOR sample
  ys=y[sam] # This is the sample data
  muhat=mean(ys) # Estimator of muy
  phat=mean(ys>=5) # Estimator of p
  
  #-----------------
  # Calculate 95% CI
  #-----------------
  mu1=muhat-1.96*sqrt((1-n/N)*var(ys)/n)
  mu2=muhat+1.96*sqrt((1-n/N)*var(ys)/n)
  p1=phat-1.96*sqrt((1-n/N)*phat*(1-phat)/n)
  p2=phat+1.96*sqrt((1-n/N)*phat*(1-phat)/n)
  
  #--------------------
  # Calcu bias, mse etc
  #--------------------
  bias[1]=bias[1]+muhat-muy
  bias[2]=bias[2]+phat-p
  mse[1]=mse[1]+(muhat-muy)^2
  mse[2]=mse[2]+(phat-p)^2
  len[1]=len[1]+mu2-mu1
  len[2]=len[2]+p2-p1
  cp[1,1]=cp[1,1]+(muy<=mu1)
  cp[1,2]=cp[1,2]+(muy>mu1)*(muy<mu2)
  cp[1,3]=cp[1,3]+(muy>=mu2)
  cp[2,1]=cp[2,1]+(p<=p1)
  cp[2,2]=cp[2,2]+(p>p1)*(p<p2)
  cp[2,3]=cp[2,3]+(p>=p2)
}

#==================
# End of simulation
#==================
bias[1]=bias[1]/(muy*nsim)
bias[2]=bias[2]/(p*nsim)
mse=mse/nsim
len=len/nsim
cp=cp/nsim

bias
mse
len
cp
