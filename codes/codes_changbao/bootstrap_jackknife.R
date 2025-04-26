####################################################
## R Example for Bootstrap and Jackknife          ##
## Variance Estimation for the Ratio Estimator    ##
##                                                ##
## By Changbao Wu, February 20, 2009              ##
####################################################

#=================================================
# R function for bootstrapping the ratio estimator
#=================================================

VBoot = function(xs, ys, mux, n, B) {
  V = rep(0, B)
  for (i in 1:B) {
    bsam = sample(n, n, replace = TRUE)
    bx = xs[bsam]
    by = ys[bsam]
    V[i] = (mean(by) / mean(bx)) * mux
  }
  VB = (B - 1) * var(V) / B
  return(VB)
}

#================================================
# R function for jackknifing the ratio estimator
#================================================

VJack = function(xs, ys, mux, n) {
  V = rep(0, n)
  for (i in 1:n) {
    V[i] = (mean(ys[-i]) / mean(xs[-i])) * mux
  }
  VJ = ((n - 1)^2 / n) * var(V)
  return(VJ)
}

##=============================================
## Create the Finite Population for Simulation
##=============================================

N = 10000        # Popu size
n = 160          # Sample size

#--------------------------------------
# The Finite Population: y = 4 + x + e
#--------------------------------------

x = 1 + rexp(N)
e = rnorm(N)
y = 4 + x + e

#----------------------
# Population Parameters
#----------------------

muy = mean(y)
mux = mean(x)

#=========================
# Take a sample by SRSWOR
#=========================

sam = sample(N, n)
ys = y[sam]
xs = x[sam]

#-----------------------------
# Compute the ratio estimator
#-----------------------------

muyhat = (mean(ys) / mean(xs)) * mux

#---------------------------------------------
# Bootstrap and Jackknife Variance estimators
#---------------------------------------------

VB = VBoot(xs, ys, mux, n, 1000)
VJ = VJack(xs, ys, mux, n)

# Display results
muyhat
muy
VB
VJ
