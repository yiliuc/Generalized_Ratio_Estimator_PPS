piij <- function(x, s) {
  N <- length(x)
  n <- length(s)
  p <- matrix(0, n, n)
  
  for (k in 1:10000) {
    ss <- syspps(x, n)
    
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (min(abs(ss - s[i])) + min(abs(ss - s[j])) == 0) {
          p[i, j] <- p[i, j] + 1
        }
      }
    }
    if (k %% 1000 == 0)
      message(" iterations completed: ", k, "/", B)
  }
  
  p <- (p + t(p)) / 10000
  return(p)
}

results <- piij(x, sam)
