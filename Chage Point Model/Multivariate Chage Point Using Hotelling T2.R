require(plyr)
require(MASS)
require(methods)

## HotellingT2- return Hotelling T2 value for given Y and W matrix

HotellingT2 <- function(y, w){
  y = matrix(y, nrow = 1)
  w = as.matrix(w)
  y %*% solve(w) %*% t(y)
}


## T2ChangePoint - return two values (1) Maximum T2 and (2) the value of k

T2ChangePoint <- function(x, keepobs = NULL){
  x = as.matrix(x)
  #print(x)
  N = nrow(x)
  ncol = ncol(x)
  if (ncol > 1) {ccm = getFunction(colMeans)
  }else {ccm = get0("mean", as.environment("package:base"))}
  
  if (is.null(keepobs)){
    keepobs = ceiling(N*0.1)
    if (keepobs < ncol) keepobs = ncol + 1
  }
  
  aa <- function(z, n, nobs){
    x = z[1:n,]
    y = z[(n+1):nobs,]
    wk = (var(x)*n + var(y)*(nobs - n)) / (nobs - 2)
    yk = sqrt(n *(nobs - n) / nobs) * (ccm(x) - ccm(y))
    T2k = HotellingT2(yk, wk)
    return(T2k)
  }
  
  ntest = c((keepobs+1):(N-keepobs))
  T2 = sapply(ntest, aa, z = x, nobs = N)
  maxk = which.max(T2)
  return(c(T2Max = T2[maxk], Change_Point = maxk +keepobs, T2 = T2))
}

## Repeat Max T^2 for N (SimNum) times and return as a Data Frame containing two columns (T2Max and Change_Point (position))

set.seed(2345)
alpha = 0.05
p = 1
SimNum = 100
n =8
mu = rep(0, p)
sigma =  diag(x = 1, nrow = p, ncol = p)
a = rdply(SimNum, T2ChangePoint (mvrnorm(n = n, mu = mu, Sigma = sigma), keepobs = p+1))
b = a[,4:ncol(a)]

quant = list(c(NULL, NULL))
for (i in 1:ncol(b)){
  quant[[i]] = c(as.integer(p+1+i), quantile(b[,i], probs = 0.95, names = FALSE, type = 1))
  b = b[-which(unlist(b[,i]) >= quant[[i]][2]),]
}
bb = t(as.data.frame(quant))
row.names(bb) = NULL


