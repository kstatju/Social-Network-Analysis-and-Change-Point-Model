require(kernlab)
require(plyr)
require(MASS)
x <- matrix(rnorm(100, 10, 5),50)
y <- matrix(rnorm(200, 20, 5),100)
z <- matrix(rnorm(9, 20, 5),3, 3)

colSums(sweep(x, 2,colMeans(x)))
############################################################

HotellingT2 <- function(y, w){
  y = matrix(y, nrow = 1)
  w = as.matrix(w)
  y %*% solve(w) %*% t(y)
}



cpkmmd <- function(x, keepobs = NULL){
  x = as.matrix(x)
  N = nrow(x)
  if (is.null(keepobs)){
    keepobs = ceiling(N*0.1)
  }
  
  aa <- function(z, n, nobs){
    x = z[1:n,]
    y = z[(n+1):nobs,]
    wk = (var(x)*n + var(y)*(nobs - n)) / (nobs - 2)
    yk = sqrt(n *(nobs - n) / nobs) * (colMeans(x) - colMeans(y))
    T2k = HotellingT2(yk, wk)
    return(T2k)
  }
  
  
  ntest = c((keepobs+1):(N-keepobs))
  stat = sapply(ntest, aa, z = x, nobs = N)
  mk = which.max(stat)
  return(c(T2Max = stat[mk], change_Point = mk+keepobs))
}

a = rdply(100, cpkmmd(mvrnorm(n = 20, mu = c(11,15,32), 
                              Sigma =  matrix(c(1, .6, .7, .6, 1, .75, .7, .75,1), 
                                              nrow = 3, byrow = TRUE))))




# KMMD function 
# only return KMMD value

kmmdfunconly <- function(Kxx,Kyy, Kxy)
{
  
  m <- dim(Kxx)[1]
  n <- dim(Kyy)[1]
  
  sumKxx <- sum(Kxx)
  
  sumKyy <- sum(Kyy)
 
  sumKxy <- sum(Kxy)
  
  mmd1 <- sqrt(max(0,sumKxx/(m*m) + sumKyy/(n*n) - 2/m/n* sumKxy))

  return(list(mmd=mmd1))
}


############################################################

# KMMD function 
# return KMMD with 1st and 3rd order 

kmmdfunc <- function(Kxx,Kyy, Kxy, alpha = 0.05)
{
  
  m <- dim(Kxx)[1]
  n <- dim(Kyy)[1]
  
  N <- max(m,n)
  M <- min(m,n)
  
  sumKxx <- sum(Kxx)
  
  if(m!=n)
    sumKxxM <- sum(Kxx[1:M,1:M])
  else
    sumKxxM <- sumKxx
  
  dgxx <- diag(Kxx)
  
  sumKxxnd <- sumKxx - sum(dgxx)
  R <- max(dgxx)
  RM <- max(dgxx[1:M])
  hu <- colSums(Kxx[1:M,1:M]) - dgxx[1:M]
  
  sumKyy <- sum(Kyy)
  if(m!=n)
    sumKyyM <- sum(Kyy[1:M,1:M])
  else
    sumKyyM <- sumKyy
  
  dgyy <- diag(Kyy)
  
  sumKyynd <- sum(Kyy) - sum(dgyy)
  R <- max(R,dgyy)
  RM <- max(RM,dgyy[1:M]) # RM instead of R in original
  hu <- hu + colSums(Kyy[1:M,1:M]) - dgyy[1:M]
  
  sumKxy <- sum(Kxy)
  if (m!=n)
    sumKxyM <- sum(Kxy[1:M,1:M])
  else
    sumKxyM <- sumKxy
  
  dg <- diag(Kxy) # up to M only
  hu <- hu - colSums(Kxy[1:M,1:M]) - colSums(t(Kxy[1:M,1:M])) + 2*dg # one sided sum
  
  mmd1 <- sqrt(max(0,sumKxx/(m*m) + sumKyy/(n*n) - 2/m/n* sumKxy))
  mmd3 <- sum(hu)/M/(M-1)
  D1 <- 2*sqrt(RM/M)+sqrt(log(1/alpha)*4*RM/M)
  
  return(list(mmd1=mmd1,mmd3=mmd3,D1=D1))
}


## Find Change point location by calculating maximum KMMD for data

cpkmmd <- function(x, keepobs = NULL){
  x = as.matrix(x)
  N = nrow(x)
  if (is.null(keepobs)){
    keepobs = ceiling(N*0.1)
  }
  aa <- function(z, n, rbf = rbf){
    nobs = nrow(z)
    ncoll = ncol(z)
    x = z[1:n,1:ncoll]
    y = z[(n+1):nobs,1:ncoll]
    stat <- kmmdfunconly(Kxx=kernelMatrix(rbf, x), Kyy =kernelMatrix(rbf, y), Kxy=kernelMatrix(rbf, x,y))
    return(stat$mmd)
  }
  
  sig.opt <- sigest(x, scaled = FALSE)[2]
  rbf <- rbfdot(sigma = sig.opt)
  
  ntest = c((keepobs+1):(N-keepobs))
  stat = sapply(ntest, aa, z = x, rbf = rbf)
  mk = which.max(stat)
  return(list(stat = stat, change.Point = mk+keepobs))
}

bb = cpkmmd(rbind(x,y))

newobscpkmmd <- function(newobs, x, y){
  ncoll = ncol(x)
  newobs = matrix(newobs, ncol = ncoll)
  x = as.matrix(x)
  y = as.matrix(y)
  nnew = nrow(newobs)
  ns = c(1:nnew)
  sig.opt <- sigest(rbind(x,y), scaled = FALSE)[2]
  rbf <- rbfdot(sigma = sig.opt)
  
  aa <- function(i, z, x, y, rbf = rbf){
    a = z[i,]
    s1 = kmmdfunconly(Kxx=kernelMatrix(rbf, rbind(x,a)), Kyy =kernelMatrix(rbf, y), Kxy=kernelMatrix(rbf, rbind(x,a),y))
    s2 = kmmdfunconly(Kxx=kernelMatrix(rbf, x), Kyy =kernelMatrix(rbf, rbind(y,a)), Kxy=kernelMatrix(rbf, x,rbind(y,a)))
    if (s1$mmd > s2$mmd) c = 1
    else c = 2
    return(c)
  }
  
  co = sapply(ns, aa, z = newobs, x = x, y = y, rbf = rbf)
  
  return(co)
  
  
}



newobscpkmmd(z[1,], x,y)

