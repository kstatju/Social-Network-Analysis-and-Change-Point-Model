#
# Compute kernel MMD test
f.kmmd<-function(k1,k2,k3){
  
  kmmd2<-0
  m<-nrow(k1)
  for(i in 1:m){
    t<-0
    for(j in 1:m){
      if( i != j){
        
        t <- t + k1[i,j]+k2[i,j]-k3[i,j]-k3[j,i]
        
      }# end if
      
    } #end loop j
    
    kmmd2<-kmmd2+t/(m*(m-1))
    
  } # end loop i
  
  return(kmmd2)
  
}


#Kernel function

kernel11 <- function(x1, y = NULL){
  if (!is.null(y)){  x1 = x1-y}
  if(all(x1 == 0)){return(0)}
  K<- sum(x1 / base::norm(as.matrix(x1),type = 'f'))
  return(K)
}

# Permutation test
PermTest.knl1 <- function(x, y, R=499, testfun=f.kmmd) {
  z <- rbind(x, y)  # pooled sample
  # Create kernel matrix
  sig.opt <- sigest(z, scaled = FALSE)[2]
  #rbf <- rbfdot(sigma = sig.opt)
  rbf <- getFunction("kernel11")
  
  myfun <- function(a, b,c) suppressWarnings(unname(testfun(a, b,c)))
  #set.seed(123) 
  DoIt <- function() {
    i <- sample(nrow(z), nrow(x))
    myfun( a=kernelMatrix(rbf, z[i,]), b =kernelMatrix(rbf, z[-i,]), c=kernelMatrix(rbf, z[i,],z[-i,]))
  }
  pstats <- replicate(R, DoIt())
  stat <- myfun(a=kernelMatrix(rbf, x), b =kernelMatrix(rbf, y), c=kernelMatrix(rbf, x,y))
  hist(pstats, col="azure", main="Empirical distribution under H0")
  abline(v=stat, col="red", lty=2)
  p.v <-mean(c(stat, pstats) >= stat)
  res <- list(p.v=p.v, sig.opt=sig.opt, R=R, stat=stat, pstats=pstats, sum=sum(c(stat, pstats) >= stat))
  return(res)
}







