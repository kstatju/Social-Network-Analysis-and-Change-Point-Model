datapath = 'C:/Users/ka746940/Desktop/UCF/STA 6908 - Edgard Maboudou/Data/'
df = read.table(paste(datapath, 'pb2.txt', sep = ''))
x = df[df$V1 == 1, 2:5]
y = df[df$V1 == 2, 2:5]


kernel11 <- function(x1, y = NULL){
  if (!is.null(y)){  x1 = x1-y}
  if(all(x1 == 0)){return(0)}
  K<- sum(x1 / base::norm(as.matrix(x1),type = 'f'))
  return(K)
}

kcalculator <- function(X, Y = NULL){
  if (is.null(Y)){
  N<-dim(X)[1]
  K<-matrix(0,N,N)
    for(i in 1:N){
        K[i,]<-apply(sweep(X, 2, unlist(X[i,])), 1, rbf_kernel)
    }
  return(K)
  }else{
    N1<-dim(X)[1]
    N2 = dim(Y)[1]
    K<-matrix(0,N1,N2)
    for(i in 1:N2){
      K[,i]<-apply(sweep(X, 2, unlist(Y[i,])), 1, rbf_kernel)
    }
    return(K)
    
  }
}

x <- matrix(runif(300),100)
y <- matrix(runif(300)+10,100)
a = kmmd(x, y)
xmatrix(a)
H0(kmmd(x, y))

rbf <- rbf_kernel(x-y)
aa = kernelMatrix(kernel11, as.matrix(x))

k = kcalculator(x, y)
library(kernlab)
kmmd(x, y, Kxy = "kernel11")

