





n <- 50 #99
Tv <- T <- rep(0,1000)
set.seed(1234567)
#59
for(k in (1:1000)){
 while(TRUE){
   X <- rnorm(n*4)^2;dim(X)=c(n,4)
   X[1,] <- rnorm(4,sd=.2)
   T[k] <- ddalpha::depth.halfspace(X[1,],X,exact=TRUE)
  if(T[k]*n  < n/5-1 & T[k]>2/n){break}
 }
  
#X_rational <- safe_double_to_rational(as.vector(X));dim(X_rational) <- dim(X)
basic_model <- Tverberg_depth_basic_model4d(X)
model <- update_Tverberg_depth_model4d(basic_model,X,1)
Tv[k] <- gurobi::gurobi(model)$objval

par(mfrow=c(1,2))
plot(X)
points(X[c(1,1),],col="red",pch=16)
plot(T[(1:k)]*n,Tv[(1:k)])
lines((1:n),(1:n))
lines((1:n),(1:n)-1,col="grey")
lines((1:n),(1:n)-2,col="grey")
lines((1:n),(1:n)-3,col="grey")
lines((1:n),(1:n)-4,col="grey")
lines((1:n),(1:n)-5,col="grey")


if(T[k]*n-Tv[k] >=3){break}
print(k)
}

X <- runif(n*3,min=0,max=1);dim(X) <- c(n,3)
X[seq_len(n/2),3] <- X[seq_len(n/2),2]+X[seq_len(n/2),1]+rnorm(n/2,sd=0.1)#cbind(X,X[,1]^2+X[,2]^2+rnorm(n,sd=0.5))
#X[,2] <- 0.1*X[,2]+X[,1]^2
T2 <- Tv2 <- rep(0,n)

o <- order(-ddalpha::depth.halfspace(X,X,exact=TRUE))
X <- X[o,]


n <- 35
set.seed(1234567)
X <- runif(n*3);dim(X) <- c(n,3)
for(k in seq_len(nrow(X))){X[k,] <- X[k,]/sum(X[k,])}
X[(1:n/2),] <- -X[(1:n/2),]
#X <- X[,(1:2)]
basic_model <- Tverberg_depth_basic_model3d(X)
for(k in seq_len(n)){
  model <- update_Tverberg_depth_model3d(basic_model,X,k)
  Tv2[k] <- gurobi::gurobi(model)$objval
  T2[k] <-  ddalpha::depth.halfspace(X[k,],X,exact=TRUE)
  plot(T2[(1:k)]*n,Tv2[(1:k)])
  lines((1:n),(1:n))
  lines((1:n),(1:n)-1,col="grey")
  lines((1:n),(1:n)-2,col="grey")
  lines((1:n),(1:n)-3,col="grey")
  lines((1:n),(1:n)-4,col="grey")
  lines((1:n),(1:n)-5,col="grey")
  
  
  
  print(k)
}


