
#illustration of theorem 5.6. for d=4

d <- 4
n <- 60
alpha <- 1/3
K <- floor(alpha*n)

unitsphere <- rnorm(n*(d-1));dim(unitsphere) <- c(n,d-1)
for(k in seq_len(n)){unitsphere[k,] <- unitsphere[k,]/pracma::Norm(unitsphere[k,])}

unitsphere1 <- cbind(unitsphere[(1:K),],0)
unitsphere2 <- cbind(unitsphere[-(1:K),],1)

X <- rbind(0,c(rep(0,d-1),1),unitsphere1,unitsphere2)
n <- n+2
T <-  ddalpha::depth.halfspace(X,X,exact=TRUE)*n

Tv <- rep(0,n)
basic_model <- Tverberg_depth_basic_model4d(X)
for(k in (1:n)){
  model <- update_Tverberg_depth_model4d(basic_model,X,k)
  Tv[k] <- gurobi::gurobi(model)$objval
    print(k)
	plot(T[(1:k)],Tv[(1:k)])
  }
  
  plot(T,Tv)



#illustration of theorem 5.6. for d=3
set.seed(1234567)
while(TRUE){
d <- 4
n <- 90
alpha <- 1/2
K <- floor(alpha*n)

unitsphere <- rnorm(n*(d-1));dim(unitsphere) <- c(n,d-1)
for(k in seq_len(n)){unitsphere[k,] <- unitsphere[k,]/pracma::Norm(unitsphere[k,])}
eps <- 10^-5
unitsphere1 <- cbind(unitsphere[(1:K),],-eps*runif(K))
unitsphere2 <- cbind(unitsphere[-(1:K),],1+eps*runif(n-K))

#eps <- 10^-12
X <- rbind(unitsphere1,unitsphere2)
#alpha <-  runif(nrow(X));alpha[1] <- 30*alpha[1]; alpha <- alpha/sum(alpha)
#test_point <-alpha%*%X

#X <- rbind(test_point,X)
#n <- nrow(X)
T <-  ddalpha::depth.halfspace(rep(0,d),X,exact=TRUE)*n
if(T>=K/2-5){break}
}

n_rep=10000
T <- Tv <- rep(0,n_rep)
X <- rbind(0,unitsphere1,unitsphere2)
basic_model <- Tverberg_depth_basic_model4d(X)

for(k in (1:n_rep)){
    

X <- rbind(0,unitsphere1,unitsphere2)
alpha <-  runif(nrow(X));alpha[1] <- 1000*abs(rcauchy(1))*alpha[1]; alpha <- alpha/sum(alpha)
test_point <-alpha%*%X	
X[1,] <- test_point +rnorm(d,sd=sample(c(0.1,0.25,0.05),size=1))
if(k==1){X[1,] <- 0}
  model <- update_Tverberg_depth_model4d(basic_model,X,1,exclude_point_index=TRUE)
  temp <-gurobi::gurobi(model)
  Tv[k] <- temp$objval
  T[k] <- ddalpha::depth.halfspace(X[1,],X[-1,],exact=TRUE)*n
    print(k)
	plot(T[(1:k)],Tv[(1:k)])
  }
  
  plot(T,Tv)



#plotting

# 1. Load the ggplot2 library
library(ggplot2)

# 2. Create a sample data frame
depths <- t(unique(cbind(T,Tv)))
 
df <- data.frame(Tukey_depth=depths[1,],Tverberg_depth=depths[2,])
 


# 3. Generate the scatter plot with custom axis text
ggplot(df, aes(x = Tukey_depth, y = Tverberg_depth)) +
  geom_point(size = 3, color = "blue") + # Draws the points
  
  # Gitterlinien für X auf Abstand 1 setzen
  #scale_x_continuous(breaks = scales::breaks_width(1)) +
  # Gitterlinien für Y auf Abstand 1 setzen
  #scale_y_continuous(breaks = scales::breaks_width(1)) + 
  labs(
    x = "Tukey depth",     # Custom X-axis label
    y = "Tverberg depth",     # Custom Y-axis label
  )
#  title = "Main Plot Title (Optional)" # Optional title
 # )


n <- 60
K <- 20

eps <- delta <- jota <- 1/5
delta <- .1
set.seed(1234567)


angles1 <- seq(0,2*pi,length.out=K) +jota*runif(K,min=-1,max=+1)
angles2 <- seq(0,2*pi,length.out=n-K) +jota*runif(n-K,min=-1,max=+1)

norms1 <- runif(K,min=1-delta,max=1+delta)
norms2 <- runif(n-K,min=1-delta,max=1+delta)


X1 <- cbind(norms1*cos(angles1),norms1*sin(angles1),-eps*runif(K,min=-1,max=+1))
X2 <- cbind(norms2*cos(angles2),norms2*sin(angles2),1+eps*runif(n-K,min=-1,max=+1))
X <- rbind(0,c(0,0,1),X1,X2)
X <- rbind(X,colMeans(X))
n <- n+3
T <- Tv <- rep(0,n)
T <-  ddalpha::depth.halfspace(X,X,exact=TRUE)*n
print(sort(unique(T)))
basic_model <- Tverberg_depth_basic_model3d(X)
for(k in (1:n)){
  model <- update_Tverberg_depth_model3d(basic_model,X,k)
  Tv[k] <- gurobi::gurobi(model)$objval
    print(k)
	plot(T[(1:k)],Tv[(1:k)])
  }
  
  plot(T,Tv)