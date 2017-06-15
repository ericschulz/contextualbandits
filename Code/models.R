#Models
##############################################################################################################
#KERNELS
##############################################################################################################
#Radial Basis Kernel
#X1, X2 is  x and x` in the kernel function k(x,x')
#theta is a parameter vector (lambda, signal variance, noise variance), where lambda is identical for all dimensions
rbf <- function(X1,X2, theta){
  #transfer to matrices
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  #check dimensions
  if(ncol(X1) != ncol(X2)){
    stop("X1 and X2 must contain input values of the same dimension.")
  } else if(!all(theta>=0)){
    stop("All parameters must be >= 0.")
  }
  #get dimensions
  N1 <- nrow(X1)
  N2 <- nrow(X2)
  d <- ncol(X1)
  #initialize sigma
  sigma <-  matrix(rep(0, N1*N2),nrow=N1)
  #Length scale parameter
  lambda <- theta[1]
  sn <- theta[2]
  #I've fixed sf to 1 now
  #loop through
  for(i in 1:d){
    #x-diff
    xdiff <- (outer(X1[,i],X2[,i],function(x,y) x - y)/lambda)^2
    sigma <- sigma + xdiff
  }
  #RBF function
  if(identical(X1,X2)){
    id <- diag(rep(1,N1))
    sigma.final <- exp(-0.5*sigma) + sn*id
  } else {
    sigma.final <- exp(-0.5*sigma)
  }
  #return final covariance matrix
  return(sigma.final)
}

##############################################################################################################
#MATRIX INVERSION
##############################################################################################################
#calculate inverse of the cov-function using sigular value decomposition
cov.inverse.svd <- function(X, tol = sqrt(.Machine$double.eps)){
  # Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  #singular value decomposition
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  #inverse
  K.inv <- structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
  #logarithm of determinant.
  log.K.det <- sum(log(s$d))
  #return inverse plus log determinant
  return(list(Inv = K.inv,lDet = log.K.det))
}

#calculate inverse of the cov-function using Cholesky
cov.inverse.chol <- function(X){
  #cholseky decomposition
  R <- chol(X)
  #complex conjugate
  Rt <- Conj(t(R))
  #invert
  R.inv <- solve(R)
  #invert
  Rt.inv <- solve(Rt)
  #multiply matrices
  X.inv <- R.inv %*% Rt.inv
  #log determinant
  log.X.det <- 2*sum(log(diag(R))) 
  #return both
  return(list(Inv = X.inv, lDet = log.X.det))
}

##############################################################################################################
#GAUSSIAN PROCESS
##############################################################################################################
#Gaussian Process function
#X.test: matrix for predcitions
#theta: vector of hyper-parameters (lambda, Sf, Sn)
#X; matrix of observations
#y: vector of observed outcomes
gp<- function(X, Y, X.test, theta){
  #if no observations
  if (dim(X)[1]==0){
    prediction<-as.data.frame(t(c(0,theta[3]))) #prior is just N(0, noise variance)
    colnames(prediction) <- c("mu", "sig")
  }else{ #General case
    #make it a matrix
    Xstar <- as.matrix(X.test)
    #dimensions
    d <- ncol(as.matrix(X))
    #calculate capital K
    K <- rbf(X,X,theta) 
    #Check if matrix is positive semi-definite
    if (is.positive.definite(K)){
      KK <- cov.inverse.chol(K) #use Cholesky
    } else {
      KK <- cov.inverse.svd(K) #use SVD
    }
    #times y
    Ky <- KK$Inv %*% Y
    #apply the kernel
    result <- apply(Xstar, 1, function(x){
      XX <- matrix(x,nrow=1) 
      Kstar <- rbf(X, XX, theta)
      Kstarstar <- rbf(XX,XX,theta)
      #get mean vector
      mu <- t(Kstar) %*% Ky
      #get covariance
      cv <- Kstarstar - (t(Kstar) %*% KK$Inv %*% Kstar) #BUG: sometimes cv<0, leading to NaN when we return sqrt(cv)
      #DEBUG
      if (cv<0){ cv <- 0} #TEMPORARY SOLUTION: MANUALLY SET CV TO BE 0
      #return means and standard deviations
      return(c(mu, sqrt(cv)))
    })
    #as a data frame with names mu and sig
    prediction <- as.data.frame(t(result))
    colnames(prediction) <- c("mu", "sig")
  }
  #return it
  return(prediction)
}
class(gp)<- c(class(gp), "GP")

library(MCMCpack)

lin<- function(X, Y, X.test){
  X<-matrix(X, ncol=3)
  Y<-as.matrix(Y, ncol=1)
  X1<-matrix(runif(30),nrow=10, ncol=3)
  Y1<-matrix(runif(10), nrow=10, ncol=1)
  m <- MCMCregress(rbind(Y1,Y)~rbind(X1,X), burnin=100, MCMC=100)
  b <- m[seq(1,100,1),1]+m[seq(1,100,1),2]*X.test[,1]+m[seq(1,100,1),3]*X.test[,2]+m[seq(1,100,1),4]*X.test[,3]
  prediction<-data.frame(mu=mean(b), sig=sd(b))
  #return it
  return(prediction)
}
class(lin)<- c(class(lin), "LIN")

#Faster version of the kalman filter that computes a single update
kalman <- function(theta=c(1,1), chosen, y, x, prevPost=NULL){
  #parameters
  mu0 <- 0 #prior mean
  var0 <- 5 #prior variance
  vart <- theta[1] #innovation variance
  vare <- theta[2] #error varriance
  if (is.null(prevPost)){#if no posterior prior, assume it is the first observation
    predictions <- data.frame(mu=rep(mu0,4), sig=rep(var0,4))
  }else{#if previous posterior is provided, update
    predictions <- prevPost
  }
  #Kalman gain
  kGain <- (predictions$sig[chosen] + vart^2) / (predictions$sig[chosen] + vart^2 + vare^2)
  #update mean
  predictions$mu[chosen] <- predictions$mu[chosen] + (kGain * (y-predictions$mu[chosen]))
  #update variance
  predictions$sig <- predictions$sig + vart^2
  predictions$sig[chosen] <- predictions$sig[chosen] * (1 - kGain)
  #return output
  return(predictions)
}
class(kalman)<- c(class(kalman), "KAL")

#Faster version of the kalman filter that computes a single update
meanbayes <- function(theta=c(1), chosen, y, x, prevPost=NULL){
  #parameters
  mu0 <- 0 #prior mean
  var0 <- 5 #prior variance
  vart <- 0 #innovation variance
  vare <- theta[1] #error varriance
  if (is.null(prevPost)){#if no posterior prior, assume it is the first observation
    predictions <- data.frame(mu=rep(mu0,4), sig=rep(var0,4))
  }else{#if previous posterior is provided, update
    predictions <- prevPost
  }
  #Kalman gain
  kGain <- (predictions$sig[chosen] + vart^2) / (predictions$sig[chosen] + vart^2 + vare^2)
  #update mean
  predictions$mu[chosen] <- predictions$mu[chosen] + (kGain * (y-predictions$mu[chosen]))
  #update variance
  predictions$sig <- predictions$sig + vart^2
  predictions$sig[chosen] <- predictions$sig[chosen] * (1 - kGain)
  #return output
  return(predictions)
}
class(meanbayes)<- c(class(meanbayes), "MEB")

