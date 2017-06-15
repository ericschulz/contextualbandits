#Gaussian Process Bandits
#Eric Schulz, Feb 2017, e.schulz@cs.ucl.ac.uk
#Schulz, Konstantinidis, Speekenbrink, 2017.
##############################################################################################################

##############################################################################################################
#PREAMPLE
##############################################################################################################

#House keeping
rm(list=ls())

#required packages packages
packages <- c('RCurl', 'devtools', 'matrixcalc', 'DEoptim', 'jsonlite')
#load them
lapply(packages, library, character.only = TRUE)

setwd("/home/hanshalbe/Desktop/cmabpred")

#get the data
source("getdata.R")
data<-getdata(3)


source("models.R")



##############################################################################################################
#GP BANDITS
##############################################################################################################
#GP multi-armed bandit function
#par: vector of parameters for softmax and GP-fitting
#X: matrix of observation participant saw
#chosen: arm participant chose
#y: outcome of chosen arm
#kernel: kernel used for fitting
#acquisition: acquisition function
gpbandit<-function(par, X, chosen, y, kernel, acquisition){
  #exponentiate parameters as par>0
  par<-exp(par)
  #gp parameters are the ones up to the penultimate
  pargp<-par[1:(length(par-1))]
  #last parameter is inverse temperature for softmax
  beta<-par[length(par)]
  #initialize ys
  y1<-y2<-y3<-y4<-c(0)
  #intialize Xs
  X1<-X2<-X3<-X4<-matrix(0, ncol=3, nrow=1)
  #make suer X is a matrix
  X<-as.matrix(X)
  #intiailize out frames
  out1<-out2<-out3<-out4<-data.frame(mu=numeric(), sig=numeric())
  #loop through observations
  for (i in 1:nrow(X)){
    #new observation
    Xnew<-matrix(X[i,], ncol=3)  
    #prediction for 1st arm based on new observations
    out1<-rbind(out1, gpr(X.test=Xnew, theta=pargp, X=X1, Y=y1, kernel=kernel))    
    #prediction for 2nd arm based on new observations
    out2<-rbind(out2, gpr(X.test=Xnew, theta=pargp, X=X2, Y=y2, kernel=kernel))    
    #prediction for 3rd arm based on new observations
    out3<-rbind(out3, gpr(X.test=Xnew, theta=pargp, X=X3, Y=y3, kernel=kernel))    
    #prediction for 4th arm based on new observations
    out4<-rbind(out4, gpr(X.test=Xnew, theta=pargp, X=X4, Y=y4, kernel=kernel))  
    #if participannt chose arm 1, append to 1-frames
    if (chosen[i]==1){
      X1<-rbind(X1, Xnew)
      y1<-c(y1, y[i])
    }
    #if participannt chose arm 2, append to 2-frames
    if (chosen[i]==2){
      X2<-rbind(X2, Xnew)
      y2<-c(y2, y[i])
    }
    #if participannt chose arm 3, append to 3-frames
    if (chosen[i]==3){
      X3<-rbind(X3, Xnew)
      y3<-c(y3, y[i])  
    }
    #if participannt chose arm 4, append to 4-frames
    if (chosen[i]==4){
      X4<-rbind(X4, Xnew)
      y4<-c(y4, y[i])
    }
  }
  utilities<-acquisition(out1, out2, out3, out4)
  utilities <- utilities/max(utilities) #scale to max value of one
  p <- exp(beta*utilities)
  p <- p/rowSums(p)
  #avoid underflow by setting a floor and a ceiling
  p <- (pmax(p, 0.00001))
  p <- (pmin(p, 0.99999))
  #Calculate Negative log likelihood
  nLL <- -sum(log(p[cbind(c(1:nrow(X)),chosen)]))
  return(nLL)  
}


##############################################################################################################
#ACQUISITION FUNCTIONS
##############################################################################################################

#Upper Confidence Bound Sampling
ucb<-function(out1, out2, out3, out4){
  #calulate all the upper confidence bounds
  ucb1<-out1$mu+2*out1$sig
  ucb2<-out2$mu+2*out2$sig
  ucb3<-out3$mu+2*out3$sig
  ucb4<-out4$mu+2*out4$sig
  #bind them together
  outtotal<-cbind(ucb1,ucb2,ucb3,ucb3)
  #return them
  return(outtotal)
}

#Thompson sampling
thompson<-function(out1, out2, out3, out4){
  #initialize
  outtotal<-matrix(0, nrow=nrow(out1), ncol=4)
  mus<-cbind(out1$mu, out2$mu, out3$mu, out4$mu)
  sigs<-cbind(out1$sig, out2$sig, out3$sig, out4$sig)
  #loop through
  for (i in 1:nrow(out1)){
    #simulate from normals given the means and variances
    simout<-apply(cbind(mus[i,], sigs[i,]),1, function(x){rnorm(n=1000,mean=x[1],sd=x[2])})
    #collect mean simulated wins per arm
    outtotal[i,]<-apply(simout==apply(simout,1,max),2,mean)
  }
  #return them
  return(outtotal)
}

#Probability of Improvement
probofimp<-function(out1, out2, out3, out4){
  #We'll compare against the middle value 0.5
  y.star<-0.5
  #intialize
  pi1<-pi2<-pi3<-pi4<-rep(0, nrow(out1))
  #calulate the probabilities of improvement
  pi1<-pnorm((out1$mu-y.star)/out1$sig)
  pi2<-pnorm((out2$mu-y.star)/out2$sig)
  pi3<-pnorm((out3$mu-y.star)/out3$sig)
  pi4<-pnorm((out4$mu-y.star)/out4$sig)
  #bind them
  outtotal<-cbind(pi1, pi2, pi3, pi4)
  #return them
  return(outtotal)
}

#Expected improvement
exofimp<-function(out1, out2, out3, out4){
  ##We'll compare against the middle value 0.5
  y.star<-0.5
  #initialize
  ei1<-ei2<-ei3<-ei4<-rep(0, nrow(out1))
  #calulate z-scores first, then expected improvements
  z1<-(out1$mu-y.star)/out1$sig
  ei1 <-(out1$mu-y.star)*pnorm(z1)+out1$sig*dnorm(z1)
  z2<-(out2$mu-y.star)/out2$sig
  ei2 <-(out2$mu-y.star)*pnorm(z2)+out2$sig*dnorm(z2)
  z3<-(out3$mu-y.star)/out3$sig
  ei3 <-(out3$mu-y.star)*pnorm(z3)+out3$sig*dnorm(z3)
  z4<-(out4$mu-y.star)/out4$sig
  ei4 <-(out4$mu-y.star)*pnorm(z4)+out4$sig*dnorm(z4)
  #adjust if they are bigger than 100
  ei1<-ifelse(ei1>1,1,ei1)
  ei2<-ifelse(ei2>1,1,ei2)
  ei3<-ifelse(ei3>1,1,ei3)
  ei4<-ifelse(ei4>1,1,ei4)
  #bind them
  outtotal<-cbind(ei1, ei2, ei3, ei4)
  #return them
  return(outtotal)
}

##############################################################################################################
#OPTIMIZATION FUNCTION
##############################################################################################################
#function to pluck in to the optimaztion routine
#x: scalar, indicates participant's id
#kernel, function, can be "rbf", "oru", or "matern"
#acquisition, function, can be "ucb", "probofimp", "expofimp", or "thompson"
optfun<-function(x, n, kernel, acquisition){
  #subselect participant
  d1<-subset(data, id==x & trial<n)
  #assign the needed vectors
  chosen <- d1$chosen
  y  <- d1$y
  x1 <- d1$x1
  x2 <- d1$x2
  x3 <- d1$x3
  #create observation matrix
  X<-as.matrix(cbind(x1,x2,x3))
  #bounds
  ubound <- rep(5, 3)
  lbound <- rep(-5, 3)
  #annonymous function
  fn<-function(x){gpbandit(x, X, chosen, y, kernel, acquisition)}
  #optimization
  fit <- DEoptim(fn=fn, lower = lbound, upper=ubound, DEoptim.control(itermax=100))
  #return optimized value
  output <- c(fit$optim$bestval, fit$optim$bestmem)
  return(output)
}



##############################################################################################################
#PREDICTION FUNCTION
##############################################################################################################
banditpredict<-function(x, pars, n, kernel, acquisition){
  #get participant
  d1<-subset(data, id==x)
  par<-exp(pars)
  #gp parameters are the ones up to the penultimate
  pargp<-par[1:(length(par-1))]
  #last parameter is inverse temperature for softmax
  beta<-par[length(par)]
  dlearn<-d1[1:(n-1),]
  dtest<-d1[n,]
  Xnew<-cbind(dtest$x1, dtest$x2, dtest$x3)
  out1<-with(subset(dlearn, chosen==1), {gpr(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y, kernel=kernel)})
  out2<-with(subset(dlearn, chosen==2), {gpr(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y, kernel=kernel)})
  out3<-with(subset(dlearn, chosen==3), {gpr(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y, kernel=kernel)})
  out4<-with(subset(dlearn, chosen==4), {gpr(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y, kernel=kernel)})
  
  outtotal<-acquisition(out1, out2, out3, out4)
  outtotal<-outtotal/max(outtotal)
  p <- exp(beta*outtotal)
  p <- p/rowSums(p)
  mcfad <- 1-1*sum(log(p[,dtest$chosen]))/(log(0.25))
  mcfad<-ifelse(mcfad<0,0,mcfad)
  return(mcfad)
}

##############################################################################################################
#OPTIMIZATION ROUTINE
##############################################################################################################
#create list of all kernel functions
kernellist<-list(oru,rbf)
#names of all kernel functions
kernelnames<-c("Oru","RBF")
#list of all acquisition functions
acqlist<-list(ucb,thompson, probofimp, exofimp)
#names of all acquisition functions
acqnames<-c("UCB", "Thompon", "ProbOfImp", "ExpectedImp")
#all combinations of kernels and acquisition functions will be needed
combs<-expand.grid(1:2, 1:4, seq(10,150,20), 1:max(data$id))
names(combs)<-c("kernel", "acquisition", "n", "id")

taskid <- as.integer(commandArgs(TRUE)[1])
taskid<-ifelse(is.na(taskid),1,taskid)

pars<-optfun(x=combs$id[taskid], n=combs$n[taskid], kernel=kernellist[[combs$kernel[taskid]]], acquisition = acqlist[[combs$acquisition[taskid]]])
perror<-banditpredict(x=combs$id[taskid], pars=pars[2:4], n=combs$n[taskid], kernel=kernellist[[combs$kernel[taskid]]], acquisition = acqlist[[combs$acquisition[taskid]]])

dat<-data.frame(id=combs$id[taskid], n=combs$n[taskid], kernel=kernelnames[combs$kernel[taskid]], acquisition=acqnames[combs$acquisition[taskid]], 
           logloss=pars[1], lambda=exp(pars[2]), sigma=exp(pars[3]),tau=exp(pars[4]), predloss=perror)

rownames(dat)<-NULL

write.csv(dat, paste0("/home/ucabchu/Scratch/cmab/gpdat",taskid,".csv"))
          