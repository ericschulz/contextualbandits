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
packages <- c('RCurl', 'devtools', 'matrixcalc', 'DEoptim')
#load them
lapply(packages, library, character.only = TRUE)

#Data Munging:
#read in raw data
mysource1 <- getURL("https://raw.githubusercontent.com/ericschulz/contextualbandits/master/Data/exp1raw.csv")
dat1<- read.csv(text = mysource1, header=TRUE)

#only keep people who have completed the task
dat1 <- dat1[!is.na(dat1$Answer.total),]

#create a data frame with participants' data
data<-data.frame(id= numeric(), x1=numeric(), x2=numeric(), x3=numeric(), chosen=numeric(), y=numeric())
#get into right format
for (i in 1:nrow(dat1)){
  d1<-dat1[i,]
  #chosen arm
  chosen <- as.numeric(strsplit(as.character(d1$Answer.chosendeck), ",")[[1]])
  #outcome
  y <- as.numeric(strsplit(as.character(d1$Answer.reward), ",")[[1]])
  #features
  x1 <- as.numeric(strsplit(as.character(d1$Answer.sunseen), ",")[[1]])
  x2 <- as.numeric(strsplit(as.character(d1$Answer.tempseen), ",")[[1]])
  x3 <- as.numeric(strsplit(as.character(d1$Answer.rainseen), ",")[[1]])
  #dummy frame
  dummy<-data.frame(id=rep(i, length(x1)), x1=x1, x2=x2, x3=x3, chosen=chosen, y=y)
  #concatenate
  data<-rbind(data, dummy)
}

#values bigger than 100 were presented as 100
data$y<-ifelse(data$y>100, 100, data$y)
#values smaller than 0 were presented as 0
data$y<-ifelse(data$y<0, 0, data$y)
#standardize to be between -1 and 1
data$y<-(data$y-50)/100
#trial number
data$trial<-rep(1:150, nrow(data)/150)


#Kalman Filter
#X.test: matrix for predcitions
#theta: vector of hyper-parameters
#X; matrix of observations
#y: vector of observed outcomes
KalmanFilter<- function(theta, chosen, Y){
  library('FKF')
  #set of choices to make predictions about
  # parameters
  mu0 <- 0 #par[1]
  var0 <- 5 #exp(par[1])
  vart <- theta[1] #innovation variance
  vare <- theta[2] #error varriance
  # setup for fkf
  a0 <- rep(mu0,4) #initial estimates of state variable
  P0 <- var0*diag(4) #covariance matrix of a0
  dt <- matrix(rep(0,4),ncol=1) #intercept of the transition equation
  ct <- matrix(rep(0,4),ncol=1) #intercept of the measurement equation
  Tt <- array(diag(4),dim=c(4,4,1)) # Factor of the transition equation
  Zt <- array(diag(4),dim=c(4,4,1)) # factor of  the measurement equation
  HHt <- array(vart*diag(4),dim=c(4,4,1)) #variance of the innovation and transition equatinos
  GGt <- array(vare*diag(4),dim=c(4,4,1)) #disturbances of the measurement equation
  yt <- matrix(NA,ncol=4,nrow=length(Y)) #observations
  #Loop over observations and add payoffs to yt
  for(i in 1:length(Y)) {
    #for each round, fill in the yt that was chosen with the payoff in Y
    yt[i,chosen[i]] <- Y[i]
  } 
  filt <- fkf(a0,P0,dt,ct,Tt,Zt,HHt,GGt,t(yt))
  E <- t(filt$at) #expectation
  S <- sapply(1:nrow(E), FUN=function(x) sqrt(diag(filt$Pt[,,x])))
  result <- rbind(E[length(Y) +1,], t(S)[length(Y)+1,]) #return expectation and variance for each of the 121 options
  #as a data frame with names mu and sig
  prediction <- as.data.frame(t(result))
  colnames(prediction) <- c("mu", "sig")
  return(prediction)
}



kalmanbandit<-function(par, chosen, y, acquisition){
  #exponentiate parameters as par>0
  par<-exp(par)
  #gp parameters are the ones up to the penultimate
  parkal<-par[1:(length(par-1))]
  #last parameter is inverse temperature for softmax
  beta<-par[length(par)]
  #intiailize out frames
  out1<-out2<-out3<-out4<-data.frame(mu=numeric(), sig=numeric())
  #loop through observations
  for (i in 1:length(y)){
    #new observation
    out<-KalmanFilter(parkal, chosen[1:i], y[1:i])
    #prediction for 1st arm based on new observations
    out1<-rbind(out1, out[1,])    
    #prediction for 2nd arm based on new observations
    out2<-rbind(out2, out[2,])    
    #prediction for 3rd arm based on new observations
    out3<-rbind(out3, out[3,])    
    #prediction for 4th arm based on new observations
    out4<-rbind(out4, out[4,])  
    }
  utilities<-acquisition(out1, out2, out3, out4)
  utilities <- utilities/max(utilities) #scale to max value of one
  p <- exp(beta*utilities)
  p <- p/rowSums(p)
  #avoid underflow by setting a floor and a ceiling
  p <- (pmax(p, 0.00001))
  p <- (pmin(p, 0.99999))
  #Calculate Negative log likelihood
  nLL <- -sum(log(p[cbind(c(1:length(y)),chosen)]))
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
optfun<-function(x, n, acquisition){
  #subselect participant
  d1<-subset(data, id==x & trial<n)
  #assign the needed vectors
  chosen <- d1$chosen
  y  <- d1$y
  #bounds
  ubound <- rep(5, 3)
  lbound <- rep(-5, 3)
  #annonymous function
  fn<-function(x){kalmanbandit(x, chosen, y, acquisition)}
  #optimization
  fit <- DEoptim(fn=fn, lower = lbound, upper=ubound, DEoptim.control(itermax=100))
  #return optimized value
  output <- c(fit$optim$bestval, fit$optim$bestmem)
  return(output)
}



##############################################################################################################
#PREDICTION FUNCTION
##############################################################################################################
banditpredict<-function(x, pars, n, acquisition){
  #get participant
  d1<-subset(data, id==x)
  par<-exp(pars)
  #gp parameters are the ones up to the penultimate
  parkal<-par[1:(length(par-1))]
  #last parameter is inverse temperature for softmax
  beta<-par[length(par)]
  dlearn<-d1[1:(n-1),]
  dtest<-d1[n,]
  out<-KalmanFilter(parkal, dlearn$chosen, dlearn$y)
  out1<-out[1,]
  out2<-out[2,]
  out3<-out[3,]
  out4<-out[4,]
  
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
acqlist<-list(ucb,thompson, probofimp, exofimp)
#names of all acquisition functions
acqnames<-c("UCB", "Thompon", "ProbOfImp", "ExpectedImp")
#all combinations of kernels and acquisition functions will be needed
combs<-expand.grid(1:4, seq(10,150,20), 1:max(data$id))
names(combs)<-c("acquisition", "n", "id")

taskid <- as.integer(commandArgs(TRUE)[1])
taskid<-ifelse(is.na(taskid),1,taskid)

parallelizeme<-function(taskid){
  pars<-optfun(x=combs$id[taskid], n=combs$n[taskid], acquisition = acqlist[[combs$acquisition[taskid]]])
  perror<-banditpredict(x=combs$id[taskid], pars=pars[2:4], n=combs$n[taskid], acquisition = acqlist[[combs$acquisition[taskid]]])
  dat<-data.frame(id=combs$id[taskid], n=combs$n[taskid], acquisition=acqnames[combs$acquisition[taskid]], 
                logloss=pars[1], varino=exp(pars[2]), varerro=exp(pars[3]), tau=exp(pars[4]), predloss=perror)
  rownames(dat)<-NULL
  write.csv(dat, paste0("/home/hanshalbe/Desktop/cmabpred/data/exp1kalman",taskid,".csv"))
  cat(taskid, "is done")
}

mclapply(as.list(1:nrow(combs)), parallelizeme, mc.preschedule=FALSE, mc.cores=3)
