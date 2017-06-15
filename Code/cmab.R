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

source("/home/hanshalbe/Desktop/cmabpred/getdata.R")

data1<-getdata(1)
data1$data<-"discrete"
data2<-getdata(2)
data2$data<-"linear"
data3<-getdata(3)
data3$data<-"gp"
dim(data1)
dim(data2)
dim(data3)
data<-rbind(data1, data2, data3)
data$id<-rep(1:(nrow(data)/150), each=150)

source("/home/hanshalbe/Desktop/cmabpred/models.R")
source("/home/hanshalbe/Desktop/cmabpred/acquisitionfunctions.R")



##############################################################################################################
#BANDIT
##############################################################################################################
bandit<-function(par, X, chosen, y, model, acquisition){
  
  #exponentiate parameters as par>0
  par<-exp(par)
  
  #if UCB is used, then the last parameter is beta
  if(class(acquisition)[2]=="UCB"){
    beta<-par[length(par)]
    tau<-par[(length(par)-1)]}
  #if UCB isn't
  if(class(acquisition)[2]!="UCB"){
    tau<-par[length(par)]
    }
  
  if(class(model)[2]=="GP"){
    #gp parameters are the ones up to the penultimate
    pargp<-par[1:2]
    #last parameter is inverse temperature for softmax
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
      out1<-rbind(out1, gp(X.test=Xnew, theta=pargp, X=X1, Y=y1))    
      #prediction for 2nd arm based on new observations
      out2<-rbind(out2, gp(X.test=Xnew, theta=pargp, X=X2, Y=y2))    
      #prediction for 3rd arm based on new observations
      out3<-rbind(out3, gp(X.test=Xnew, theta=pargp, X=X3, Y=y3))    
      #prediction for 4th arm based on new observations
      out4<-rbind(out4, gp(X.test=Xnew, theta=pargp, X=X4, Y=y4))  
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
  }
  
  if(class(model)[2]=="KAL"){
    
    parkal<-par[1:2]
    out1<-out2<-out3<-out4<-data.frame(mu=numeric(), sig=numeric())
    #loop through observations
    for (i in 1:length(y)){
      #new observation
      out<-model(theta=parkal, chosen=chosen[1:i], y=y[1:i], x=X)
      #prediction for 1st arm based on new observations
      out1<-rbind(out1, out[1,])    
      #prediction for 2nd arm based on new observations
      out2<-rbind(out2, out[2,])    
      #prediction for 3rd arm based on new observations
      out3<-rbind(out3, out[3,])    
      #prediction for 4th arm based on new observations
      out4<-rbind(out4, out[4,])  
    }
  }
  
  if(class(model)[2]=="MEB"){
    parmb<-par[1]
    out1<-out2<-out3<-out4<-data.frame(mu=numeric(), sig=numeric())
    #loop through observations
    for (i in 1:length(y)){
      #new observation
      out<-model(theta=parmb, chosen=chosen[1:i], y=y[1:i], x=X)
      #prediction for 1st arm based on new observations
      out1<-rbind(out1, out[1,])    
      #prediction for 2nd arm based on new observations
      out2<-rbind(out2, out[2,])    
      #prediction for 3rd arm based on new observations
      out3<-rbind(out3, out[3,])    
      #prediction for 4th arm based on new observations
      out4<-rbind(out4, out[4,])  
    }
  }
  
  if(class(model)[2]=="LIN"){
    #initialize ys
    y1<-y2<-y3<-y4<-c(0)
    #intialize Xsscp -r user@your.server.example.com:/path/to/foo /home/user/Desktop/
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
      out1<-rbind(out1, lin(X.test=Xnew, X=X1, Y=y1))    
      #prediction for 2nd arm based on new observations
      out2<-rbind(out2, lin(X.test=Xnew, X=X2, Y=y2))    
      #prediction for 3rd arm based on new observations
      out3<-rbind(out3,lin(X.test=Xnew, X=X3, Y=y3))    
      #prediction for 4th arm based on new observations
      out4<-rbind(out4, lin(X.test=Xnew, X=X4, Y=y4))  
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
  }
  
  if(class(acquisition)[2]=="UCB"){utilities<-acquisition(out1, out2, out3, out4, beta)}
  if(class(acquisition)[2]!="UCB"){utilities<-acquisition(out1, out2, out3, out4)}
  
  utilities <- utilities/max(utilities) #scale to max value of one
  p <- exp(tau*utilities)
  p <- p/rowSums(p)
  #avoid underflow by setting a floor and a ceiling
  p <- (pmax(p, 0.00001))
  p <- (pmin(p, 0.99999))
  #Calculate Negative log likelihood
  nLL <- -sum(log(p[cbind(c(1:nrow(X)),chosen)]))
  return(nLL)  
}

##############################################################################################################
#OPTIMIZATION FUNCTION
##############################################################################################################
#function to pluck in to the optimaztion routine
#x: scalar, indicates participant's id
#kernel, function, can be "rbf", "oru", or "matern"
#acquisition, function, can be "ucb", "probofimp", "expofimp", or "thompson"
optfun<-function(x, n, model, acquisition){
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
  if(class(model)[2]=="LIN"){
  ubound <- rep(5, 1)
  lbound <- rep(-5, 1)
  }
  if(class(model)[2]=="GP"){
    ubound <- rep(5, 3)
    lbound <- rep(-5, 3)
  }
  if(class(model)[2]=="KAL"){
    ubound <- rep(5, 3)
    lbound <- rep(-5, 3)
  }
  if(class(model)[2]=="MEB"){
    ubound <- rep(5, 2)
    lbound <- rep(-5, 2)
  }
  
  if(class(acquisition)[2]=="UCB"){
    ubound <- c(ubound,5)
    lbound <- c(lbound,-5)
  }
  #annonymous function
  fn<-function(x){bandit(x, X, chosen, y, model, acquisition)}
  #optimization
  fit <- DEoptim(fn=fn, lower = lbound, upper=ubound, DEoptim.control(itermax=10))
  #return optimized value
  output <- c(fit$optim$bestval, fit$optim$bestmem)
  return(output)
}


##############################################################################################################
#PREDICTION FUNCTION
##############################################################################################################
banditpredict<-function(x, par, n, model, acquisition){
  #get participant
  d1<-subset(data, id==x)
  #last parameter is inverse temperature for softmax
  dlearn<-d1[1:(n-1),]
  dtest<-d1[n,]
  #exponentiate parameters as par>0
  par<-exp(par)
  
  #if UCB is used, then the last parameter is beta
  if(class(acquisition)[2]=="UCB"){
    beta<-par[length(par)]
    tau<-par[(length(par)-1)]}
  #if UCB isn't
  if(class(acquisition)[2]!="UCB"){
    tau<-par[length(par)]
  }
  
  if(class(model)[2]=="GP"){
    #gp parameters are the ones up to the penultimate
    pargp<-par[1:2]
    Xnew<-cbind(dtest$x1, dtest$x2, dtest$x3)
    out1<-with(subset(dlearn, chosen==1), {gp(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y)})
    out1$sig<-ifelse(is.na(out1$sig),pargp[2],out1$sig)
    out2<-with(subset(dlearn, chosen==2), {gp(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y)})
    out2$sig<-ifelse(is.na(out2$sig),pargp[2],out2$sig)
    out3<-with(subset(dlearn, chosen==3), {gp(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y)})
    out3$sig<-ifelse(is.na(out3$sig),pargp[2],out3$sig)
    out4<-with(subset(dlearn, chosen==4), {gp(X.test=Xnew, theta=pargp, X=cbind(x1,x2,x3), Y=y)})
    out4$sig<-ifelse(is.na(out4$sig),pargp[2],out4$sig)
    
  }
  if(class(model)[2]=="KAL"){
   
    parkal<-par[1:2]
    out<-model(theta=parkal, chosen=dlearn$chosen, y=dlearn$y, x=cbind(dlearn$x1,dlearn$x2,dlearn$x3))
    #prediction for 1st arm based on new observations
    out1<-out[1,]    
    #prediction for 2nd arm based on new observations
    out2<-out[2,]
    #prediction for 3rd arm based on new observations
    out3<-out[3,]    
    #prediction for 4th arm based on new observations
    out4<-out[4,]  
  }
  if(class(model)[2]=="MEB"){
    
    parmb<-par[1]
    out<-model(theta=parmb, chosen=dlearn$chosen, y=dlearn$y, x=cbind(dlearn$x1,dlearn$x2,dlearn$x3))
    #prediction for 1st arm based on new observations
    out1<-out[1,]    
    #prediction for 2nd arm based on new observations
    out2<-out[2,]
    #prediction for 3rd arm based on new observations
    out3<-out[3,]    
    #prediction for 4th arm based on new observations
    out4<-out[4,]  
  }
  if(class(model)[2]=="LIN"){
    #gp parameters are the ones up to the penultimate
    Xnew<-cbind(dtest$x1, dtest$x2, dtest$x3)
    out1<-with(subset(dlearn, chosen==1), {lin(X.test=Xnew, X=cbind(x1,x2,x3), Y=y)})
    out2<-with(subset(dlearn, chosen==2), {lin(X.test=Xnew, X=cbind(x1,x2,x3), Y=y)})
    out3<-with(subset(dlearn, chosen==3), {lin(X.test=Xnew, X=cbind(x1,x2,x3), Y=y)})
    out4<-with(subset(dlearn, chosen==4), {lin(X.test=Xnew, X=cbind(x1,x2,x3), Y=y)})
  }
  
  if(class(acquisition)[2]=="UCB"){utilities<-acquisition(out1, out2, out3, out4, beta)}
  if(class(acquisition)[2]!="UCB"){utilities<-acquisition(out1, out2, out3, out4)}
  utilities <- utilities/max(utilities) #scale to max value of one
  p <- exp(tau*utilities)
  p <- p/rowSums(p)
  mcfad <- 1-1*sum(log(p[,dtest$chosen]))/(log(0.25))
  mcfad<-ifelse(mcfad<0,0,mcfad)
  return(mcfad)
}


##############################################################################################################
#OPTIMIZATION ROUTINE
##############################################################################################################
#create list of all kernel functions
modellist<-list(meanbayes,kalman, lin, gp)
#names of all kernel functions
modelnames<-c("MEB","KAL", "LIN", "GP")
#list of all acquisition functions
acqlist<-list(ucb,thompson, probofimp, exofimp)
#names of all acquisition functions
acqnames<-c("UCB", "Thompon", "ProbOfImp", "ExpectedImp")
#all combinations of kernels and acquisition functions will be needed
combs<-expand.grid(1:4, 1:4, seq(10,150,20), 1:max(data$id))
names(combs)<-c("model", "acquisition", "n", "id")

taskid <- as.integer(commandArgs(TRUE)[1])
taskid<-4

pars<-optfun(x=combs$id[taskid], n=combs$n[taskid], model=modellist[[combs$model[taskid]]], acquisition = acqlist[[combs$acquisition[taskid]]])
perror<-banditpredict(x=combs$id[taskid], par=pars[-1], n=combs$n[taskid], model=modellist[[combs$model[taskid]]], acquisition = acqlist[[combs$acquisition[taskid]]])

if(class(acqlist[[combs$acquisition[taskid]]])[2]=="UCB"){
  beta<-as.numeric(exp(pars[length(pars)]))
  tau<-as.numeric(exp(pars[(length(pars)-1)]))}
#if UCB isn't
if(class(acqlist[[combs$acquisition[taskid]]])[2]!="UCB"){
  tau<-as.numeric(exp(pars[(length(pars))]))
  beta<-NA
}

if(class(modellist[[combs$model[taskid]]])[2]=="LIN"){
 par1<-NA
 par2<-NA
}
if(class(modellist[[combs$model[taskid]]])[2]=="KAL"){
  par1<-as.numeric(exp(pars[1]))
  par2<-as.numeric(exp(pars[2]))
}  
if(class(modellist[[combs$model[taskid]]])[2]=="GP"){
  par1<-as.numeric(exp(pars[1]))
  par2<-as.numeric(exp(pars[2]))
}  
if(class(modellist[[combs$model[taskid]]])[2]=="MEB"){
  par1<-as.numeric(exp(pars[1]))
  par2<-NA
}  

dat<-data.frame(id=combs$id[taskid], data=subset(data, id==combs$id[taskid])$data[1], n=combs$n[taskid], model=modelnames[combs$model[taskid]], acquisition=acqnames[combs$acquisition[taskid]], 
           logloss=pars[1], par1=par1, par2=par2, tau=tau, beta=beta, predloss=perror)

rownames(dat)<-NULL
sdat<-read.csv("gpdat.csv")
sdat$X<-NULL
sdat<-read.csv("/home/ucabchu/Scratch/cmab2/gpdat.csv")
sdat<-rbind(sdat, dat)
write.csv(sdat, "/home/ucabchu/Scratch/cmab2/gpdat.csv")
          