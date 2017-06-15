#Upper Confidence Bound Sampling
ucb<-function(out1, out2, out3, out4, beta){
  #calulate all the upper confidence bounds
  ucb1<-out1$mu+beta*out1$sig
  ucb2<-out2$mu+beta*out2$sig
  ucb3<-out3$mu+beta*out3$sig
  ucb4<-out4$mu+beta*out4$sig
  #bind them together
  outtotal<-cbind(ucb1,ucb2,ucb3,ucb3)
  #return them
  return(outtotal)
}
class(ucb)<- c(class(ucb), "UCB")

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
class(thompson)<- c(class(thompson), "THO")

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
class(probofimp)<- c(class(probofimp), "POI")

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
class(exofimp)<- c(class(exofimp), "EXI")