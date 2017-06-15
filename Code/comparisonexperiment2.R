#Bayesian Mean
bmei<-read.csv("/home/hanshalbe/Desktop/cmab/data2/BayesMeanExpectedImp.csv")
bmpi<-read.csv("/home/hanshalbe/Desktop/cmab/data2/BayesMeanProbOfImp.csv")
bmth<-read.csv("/home/hanshalbe/Desktop/cmab/data2/BayesMeanThompon.csv")
bmuc<-read.csv("/home/hanshalbe/Desktop/cmab/data2/BayesMeanUCB.csv")

#Kalman Filter
kfei<-read.csv("/home/hanshalbe/Desktop/cmab/data2/KalmanExpectedImp.csv")
kfpi<-read.csv("/home/hanshalbe/Desktop/cmab/data2/KalmanProbOfImp.csv")
kfth<-read.csv("/home/hanshalbe/Desktop/cmab/data2/KalmanThompon.csv")
kfuc<-read.csv("/home/hanshalbe/Desktop/cmab/data2/KalmanUCB.csv")


#Linear Regression
mLL <- function(par,chosen,V) {
  par<-exp(par)
  luce <- exp(par*V)
  luce[is.na(luce)]<-0.0001
  luce <- luce/rowSums(luce)
  luce<-(pmax(luce, 0.00001))
  luce<-(pmin(luce, 0.99999))
  liglik<- -2*sum(log(luce[cbind(1:nrow(luce),chosen)]))
  return(liglik)
}


liei<-rep(415, 59)
lipi<-rep(415, 59)
lith<-rep(415, 59)
liuc<-rep(415, 59)

for (folks in 1:59){
  
  
  
  dat<-read.csv(paste0("/home/hanshalbe/Desktop/cmab/data2/lin/p",folks, "lin.csv"))
  
  #ucb
  ucb1<-dat[,1]+2*sqrt(dat[,5])
  ucb1<-ifelse(ucb1>100,100,ucb1)
  ucb1<-ifelse(ucb1<0,0,ucb1)
  ucb2<-dat[,2]+2*sqrt(dat[,6])
  ucb2<-ifelse(ucb2>100,100,ucb2)
  ucb2<-ifelse(ucb2<0,0,ucb2)
  ucb3<-dat[,3]+2*sqrt(dat[,7])
  ucb3<-ifelse(ucb3>100,100,ucb3)
  ucb3<-ifelse(ucb3<0,0,ucb3)
  ucb4<-dat[,4]+2*sqrt(dat[,8])
  ucb4<-ifelse(ucb4>100,100,ucb4)
  ucb4<-ifelse(ucb4<0,0,ucb4)
  mucb<-cbind(ucb1, ucb2, ucb3, ucb4)
  liuc[folks]<- optim(0, mLL, method="BFGS", V=mucb, chosen=dat[,12], hessian=TRUE)$value 
  
  
  #probability of improvement
  y.star<-0
  pi1<-pi2<-pi3<-pi4<-rep(0, nrow(dat))
  for (i in 1:nrow(dat)){
    pi1[i]<-pnorm((dat[i,1]-y.star)/sqrt(dat[i,5]))
    pi2[i]<-pnorm((dat[i,2]-y.star)/sqrt(dat[i,6]))
    pi3[i]<-pnorm((dat[i,3]-y.star)/sqrt(dat[i,7]))
    pi4[i]<-pnorm((dat[i,4]-y.star)/sqrt(dat[i,8]))
    y.star<-ifelse(dat[i,13]>y.star, dat[i,13], y.star)
  }
  mpi<-cbind(pi1, pi2, pi3, pi4)
  lipi[folks]<-optim(1, mLL, method="BFGS", V=mpi, chosen=dat[,12], hessian=TRUE)$value 
  
  #expected improvement
  y.star<-0
  ei1<-ei2<-ei3<-ei4<-rep(0, nrow(dat))
  for (i in 1:nrow(dat)){
    z1<-(dat[i,1]-y.star)/sqrt(dat[i,5])
    ei1[i] <-(dat[i,1]-y.star)*pnorm(z1)+sqrt(dat[i,5])*dnorm(z1)
    z2<-(dat[i,2]-y.star)/sqrt(dat[i,6])
    ei2[i] <-(dat[i,2]-y.star)*pnorm(z2)+sqrt(dat[i,6])*dnorm(z2)
    z3<-(dat[i,3]-y.star)/sqrt(dat[i,7])
    ei3[i] <-(dat[i,3]-y.star)*pnorm(z3)+sqrt(dat[i,7])*dnorm(z3)
    z4<-(dat[i,4]-y.star)/sqrt(dat[i,8])
    ei4[i] <-(dat[i,4]-y.star)*pnorm(z4)+sqrt(dat[i,8])*dnorm(z4)
    y.star<-ifelse(dat[i,13]>y.star, dat[i,13], y.star)
  }
  ei1<-ifelse(ei1>100,100,ei1)
  ei2<-ifelse(ei2>100,100,ei2)
  ei3<-ifelse(ei3>100,100,ei3)
  ei4<-ifelse(ei4>100,100,ei4)
  mei<-cbind(ei1, ei2, ei3, ei4)
  liei[folks]<- optim(0, mLL, method="BFGS", V=mei, chosen=dat[,12], hessian=TRUE)$value 
  
  #Thompson
  mth<-mei
  for (i in 1:nrow(dat)){
    simout<-apply(cbind(t(dat[i,1:4]), t(dat[i,5:8])),1, function(x){rnorm(n=1000,mean=x[1],sd=x[2])})
    mth[i,]<-apply(simout==apply(simout,1,max),2,mean)
  }
  lith[folks]<- optim(1, mLL, method="BFGS", V=mth, chosen=dat[,12])$value 
  
}

#Ornstein Uhlenbeck
ouei<-read.csv("/home/hanshalbe/Desktop/cmab/data2/OruExpectedImp.csv")
oupi<-read.csv("/home/hanshalbe/Desktop/cmab/data2/OruProbOfImp.csv")
outh<-read.csv("/home/hanshalbe/Desktop/cmab/data2/OruThompon.csv")
ouuc<-read.csv("/home/hanshalbe/Desktop/cmab/data2/OruUCB.csv")

#Matern
maei<-read.csv("/home/hanshalbe/Desktop/cmab/data2/MaternExpectedImp.csv")
mapi<-read.csv("/home/hanshalbe/Desktop/cmab/data2/MaternProbOfImp.csv")
math<-read.csv("/home/hanshalbe/Desktop/cmab/data2/MaternThompon.csv")
mauc<-read.csv("/home/hanshalbe/Desktop/cmab/data2/MaternUCB.csv")

#RBF
rbei<-read.csv("/home/hanshalbe/Desktop/cmab/data2/RBFExpectedImp.csv")
rbpi<-read.csv("/home/hanshalbe/Desktop/cmab/data2/RBFProbOfImp.csv")
rbth<-read.csv("/home/hanshalbe/Desktop/cmab/data2/RBFThompon.csv")
rbuc<-read.csv("/home/hanshalbe/Desktop/cmab/data2/RBFUCB.csv")

dfinal<-data.frame(aic=c(bmei$x, bmpi$x, bmth$x, bmuc$x, kfei$x, kfpi$x, kfth$x, kfuc$x, liei, lipi, lith, liuc, ouei$x, oupi$x, outh$x, ouuc$x, maei$x, mapi$x, math$x, mauc$x, rbei$x, rbpi$x, rbth$x, rbuc$x),
                   model=rep(c("Bayesian Mean", "Kalman Filter", "Linear Regression", "Ornstein-Uhlenbeck", "Matérn 3/2", "Radial Basis"), each=59*4),
                   acquisition=rep(rep(c("Expected Improvement", "Probability of Improvement", "Probability of max. Utility", "Upper Confidence Bound"), each=59), 6),
                   id=rep(1:59, 24))

head(dfinal)
length(indicator)

dfinal<-subset(dfinal, id%in%indicator)

dfinal$evidence<- -0.5*dfinal$aic
m<-subset(dfinal, id==1)$evidence
for (i in 2:length(indicator)){
  m<-rbind(m, subset(dfinal, id==i)$evidence) 
}
dim(m)
m<-cbind(m, rep(150*log(0.25),55))
write.csv(m, "/home/hanshalbe/Desktop/cmab/data2/evidence.csv")

library(plyr)
dplot<-ddply(dfinal,~acquisition+model,summarise,mu=mean(aic))

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

dplot$mu1<-1-range01(dplot$mu)/0.98
dplot$mu<-round(dplot$mu)
dplot$Measure<-"AIC"

protected<-read.csv("/home/hanshalbe/Desktop/cmab/data2/protected2.csv", header=FALSE)
dummy<-data.frame(p=as.numeric(protected[1:24]), model=subset(dfinal, id==1)$model, acq=subset(dfinal, id==1)$acquisition)
dummy2<-ddply(dummy,~acq+model,summarise,mu=mean(p))

dplot2<-dplot
dplot2$Measure<-"Exceedance Probability"
dplot2$mu<-dummy2$mu
dplot2$mu1<-range01(dplot2$mu)*0.98+0.01


dfinal$best<-ave(dfinal$aic, dfinal$id, FUN=min)==dfinal$aic
dplot3<-dplot2
dplot3$mu<-ddply(dfinal,~acquisition+model,summarise,mu=sum(best))$mu
dplot3$mu1<-dplot3$mu/sum(dplot3$mu)
dplot3$mu1<-range01(dplot3$mu1)*0.98+0.01
dplot3$Measure<-"#Best"

dp<-rbind(dplot, dplot2, dplot3)  
dp$mu<-round(dp$mu*100)/100
library(ggplot2)
dodge <- position_dodge(width=0.5)
dp$model<-factor(dp$model, 
                 levels=c("Bayesian Mean", "Kalman Filter", "Linear Regression", "Ornstein-Uhlenbeck", "Matérn 3/2", "Radial Basis"))

p <- ggplot(dp, aes(y=mu1, x=acquisition, fill=Measure)) + 
  #bars
  geom_bar(position="dodge", stat="identity")+
  #golden ratio error bars
  ggtitle("Model Comparison")+theme_classic() +xlab("\nAcquisition Function")+ylab("Quality\n")+
  scale_fill_manual(values=c("grey25","grey50", "grey75"))+
  scale_y_continuous(breaks = round(seq(0, 1, 0.5), 1), limits=c(-0.1,1.2))+
  geom_text(aes(label=mu), position=position_dodge(width=0.9), vjust=-0.25)+
  #fill with Wes Anderson colors
  theme(legend.position="top",
        text = element_text(size=18, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))+
  facet_wrap(~ model, nrow=8)

pdf("/home/hanshalbe/Desktop/cmab/paper/modelcomparison2.pdf", width=10, height=10)
p
dev.off()
indicator<-tapply(data$y, data$id, function(x){t.test(x, mu=50, conf.level = 0.90)$conf.int[[1]]})
sum(indicator<48)
length(indicator)
dim(m)
indicator<-indicator>48
mnew<-m[indicator,]
dim(mnew)
write.csv(mnew, "/home/hanshalbe/Desktop/cmab/data2/evidencebest2.csv")

protectednew<-as.numeric(read.csv("/home/hanshalbe/Desktop/cmab/data2/protectedbest2.csv", header=FALSE))
dummy<-data.frame(p=as.numeric(protectednew[1:24]), model=subset(dfinal, id==1)$model, acq=subset(dfinal, id==1)$acquisition)
dummy2<-ddply(dummy,~acq+model,summarise,mu=mean(p))
