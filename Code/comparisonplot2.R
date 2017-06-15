rm(list=ls())

packages <- c('plyr', 'ggplot2', 'MASS')
lapply(packages, library, character.only = TRUE)
par<-0
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




dcollect<-data.frame(id=numeric, model=numeric(), acq=numeric(), aic=numeric())                  
models<-c("lin", "qua", "cub", "mat1", "mat3", "mat5", "sqe")


for (folks in 1:59){
  
  
  for (m in seq_along(models)){
    
    dat<-read.csv(paste0("/home/hanshalbe/Dropbox/CMAB/data/raw2/p",folks,models[m],".csv"), header=FALSE)
    
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
    ucbaic<- optim(0, mLL, method="BFGS", V=mucb, chosen=dat[,12], hessian=TRUE)$value 
    
    
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
    piaic<-optim(1, mLL, method="BFGS", V=mpi, chosen=dat[,12], hessian=TRUE)$value 
    
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
    eiaic<- optim(0, mLL, method="BFGS", V=mei, chosen=dat[,12], hessian=TRUE)$value 
    
    #Thompson
    mth<-mei
    for (i in 1:nrow(dat)){
      simout<-apply(cbind(t(dat[i,1:4]), t(dat[i,5:8])),1, function(x){rnorm(n=1000,mean=x[1],sd=x[2])})
      mth[i,]<-apply(simout==apply(simout,1,max),2,mean)
    }
    thaic<- optim(1, mLL, method="BFGS", V=mth, chosen=dat[,12], hessian=TRUE)$value 
    
    #random
    mrand<-matrix(1, nrow=nrow(mucb), ncol=ncol(mucb))
    ranaic<-optim(1, mLL, method="BFGS", V=mrand, chosen=dat[,12], hessian=TRUE)$value
    
    dummy<-data.frame(id=rep(folks, 5),
                      model=rep(models[m], 5),
                      acq=c("random", "ucb", "poi", "ei", "thomp"),
                      aic=c(ranaic , ucbaic, piaic, eiaic, thaic))
    dcollect<-rbind(dcollect, dummy)
  }
}

dcollect[is.na(dcollect)]<-ranaic
head(dcollect,30)
d<-na.omit(dcollect)
d<-subset(dcollect, acq!="random")
d[d$aic==ave(d$aic, d$id, FUN=min),]
d$evidence<- -0.5*d$aic

m<-unique(paste0(d$acq, d$model))
m<-subset(d, id==1)$evidence
for (i in 2:max(d$id)){
  m<-rbind(m, subset(d, id==i)$evidence) 
}

apply(m,2,mean)

write.csv(m, "/home/hanshalbe/evidence2.csv")
names(d)
se<-function(x){return(sd(x)/sqrt(length(x)))}
dplot<-ddply(d,~acq+model,summarise,mu=mean(aic))





dodge <- position_dodge(width=0.5)
#do the plot
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
dplot$mu<-1-range01(dplot$mu)/0.9
dplot$acq<-mapvalues(dplot$acq, unique(dplot$acq), c("ExI", "PoI", "Thompson", "UCB"))
dplot$model<-mapvalues(dplot$model, unique(dplot$model), 
                       c("Linear", "Quadratic", "Cubic", "Matérn 1/2", "Matérn 3/2", "Matérn 5/2", "RBF"))

protected<-read.csv("/home/hanshalbe/Dropbox/CMAB/data/protected2.csv", header=FALSE)
dummy<-data.frame(p=as.numeric(protected), model=subset(d, id==1)$model, acq=subset(d, id==1)$acq)
dummy2<-ddply(dummy,~acq+model,summarise,mu=mean(p))

dplot2<-dplot
dplot2$Measure<-"Exceedance Probability"
dplot2$mu<-dummy2$mu
dplot2$mu<-range01(dplot2$mu)*0.9+0.1


dplot$Measure<-"AIC"


d$best<-ave(d$aic, d$id, FUN=min)==d$aic
dplot3<-dplot2
dplot3$mu<-ddply(d,~acq+model,summarise,mu=sum(best))$mu
dplot3$mu<-dplot3$mu/sum(dplot3$mu)
dplot3$mu<-range01(dplot3$mu)*0.9+0.1
dplot3$Measure<-"#Best"

dp<-rbind(dplot, dplot2, dplot3)  

p <- ggplot(dp, aes(y=mu, x=acq, fill=Measure)) + 
  #bars
  geom_bar(position="dodge", stat="identity")+
  #golden ratio error bars
  ggtitle("Model Comparison")+theme_classic() +xlab("\nAcquisition Function")+ylab("Quality\n")+
  scale_fill_manual(values=c("grey25","grey50", "grey75"))+
  #fill with Wes Anderson colors
  theme(legend.position="top",
        text = element_text(size=18, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))+
  facet_wrap(~ model, nrow=8)

pdf("modelcomparison2.pdf", width=10, height=10)
p
dev.off()
dp
