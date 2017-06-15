d<-read.csv(file.choose(), header=TRUE)
library(plyr)
names(d)

dgp<-subset(d,model =="GP")
dgp$par1<-sqrt(log(dgp$par1))
se<-function(x){return(sd(x)/sqrt(length(x)))}
dp1<-ddply(subset(dgp, acquisition=="UCB"), ~data+n, summarize, mean=mean(predloss), se=se(predloss))
dp2<-ddply(subset(dgp, acquisition=="UCB"), ~data+n, summarize, mean=mean(tau), se=se(tau)*0.5)
dp3<-ddply(subset(dgp, acquisition=="UCB"), ~data+n, summarize, mean=mean(sqrt(beta)), se=se(sqrt(beta))*0.5)
dp4<-ddply(dgp, ~data+n, summarize, mean=mean(sqrt(par2)), se=se(sqrt(par2)))
dp<-rbind(dp1,dp2,dp3,dp4)
dp$parameter<-rep(c("R^2", "tau", "beta", "sigma"), each=nrow(dp1))
dp$parameter<-factor(dp$parameter,levels=c("R^2", "tau", "beta", "sigma"), labels=c("R^2", "tau^-1", "beta", "sigma"))

library(ggplot2)
pd <- position_dodge(.1)
dp$Experiment<-dp$data
dp$Experiment<-mapvalues(dp$Experiment, c("discrete", "gp", "linear"), c("Discrete", "Nonlinear","Linear"))
dp$Experiment<-factor(dp$Experiment, levels=c("Discrete", "Linear", "Nonlinear"))
p<-ggplot(dp, aes(x=n, y=mean, colour=Experiment)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=5, size=1, position=pd) +
  geom_line(position=pd, size=1) +
  geom_point(size=3)+
  ylab("Estimate (± standard error)")+ggtitle("Parameter estimates across trials")+xlab("Time")+
  scale_color_manual(values=c("grey35", "grey60", "black"))+
  theme_classic()+theme(text = element_text(size=24,  family="serif"), legend.position="top")+
  facet_wrap(~parameter, nrow=2, scales="free", labeller=label_parsed)
p

pdf("paracomptrials.pdf", width=10, height=10)
p
dev.off()

db1<-ddply(dgp, ~data, summarize, m=mean(predloss), se=se(predloss))
db2<-ddply(subset(dgp, acquisition=="UCB"), ~data, summarize, m=mean(sqrt(beta)), se=se(sqrt(beta))*0.5)
db3<-ddply(subset(dgp, acquisition=="UCB"), ~data, summarize, m=mean(sqrt(par2)), se=se(sqrt(par2)))
db4<-ddply(subset(dgp, acquisition=="UCB") , ~data, summarize, m=mean(par1), se=se(par1))
dplot<-rbind(db1, db2, db3, db4)
dplot$parameter<-rep(c("R^2", "beta", "sigma", "lambda"), each=nrow(db1))
dplot$parameter<-factor(dplot$parameter,levels=c("R^2", "beta", "sigma", "lambda"), labels=c("R^2", "beta", "sigma", "lambda"))
dplot$Experiment<-dplot$data
dplot$Experiment<-mapvalues(dplot$Experiment, c("discrete", "gp", "linear"), c("Discrete", "Nonlinear","Linear"))
dplot$Experiment<-factor(dplot$Experiment, levels=c("Discrete", "Linear", "Nonlinear"))

pd <- position_dodge(.1)

p2<-ggplot(dplot, aes(x=Experiment, y=m, fill=Experiment)) +
  geom_bar(position="dodge",stat="identity")+
  geom_errorbar(aes(ymin=m-se, ymax=m+se), width=0.4, size=1, position=pd) +
  ylab("Mean stimate (± standard error)")+ggtitle("Parameter comparison")+xlab("Experiment")+
  scale_fill_manual(values=c("grey35", "grey60", "grey22"))+
  theme_classic()+theme(text = element_text(size=24,  family="serif"), legend.position="top")+
  facet_wrap(~parameter, nrow=2, scales="free", labeller=label_parsed)+
  geom_point(size=3)

pdf("paracomp.pdf", width=10, height=10)
p2
dev.off()

dgp<-subset(d,model =="GP" & acquisition=="UCB")
library(BayesFactor)
d<-ddply(dgp, ~id+data, summarize, m=mean(predloss))
t.test(subset(d, data=="discrete")$m, subset(d, data!="discrete")$m, var.equal=TRUE)
bf<-ttestBF(subset(d, data=="discrete")$m, subset(d, data!="discrete")$m, posterior = TRUE, iterations = 1000)
bf[1]^{-1}

t.test(subset(d, data=="linear")$m, subset(d, data=="gp")$m, var.equal=TRUE)
bf<-ttestBF(subset(d, data=="linear")$m, subset(d, data=="gp")$m, posterior = TRUE, iterations = 1000)
bf[1]

dgp

head(dgp)

d<-ddply(dgp, ~data+id, summarize, m=mean(sqrt(beta)))
t.test(subset(d, data=="discrete")$m, subset(d, data=="gp")$m, var.equal=TRUE)
bf<-ttestBF(subset(d, data=="discrete")$m, subset(d, data=="gp")$m, posterior = TRUE, iterations = 1000)
bf[1]
dgp$num<-mapvalues(dgp$data, from=c("discrete","gp","linear"), to=c(1,3,2))
cor.test(as.numeric(dgp$num), dgp$predloss)

d<-ddply(dgp, ~data+id, summarize, m=mean(sqrt(par2)))
t.test(subset(d, data=="discrete")$m, subset(d, data!="linear")$m, var.equal=TRUE)
bf<-ttestBF(subset(d, data=="discrete")$m, subset(d, data!="discrete")$m, posterior = TRUE, iterations = 1000)
bf[1]

t.test(subset(d, data=="linear")$m, subset(d, data!="gp")$m, var.equal=TRUE)
bf<-ttestBF(subset(d, data=="linear")$m, subset(d, data!="gp")$m, posterior = TRUE, iterations = 1000)
bf[1]

d<-ddply(dgp, ~data+id, summarize, m=median(sqrt(log10(par2))))
d$m[is.na(d$m)]<-mean(d$m, na.rm=TRUE)
t.test(subset(d, data=="linear")$m, subset(d, data!="gp")$m, var.equal=TRUE)
bf<-ttestBF(subset(d, data=="linear")$m, subset(d, data!="gp")$m, posterior = TRUE, iterations = 1000)
bf[1]


d<-ddply(dgp, ~id, summarize, m=cor(n, predloss))
d$m[is.na(d$m)]<-mean(d$m, na.rm=TRUE)
t.test(subset(d, data!="discrete")$m, mu=0)
bf<-ttestBF(subset(d, data!="discrete")$m, mu=0, posterior = TRUE, iterations = 1000)
bf[1]
names(dgp)

d<-ddply(dgp, ~id+data, summarize, m=cor(n,tau))
d$m[is.na(d$m)]<-mean(d$m, na.rm=TRUE)
mean(d$m)
t.test(d$m, mu=0)
bf<-ttestBF(d$m, mu=0, posterior = TRUE, iterations = 1000)
bf[1]^(-1)


d<-ddply(dgp, ~id+data, summarize, m=cor(n, sqrt(beta)))
d$m[is.na(d$m)]<-mean(d$m, na.rm=TRUE)
t.test(d$m, mu=0)
bf<-ttestBF(d$m, mu=0, posterior = TRUE, iterations = 1000)
bf[1]^{-1}

d<-ddply(dgp, ~id+data, summarize, m=cor(n, sqrt(par2)))
d$m[is.na(d$m)]<-mean(d$m, na.rm=TRUE)
t.test(subset(d, data=="gp")$m, mu=0)
bf<-ttestBF(subset(d, data=="gp")$m, mu=0, posterior = TRUE, iterations = 1000)
bf[1]^(-1)
