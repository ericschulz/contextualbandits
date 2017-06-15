#House keeping
rm(list=ls())

#required packages packages
packages <- c('RCurl', 'devtools', 'ggplot2', 'parallel', 'plyr', 'jsonlite')
#load them
lapply(packages, library, character.only = TRUE)

#Data Munging:
#download data
myjson<-fromJSON("https://melt.firebaseio.com/.json")

k<-1
data<-data.frame(id= numeric(), x1=numeric(), x2=numeric(), x3=numeric(), 
                 k1=numeric(), k2=numeric(), k3=numeric(), chosen=numeric(), y=numeric())
for (i in 62:length(myjson)){
  d<-myjson[[i]]
  x1<-d$sunseen 
  x2<-d$tempseen 
  x3<-d$rainseen
  k1<-rep(d$k1, length(x1))
  k2<-rep(d$k2, length(x1))
  k3<-rep(d$k3, length(x1))
  y<-d$reward
  chosen<-d$chosendeck
  #dummy frame
  dummy<-data.frame(id=rep(k, length(x1)), x1=x1, x2=x2, x3=x3, k1=k1, k2=k2, k3=k3, chosen=chosen, y=y)
  #concatenate
  data<-rbind(data, dummy)
  #k increment
  k<-k+1
}

names(data)

tapply(data$x1, data$id, length)
data<-subset(data, id!=3)

indicator<-tapply(data$y, data$id, function(x){t.test(x, mu=50)$conf.int[[1]]})
sum(indicator<45)
indicator<-c(1:max(data$id))[indicator>45]
data<-subset(data, id %in% indicator)

dplot1<-ddply(data,~id,summarise, mu=mean(y))

p<-ggplot(dplot1, aes(mu)) +
  geom_histogram(bins=8, col="grey50")+
  ggtitle("Average Scores")+theme_classic() +xlab("\nMean Score")+ylab("Counts\n")+
  #fill with Wes Anderson colors
  theme(text = element_text(size=18, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))


p

pdf("/home/hanshalbe/Desktop/cmab/paper/figures/score3.pdf")
p
dev.off()


myfuncs<-fromJSON("https://functions.firebaseio.com/.json")

p1<-p2<-p3<-p4<-rep(NA, nrow(data))
for (i in 1:nrow(data)){
  p1[i]<-myfuncs[[3]][[data$k1[i]+1]][myfuncs[[1]][[1]]==data$x1[i] & myfuncs[[2]][[1]]==data$x2[i]]
  p2[i]<-myfuncs[[3]][[data$k2[i]+1]][myfuncs[[1]][[1]]==data$x1[i] & myfuncs[[2]][[1]]==data$x3[i]]
  p3[i]<-myfuncs[[3]][[data$k3[i]+1]][myfuncs[[1]][[1]]==data$x2[i] & myfuncs[[2]][[1]]==data$x3[i]]
  p4[i]<-50
}

comb<-cbind(p1,p2,p3,p4)
pchoic<-apply(comb,1,which.max)



data$best<-ifelse(pchoic==data$chosen,1,0)
data$trials<-rep(1:150, max(data$id)-1)
data$trials<-round(data$trials/10)*10
dplot2<-ddply(data, ~trials, summarize, prop=mean(y))
cor.test(dplot2$trials, dplot2$prop)

p<-ggplot(dplot, aes(x=trials, y=prop)) + 
  geom_point(col="grey",size=1)+
  stat_smooth(method="lm", col="black")+
  geom_errorbar(data=dplot2, aes(x=tens, y=mu, ymin=mu-se, ymax=mu+se), width=1.2, size=1, position=pd, col="red") +
  xlab("Trial")+ylab("Mean")+
  ggtitle("Outcomes over time")+
  theme_bw()+theme(text = element_text(size=22))+ ylim(c(40, 70))+
  annotate("text", x = 40, y = 70, label = "r=0.19, p<0.05", size=8)+
  #fill with Wes Anderson colors
  theme(text = element_text(size=18, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))
p

pdf("meanovertime3.pdf")
print(p)
dev.off()


data$trials<-rep(1:150, nrow(data)/150)
isfour<-function(x) {mean(x==4)}
#data$trials<-round(data$trials/10)*10
dplot2<-ddply(data, ~trials, summarize, prop=isfour(chosen))
dplot2$prop[1]<-0.25
cor.test(dplot2$trials, dplot2$prop)

p<-ggplot(dplot2, aes(x=trials, y=prop)) + 
  geom_point(col="grey",size=3)+
  stat_smooth(method="lm", col="black")+
  xlab("Trial")+ylab("Mean")+
  ggtitle("Proportion of 4th arm")+
  theme_classic()+theme(text = element_text(size=22))+ ylim(c(0,0.4))+
  annotate("text", x = 40, y = 0.35, label = "r=0.-0.22, p<0.05", size=8)+
  #fill with Wes Anderson colors
  theme(text = element_text(size=25, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))

pdf("propoffourth3.pdf")
p
dev.off()
###########################################################################################