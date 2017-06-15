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
data<-data.frame(id= numeric(), x1=numeric(), x2=numeric(), x3=numeric(), chosen=numeric(), y=numeric())
for (i in 1:61){
  d<-myjson[[k]]
  x1<-d$sunseen 
  x2<-d$tempseen 
  x3<-d$rainseen
  y<-d$reward
  chosen<-d$chosendeck
  #dummy frame
  dummy<-data.frame(id=rep(i, length(x1)), x1=x1, x2=x2, x3=x3, chosen=chosen, y=y)
  #concatenate
  data<-rbind(data, dummy)
  #k increment
  k<-k+1
}

names(data)

tapply(data$x1, data$id, length)
data<-subset(data, id<60)
indicator<-tapply(data$y, data$id, function(x){t.test(x[1:150], mu=50,  conf.level = 0.9)$conf.int[1]} )

sum(indicator<45)

indicator<-c(1:59)[indicator>45]
data<-subset(data, id %in% indicator)



dplot1<-ddply(data,~id,summarise, mu=mean(y))

p<-ggplot(dplot1, aes(mu)) +
  geom_histogram(bins=10, col="grey50")+
  ggtitle("Average Scores")+theme_classic() +xlab("\nMean Score")+ylab("Counts\n")+
  #fill with Wes Anderson colors
  theme(text = element_text(size=25, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))

p

pdf("/home/hanshalbe/Desktop/cmab/paper/figures/score2.pdf")
p
dev.off()

genoutput<-function(x,y,z){
  p1<-50+x*3+y*-3
  p2<-50+x*-3+z*3
  p3<-50+y*3+z*-3
  p4<-rep(50, length(x))
  comb<-cbind(p1,p2,p3,p4)
  return(apply(comb,1,which.max))
}

genoutput(data$x1, data$x2, data$x3)

data$trials<-rep(1:150, max(data$id)-2)
dplot<-ddply(data, ~trials, summarize, prop=mean(y))

cor(dplot2$trials, dplot2$prop)
p<-ggplot(dplot, aes(x=trials, y=prop)) + 
  geom_point(col="grey",size=3)+
  stat_smooth(method="lm", col="black")+
  geom_errorbar(data=dplot2, aes(x=tens, y=mu, ymin=mu-se, ymax=mu+se), width=1.2, size=1, position=pd, col="red") +
  xlab("\nTrial")+ylab("Mean\n")+
  ggtitle("Outcomes over time")+
  theme_classic()+theme(text = element_text(size=22))+ ylim(c(45,65))+
  annotate("text", x = 40, y = 65, label = "r=0.39, p<0.01", size=8)+
  #fill with Wes Anderson colors
  theme(text = element_text(size=25, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))
p


pdf("meanovertime2.pdf")
print(p)
dev.off()

data$isfour<-ifelse(data$chosen==4,1,0)
dplot2<-ddply(data, ~trials, summarize, prop=mean(isfour))
dplot2$prop[1]<-0.25
cor(dplot2$prop, dplot2$trials)
p<-ggplot(dplot2, aes(x=trials, y=prop)) + 
  geom_point(col="grey",size=3)+
  stat_smooth(method="lm", col="black")+
  xlab("\nTrial")+ylab("Mean\n")+
  ggtitle("Outcomes over time")+
  theme_classic()+theme(text = element_text(size=22))+ ylim(c(0,0.4))+
  annotate("text", x = 40, y = 0.35, label = "r=0.-0.22, p<0.05", size=8)+
  #fill with Wes Anderson colors
  theme(text = element_text(size=25, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))

pdf("propfour.pdf")
p
dev.off()
###########################################################################################