#House keeping
rm(list=ls())

#required packages packages
packages <- c('RCurl', 'devtools', 'ggplot2', 'parallel', 'plyr')
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
data
names(data)


indicator<-tapply(data$y, data$id, function(x){t.test(x, mu=50)$conf.int[[1]]>50})

dplot1<-ddply(data,~id,summarise, mu=mean(y))

p1<-ggplot(dplot1, aes(mu)) +
  geom_histogram(bins=10, col="grey50")+
  ggtitle("Average Scores")+theme_classic() +xlab("\nMean Score")+ylab("Counts\n")+
    #fill with Wes Anderson colors
  theme(text = element_text(size=25, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))


p1

pdf("/home/hanshalbe/Desktop/cmab/paper/figures/score1.pdf")
p
dev.off()

head(data)
data$isfour<-ifelse(data$chosen==4,1,0)
data$trials<-rep(1:150, max(data$id))
dplot2<-ddply(data, ~trials, summarize, prop=mean(isfour))
dplot2$prop[1]<-0.25
p<-ggplot(dplot2, aes(x=trials, y=prop)) + 
  geom_point(col="grey",size=3)+
  stat_smooth(method="lm", col="black")+
  xlab("\nTrial")+ylab("Mean\n")+
  ggtitle("Proportion of 4th arm")+
  theme_classic()+theme(text = element_text(size=22))+ ylim(c(0,0.4))+
  annotate("text", x = 40, y = 0.35, label = "r=0.-0.22, p<0.05", size=8)+
  #fill with Wes Anderson colors
  theme(text = element_text(size=25, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))

pdf("propfour.pdf")
p
dev.off()

genoutput<-function(x,y,z){
  p1<-50+x*20+y*-20
  p2<-50+x*-20+z*20
  p3<-50+y*20+z*-20
  p4<-rep(50, length(x))
  comb<-cbind(p1,p2,p3,p4)
  return(apply(comb,1,which.max))
}


data$best<-ifelse(genoutput(data$x1, data$x2, data$x3)==data$chosen,1,0)
dplot2<-ddply(data, ~trials, summarize, prop=mean(y))

cor(dplot2$trials, dplot2$prop)
p<-ggplot(dplot2, aes(x=trials, y=prop)) + 
  geom_point(col="grey",size=3)+
  stat_smooth(method="lm", col="black")+
  xlab("\nTrial")+ylab("Mean\n")+
  ggtitle("Outcomes over time")+
  theme_classic()+theme(text = element_text(size=22))+ ylim(c(40,80))+
  annotate("text", x = 40, y = 80, label = "r=0.74, p<0.01", size=8)+
  #fill with Wes Anderson colors
  theme(text = element_text(size=25, family="serif"),
        strip.background = element_blank(),
        plot.title = element_text(size=26))
p

pdf("/home/hanshalbe/Desktop/cmab/paper/figures/meanovertime1.pdf")
print(p)
dev.off()
###########################################################################################