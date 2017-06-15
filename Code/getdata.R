getdata<-function(numb){
  
  if (numb==1){
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
  }
  
  if (numb==2){
    #Data Munging:
    #download data
    myjson<-fromJSON("https://melt.firebaseio.com/.json")
    
    k<-1
    data<-data.frame(id= numeric(), x1=numeric(), x2=numeric(), x3=numeric(), chosen=numeric(), y=numeric())
    for (i in 1:59){
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
    
    #values bigger than 100 were presented as 100
    data$y<-ifelse(data$y>100, 100, data$y)
    #values smaller than 0 were presented as 0
    data$y<-ifelse(data$y<0, 0, data$y)
    #standardize to be between -1 and 1
    data$y<-(data$y-50)/100
    #trial number
    data$trial<-rep(1:150, nrow(data)/150)
  }
  
  if (numb==3){
    myjson<-fromJSON("https://melt.firebaseio.com/.json")
    
    k<-1
    data<-data.frame(id= numeric(), x1=numeric(), x2=numeric(), x3=numeric(), y=numeric())
    for (i in 62:length(myjson)){
      d<-myjson[[i]]
      x1<-d$sunseen 
      x2<-d$tempseen 
      x3<-d$rainseen
      y<-d$reward
      chosen<-d$chosendeck
      #dummy frame
      dummy<-data.frame(id=rep(k, length(x1)), x1=x1, x2=x2, x3=x3, chosen=chosen, y=y)
      #concatenate
      data<-rbind(data, dummy)
      #k increment
      k<-k+1
    }
    
    data<-subset(data, id!=3)
    data$id<-rep(1:(nrow(data)/150), each=150)
    #values bigger than 100 were presented as 100
    data$y<-ifelse(data$y>100, 100, data$y)
    #values smaller than 0 were presented as 0
    data$y<-ifelse(data$y<0, 0, data$y)
    #standardize to be between -1 and 1
    data$y<-(data$y-50)/100
    #trial number
    data$trial<-rep(1:150, nrow(data)/150)
  }
  
return(data)
}