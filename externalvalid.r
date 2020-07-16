
###Install and Load Packages
packages <- c("ranger", "ggplot2", "qs", "vtreat", "plotrix", "Amelia",
              "CORElearn","rms","timeROC","scales","pec","ggthemes","dplyr","plyr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

{
library(ranger)
library(ggplot2)
library(qs)
library(plyr)
library(dplyr)
library(vtreat)
library(plotrix)
library(Amelia)
library(CORElearn)
library(rms)
library(timeROC)
library(scales)
library(pec)
library(ggthemes)
}

list<-qread('modellist.q')
rsflist<-list$rsflist
dummies<-list$dummies
testDF<-list$testDF
###Ensure that testing data has same variables and factor levels as displayed below; then merge into this dataframe
str(testDF)

###Label data to be tested as testingDATA
testingDATA

###Remove cases with missing values from test set 
testingDATA<-na.omit(testingDATA)

##Censor data at 60 months and round up 
{
testingDATA$monthssurv<-ceiling(testingDATA$monthssurv)

for (i in 1:nrow(testingDATA)){
  ifelse(testingDATA$monthssurv[i]>=60, testingDATA$ons_death[i]<-0,testingDATA$ons_death[i]<-testingDATA$ons_death[i])
}
for (i in 1:nrow(testingDATA)){
  ifelse(testingDATA$monthssurv[i]>=60, testingDATA$monthssurv[i]<-60,testingDATA$monthssurv[i]<-testingDATA$monthssurv[i])
}
}

###Calculate predictions, can be time consuming (around 10-20s for each 100 cases)

system.time({
  data<-vtreat::prepare(dummies,testingDATA,pruneSig=c())
  estimates<-1
  survfull<-1
  for (i in 1:10){
    estimates[i]<-list(log(predict(rsflist[[i]],data,type='response')$survival))
  }
  for (i in 1:10){
  survfull[i]<-list(log(predict(rsflist[[i]],data,type='response',predict.all=TRUE)$survival))
  }
  standerx1<-1
  standerx<-as.data.frame(1:nrow(testingDATA))
  for (i in 1:10){
  for (p in 1:60){
    for (q in 1:nrow(testingDATA)){
      standerx[q,p]<-std.error(survfull[[i]][q,p,c(1:rsflist[[i]]$num.trees)])
      standerx1[i]<-list(standerx)
    }  
  }
  }
mipreds<-as.data.frame(1:nrow(data))
for (m in 1:60){
  mipreds[,m]<-as.data.frame(t(as.data.frame(mi.meld(q=cbind(estimates[[1]][,m],estimates[[2]][,m],estimates[[3]][,m],estimates[[5]][,m],
                                                             estimates[[6]][,m],estimates[[7]][,m],estimates[[8]][,m],estimates[[9]][,m],estimates[[10]][,m]),
                                                     se=cbind(standerx1[[1]][,m],standerx1[[2]][,m],standerx1[[3]][,m],standerx1[[5]][,m],
                                                              standerx1[[6]][,m],standerx1[[7]][,m],standerx1[[8]][,m],standerx1[[9]][,m],standerx1[[10]][,m]),byrow=FALSE)[1])))
}
mipreds2<-exp(mipreds)
colnames(mipreds2)<-c(1:60)
finalpredictions<-mipreds2
})

####Dataframe of final predicted probability of survival at time 1:60 months
finalpredictions
  
###Plot 60 month Time dependent ROC Curve
{
  troc<-timeROC(testingDATA$monthssurv,testingDATA$ons_death,1-finalpredictions[,59],times=59,cause=1,weighting='marginal',iid=TRUE)
  xx<-as.data.frame(1-troc$TP)
  xx$Spec<-1-troc$FP[,2]
  xx$Sens<-1-xx$`t=59`
  timeROC<-ggplot(xx,aes(Spec,Sens))+
    geom_line()+
    xlim(1,0)+
    xlab('Specificity')+
    ylab('Sensitivity')+
    ggtitle('Time Dependent ROC at 60 months')+
    geom_abline(intercept=1,linetype=3)+
    annotation_custom(grid::textGrob(paste("AUC = ",percent(troc$AUC[2]))), 
                      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
  
  
  timeROC
}

timeROC

###Calculate Harrell's C-Index (discrimination) 
{  
surv.obj=with(data,Surv(testingDATA$monthssurv,testingDATA$ons_death))
CIndex<-rcorr.cens(x=finalpredictions[,60],S=surv.obj)[1]
}

CIndex

###Calculate Integrated Brier Score
{
new<-as.data.frame(1:nrow(testingDATA))
new$`1:nrow(testingDATA)`<-1
finalpredictions1<-as.data.frame(c(new,finalpredictions))
pec<-pec(list(predictions=as.matrix(finalpredictions1)),Hist(monthssurv,ons_death)~1,data=testingDATA)
ibrier<-crps(pec)
ibrier
}

ibrier

###Plot Calibration charts from 1-5 years

{
r12<-calPlot(finalpredictions[,12],time=12,formula=Surv(monthssurv,ons_death)~.,data=testingDATA,type='survival',legend.legend=c("12 months"),plot=FALSE)
r24<-calPlot(finalpredictions[,24],time=24,formula=Surv(monthssurv,ons_death)~.,data=testingDATA,type='survival',legend.legend=c("24 months"),plot=FALSE)
r36<-calPlot(finalpredictions[,36],time=36,formula=Surv(monthssurv,ons_death)~.,data=testingDATA,type='survival',legend.legend=c("36 months"),plot=FALSE)
r48<-calPlot(finalpredictions[,48],time=48,formula=Surv(monthssurv,ons_death)~.,data=testingDATA,type='survival',legend.legend=c("48 months"),plot=FALSE)
r60<-calPlot(finalpredictions[,60],time=60,formula=Surv(monthssurv,ons_death)~.,data=testingDATA,type='survival',legend.legend=c("60 months"),plot=FALSE)
}

plot(r12)
plot(r24)
plot(r36)
plot(r48)
plot(r60)
  
####Quintile Calibration Plot
{

  calib<-finalpredictions
  calib$id<-1:nrow(calib)
  calib$death<-testingDATA$ons_death
  calib$months<-testingDATA$monthssurv
  finalpredictions$`60`
  calib$decile<-with(finalpredictions,cut(`60`,
                                  breaks = quantile(`60`,probs=seq(0,1,by=0.2)),
                                  include.lowest = TRUE))
  levels(calib$decile)<-c('0-20','20-40','40-60','60-80','80-100')
  ###calib$decile<-with(calib,cut(V1.50,breaks = quantile(V1.50,probs=seq(0,1,by=0.25)),include.lowest = TRUE))
  ###levels(calib$decile)<-c('0-10','10-20','20-30','30-40')
  count(calib$decile)
  mts<-as.data.frame(cbind(testingDATA$monthssurv,testingDATA$ons_death))
  names(mts)<-c('monthssurv','ons_death')
  predscalib2<-finalpredictions
  calib1<-finalpredictions
  calib1$ons_death<-testingDATA$ons_death
  calib1$monthssurv<-testingDATA$monthssurv
  
  {
    predscalib2$decile<-calib$decile
    predscalib2$death<-calib1$ons_death
    predscalib2$months<-calib1$monthssurv
    dec1<-subset(predscalib2,predscalib2$decile=='0-20')
    dec1a<-select(dec1,-decile,-months,-death)
    estimatesdec1<-as.data.frame(1:60)
    estimatesdec1$Time<-estimatesdec1[,1]
    estimatesdec1$survival<-colMeans(dec1a)
    OSKMa <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
    OSKM2a<-as.data.frame(cbind(OSKMa$time,OSKMa$surv))
    names(OSKM2a)<-c('time','surv')  
    
    predscalib2$decile<-calib$decile
    predscalib2$death<-calib1$ons_death
    predscalib2$months<-calib1$monthssurv
    dec1<-subset(predscalib2,predscalib2$decile=='20-40')
    dec1a<-select(dec1,-decile,-months,-death)
    estimatesdec1b<-as.data.frame(1:60)
    estimatesdec1b$Time<-estimatesdec1b[,1]
    estimatesdec1b$survival<-colMeans(dec1a)
    OSKMb <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
    OSKM2b<-as.data.frame(cbind(OSKMb$time,OSKMb$surv))
    names(OSKM2b)<-c('time','surv')  
    
    predscalib2$decile<-calib$decile
    predscalib2$death<-calib1$ons_death
    predscalib2$months<-calib1$monthssurv
    dec1<-subset(predscalib2,predscalib2$decile=='40-60')
    dec1a<-select(dec1,-decile,-months,-death)
    estimatesdec1c<-as.data.frame(1:60)
    estimatesdec1c$Time<-estimatesdec1c[,1]
    estimatesdec1c$survival<-colMeans(dec1a)
    OSKMc <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
    OSKM2c<-as.data.frame(cbind(OSKMc$time,OSKMc$surv))
    names(OSKM2c)<-c('time','surv')  
    
    predscalib2$decile<-calib$decile
    predscalib2$death<-calib1$ons_death
    predscalib2$months<-calib1$monthssurv
    dec1<-subset(predscalib2,predscalib2$decile=='60-80')
    dec1a<-select(dec1,-decile,-months,-death)
    estimatesdec1d<-as.data.frame(1:60)
    estimatesdec1d$Time<-estimatesdec1d[,1]
    estimatesdec1d$survival<-colMeans(dec1a)
    OSKMd <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
    OSKM2d<-as.data.frame(cbind(OSKMd$time,OSKMd$surv))
    names(OSKM2d)<-c('time','surv') 
    
    predscalib2$decile<-calib$decile
    predscalib2$death<-calib1$ons_death
    predscalib2$months<-calib1$monthssurv
    dec1<-subset(predscalib2,predscalib2$decile=='80-100')
    dec1a<-select(dec1,-decile,-months,-death)
    estimatesdec1e<-as.data.frame(1:60)
    estimatesdec1e$Time<-estimatesdec1e[,1]
    estimatesdec1e$survival<-colMeans(dec1a)
    OSKMe <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
    OSKM2e<-as.data.frame(cbind(OSKMe$time,OSKMe$surv))
    names(OSKM2e)<-c('time','surv') 
    estimatesdec1$col='Predicted'
    estimatesdec1b$col='Predicted'
    estimatesdec1c$col='Predicted'
    estimatesdec1d$col='Predicted'
    estimatesdec1e$col='Predicted'
    estimatesdec1$lt='1'
    estimatesdec1b$lt='2'
    estimatesdec1c$lt='3'
    estimatesdec1d$lt='4'
    estimatesdec1e$lt='5'
    estimatesdec1<-select(estimatesdec1,-1)
    estimatesdec1b<-select(estimatesdec1b,-1)
    estimatesdec1c<-select(estimatesdec1c,-1)
    estimatesdec1d<-select(estimatesdec1d,-1)
    estimatesdec1e<-select(estimatesdec1e,-1)
    colnames(OSKM2a)<-c('Time','survival')
    colnames(OSKM2b)<-c('Time','survival')
    colnames(OSKM2c)<-c('Time','survival')
    colnames(OSKM2d)<-c('Time','survival')
    colnames(OSKM2e)<-c('Time','survival')
    OSKM2a$col='Observed'
    OSKM2b$col='Observed'
    OSKM2c$col='Observed'
    OSKM2d$col='Observed'
    OSKM2e$col='Observed'
    OSKM2a$lt='1'
    OSKM2b$lt='2'
    OSKM2c$lt='3'
    OSKM2d$lt='4'
    OSKM2e$lt='5'
    colnames(estimatesdec1)
    colnames(OSKM2a)
    full<-rbind(estimatesdec1,estimatesdec1b,estimatesdec1c,estimatesdec1d,estimatesdec1e,
                OSKM2a,OSKM2b,OSKM2c,OSKM2d,OSKM2e)

    quintileplot<-ggplot(full,aes(Time,survival))+
      geom_line(data=full,mapping=aes(colour=col,linetype=lt))+
      xlim(0,60)+ylim(0,1)+ylab('Survival')+
      scale_color_discrete(name='')+
      scale_linetype_discrete(name='Probability Quintile',
                              breaks=c('5','4','3','2','1'),
                              labels=c('80-100','60-80','40-60','20-40','0-20'))+
      ggtitle('Observed Survival by Probability Quintiles')
    
    
  }
}

quintileplot

