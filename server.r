histdata <- 1
{
library(shiny)
library(scales)
library(shinyBS)
library(plotly)
library(lattice)
library(plyr)
library(shinyWidgets)
library(caret)
library(xgboost)
library(glmnet)
library(survival)
library(rms)
library(MKmisc)
library(iml)
library(formattable)
library(givitiR) 
library(pec)
library(ranger)
library(vtreat)
library(gbm)
library(pROC)
library(personograph)
library(caretEnsemble)
library(e1071)
library(qs)
library(ranger)
library(ggplot2)
library(vtreat)
library(plotrix)
library(Amelia)
library(CORElearn)
library(rms)
library(tableone)
  library(summarytools)
  library(ggthemes)
  library(scales)
  
}
list<-qread('modellist.q')
rsflist<-list$rsflist
finalvar<-list$finalvar
dummies<-list$dummies
deathtimes<-list$deathtimes
cab<-list$cab
selvar<-qread('selvar.q')
sumvar<-qread('sumvar.q')
dfsumm<-qread('dfsumm.q')
s<-qread('s.q')
s1<-qread('s1.q')
calib<-qread('calibrsf2.q')
censordata<-qread('censordata2.q')
condplot<-qread('condplot.q')
km<-qread('kmgps.q')
oskm<-qread('oskmstg.q')
quantiles<-c(0.783,0.583,0.28,0.131)
predict_fct = function(model, data){
  data<-vtreat::prepare(dummies,data,pruneSig=c())
  estimates<-1
  survfull<-1
  for (i in 1:10){
    estimates[i]<-list(as.data.frame(log(predict(model[[i]],data,type='response')$survival)))
  }
  for (i in 1:10){
    survfull[i]<-list(log(predict(model[[i]],data,type='response',predict.all=TRUE,threads=3)$survival))
  }
  standerx1<-1
  standerx<-as.data.frame(1:nrow(data))
  for (i in 1:10){
    for (p in 1:60){
      for (q in 1:nrow(data)){
        standerx[q,p]<-std.error(survfull[[i]][q,p,c(1:rsflist[[i]]$num.trees)])
        standerx1[i]<-list(as.data.frame(standerx))
      }  
    }
  }
  mipreds<-as.data.frame(1:nrow(data))
  if(nrow(mipreds)>1){
    for (u in 1:60){
      mipreds[,u]<-as.data.frame(t(as.data.frame(mi.meld(q=cbind(estimates[[1]][,u],estimates[[2]][,u],estimates[[3]][,u],estimates[[5]][,u],
                                                                 estimates[[6]][,u],estimates[[7]][,u],estimates[[8]][,u],estimates[[9]][,u],estimates[[10]][,u]),
                                                         se=cbind(standerx1[[1]][,u],standerx1[[2]][,u],standerx1[[3]][,u],standerx1[[5]][,u],
                                                                  standerx1[[6]][,u],standerx1[[7]][,u],standerx1[[8]][,u],standerx1[[9]][,u],standerx1[[10]][,u]),byrow=FALSE)[1])))
    }
  } else
    for (u in 1:60){
      mipreds[,u]<-as.data.frame(t(as.data.frame(mi.meld(q=cbind(estimates[[1]][u,],estimates[[2]][u,],estimates[[3]][u,],estimates[[5]][u,],
                                                                 estimates[[6]][u,],estimates[[7]][u,],estimates[[8]][u,],estimates[[9]][u,],estimates[[10]][u,]),
                                                         se=cbind(standerx1[[1]][,u],standerx1[[2]][,u],standerx1[[3]][,u],standerx1[[5]][,u],
                                                                  standerx1[[6]][,u],standerx1[[7]][,u],standerx1[[8]][,u],standerx1[[9]][,u],standerx1[[10]][,u]),byrow=FALSE)[1])))
    }
  mipreds2<-exp(mipreds)
  predscalib<-mipreds2
}
stagingfunction<-function(data){
  data$pT.ypT<-as.numeric(substring(paste(data$pT.ypT),2))
  data$pM.ypM<-0
  stage<-ifelse(data$pT.ypT==0&data$Total.LN.positive==0,'0',
                ifelse(data$Neoadjuvant.Treatment!='None'& data$pT.ypT<3 & data$Total.LN.positive==0 & data$pM.ypM==0,'1',
                       ifelse(data$Neoadjuvant.Treatment!='None'& data$pT.ypT==3 & data$Total.LN.positive==0 & data$pM.ypM==0,'2',
                              ifelse(data$Neoadjuvant.Treatment!='None'& data$pT.ypT<3 & data$Total.LN.positive<3 & data$pM.ypM==0,'3a',
                                     ifelse(data$Neoadjuvant.Treatment!='None'& data$pT.ypT==4 & data$Total.LN.positive==0 & data$pM.ypM==0,'3b',
                                            ifelse(data$Neoadjuvant.Treatment!='None'& data$Total.LN.positive<7 & data$pM.ypM==0,'3b',
                                                   ifelse(data$Neoadjuvant.Treatment!='None'& data$pM.ypM==0, '4a',
                                                          ifelse(data$Neoadjuvant.Treatment!='None'& data$pM.ypM==1,'4b',
                                                                 ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT<2&data$Total.LN.positive==0,'1',
                                                                        ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==2&data$Total.LN.positive==0&data$Grade=='G1','1',
                                                                               ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==2&data$Total.LN.positive==0&data$Grade=='G2','1',
                                                                                      ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==2&data$Total.LN.positive==0,'2',
                                                                                             ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==1&data$Total.LN.positive<3,'2',
                                                                                                    ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==3&data$Total.LN.positive==0,'2',
                                                                                                           ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==1&data$Total.LN.positive<7,'3a',
                                                                                                                  ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==2&data$Total.LN.positive<3,'3a',
                                                                                                                         ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==2&data$Total.LN.positive<7,'3b',
                                                                                                                                ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==3&data$Total.LN.positive<7,'3b',
                                                                                                                                       ifelse(data$Neoadjuvant.Treatment=='None'&data$pT.ypT==4&data$Total.LN.positive<3,'3b','4a'
                                                                                                                                       )))))))))))))))))))
  stage
}
library(mice)

bsModalNoClose <-function(...) {
  b = bsModal(...)
  b[[2]]$`data-backdrop` = "static"
  b[[2]]$`data-keyboard` = "false"
  return(b)
}



shinyServer
(
  function(input, output)
  {
    
    
    output$downloadData <- downloadHandler(
      filename <- function() {
        paste("modellist", "q", sep=".")
      },
      
      content <- function(file) {
        file.copy("modellist.q", file)
      },
      contentType = "application/zip"
    )
    percent <- function(x, digits = 1, format = "f", ...) {
      paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
    }
    inputdata <- reactive({
      xcomp<-ifelse(input$Any.Complication=='Yes',1,0)
      xhist<-ifelse(input$Histopathology=='SCC','SCC','Adenocarcinoma')
      data0 <- data.frame(
        Age = as.numeric(input$Age),
        Histopathology = as.factor(xhist),
        'Neoadjuvant.Treatment' = as.factor(input$Neoadjuvant.Treatment),
        'Completion.of.Neoadjuvant.Treatment'= as.factor(input$Completion.of.Neoadjuvant.Treatment),
        'Involved.Longitudinal Margin' = as.factor(input$Long),
        'Involved.Circumferential Margin' = as.factor(input$Circ),
        'Total.LN.positive' = as.numeric(input$Total.LN.Positive),
        'Any.Complication' = as.factor(xcomp),
        'Site.of.Tumour' = as.factor(input$Site.of.Tumour),
        'Gender' = as.factor(input$Gender),
        'Grade'=as.factor(input$Grade),
        pT.ypT = as.ordered(input$pT.ypT),
        cT = as.ordered(input$cT),
        cN = as.ordered(input$cN)
      )
predscalib<-predict_fct(rsflist,data0)
predscalib
    })
    inputdata2 <- reactive({
      xcomp<-ifelse(input$Any.Complication=='Yes',1,0)
      xhist<-ifelse(input$Histopathology=='SCC','SCC','Adenocarcinoma')
      data0 <- data.frame(
        Age = as.numeric(input$Age),
        Histopathology = as.factor(xhist),
        'Neoadjuvant.Treatment' = as.factor(input$Neoadjuvant.Treatment),
        'Completion.of.Neoadjuvant.Treatment'= as.factor(input$Completion.of.Neoadjuvant.Treatment),
        'Involved.Longitudinal Margin' = as.factor(input$Long),
        'Involved.Circumferential Margin' = as.factor(input$Circ),
        'Total.LN.positive' = as.numeric(input$Total.LN.Positive),
        'Any.Complication' = as.factor(xcomp),
        'Site.of.Tumour' = as.factor(input$Site.of.Tumour),
        'Gender' = as.factor(input$Gender),
        'Grade'=as.factor(input$Grade),
        pT.ypT = as.ordered(input$pT.ypT),
        cT = as.ordered(input$cT),
        cN = as.ordered(input$cN)
      )
      data0
    })
    output$result <- renderTable({
      
      predscalib = inputdata()
      
      {
        pred1yrpct<-predscalib[,12]
        pred2yrpct<-predscalib[,24]
        pred5yrpct<-predscalib[,60]
        predictionstab<-rbind(pred1yrpct,pred2yrpct,pred5yrpct)
        predictionstab1<-as.data.frame(percent(predictionstab))
        rownames(predictionstab1)<-c("Alive at 12 months",'Alive at 24 months','Alive at 60 months')
        colnames(predictionstab1)<-c('Probability')
        }
      resultTable = (predictionstab1)
      resultTable
    },rownames = TRUE)
    plotInput <- reactive({
      predscalib = inputdata()
      data=inputdata2()
      {
        x3<-stagingfunction(data)
        stg<-ifelse(data$Neoadjuvant.Treatment=='None',paste('pTNM',x3),
                    ifelse(!is.na(data$Neoadjuvant.Treatment),paste('ypTNM',x3)))
        gp<-ifelse((x3==0|x3==1)&predscalib[,60]>quantiles[1],'Low',
                   ifelse((x3==0|x3==1)&predscalib[,60]<=quantiles[1],'High',
                          ifelse((x3==2|x3=='3a')&predscalib[,60]>quantiles[2],'Low',
                                 ifelse((x3==2|x3=='3a')&predscalib[,60]<=quantiles[2],'High',
                                        ifelse((x3=='3b')&predscalib[,60]>quantiles[3],'Low',
                                               ifelse((x3=='3b')&predscalib[,60]<=quantiles[3],'High',
                                                      ifelse((x3=='4a')&predscalib[,60]>quantiles[4],'Low',
                                                             ifelse((x3=='4a')&predscalib[,60]<=quantiles[4],'High',NA))))))))
        x4<-x3
        x4[x4==0]<-'0-1'
        x4[x4==1]<-'0-1'
        x4[x4==2]<-'2-3a'
        x4[x4=='3a']<-'2-3a'
        fingp<-paste('TNM',x4,gp)
        survs<-as.data.frame(1)
        ifelse(x3==0,survs<-cbind(oskm$surv[1:76],oskm$time[1:76]),
               ifelse(x3==1,survs<-cbind(oskm$surv[77:157],oskm$time[77:157]),
                      ifelse(x3==2,survs<-cbind(oskm$surv[158:238],oskm$time[158:238]),
                             ifelse(x3=='3a',survs<-cbind(oskm$surv[239:313],oskm$time[239:313]),
                                    ifelse(x3=='3b',survs<-cbind(oskm$surv[314:394],oskm$time[314:394]),
                                           ifelse(x3=='4a',survs<-cbind(oskm$surv[395:463],oskm$time[395:463]),NA))))))
        
        
        ifelse((x3==0|x3==1)&gp=='High',survs2<-cbind(km$surv[1:60],km$time[1:60]),
               ifelse((x3==2|x3=='3a')&gp=='High',survs2<-cbind(km$surv[61:119],km$time[61:119]),
                      ifelse((x3=='3b')&gp=='High',survs2<-cbind(km$surv[120:177],km$time[120:177]),
                             ifelse((x3=='4a')&gp=='High',survs2<-cbind(km$surv[178:220],km$time[178:220]),
                                    ifelse((x3==0|x3==1)&gp=='Low',survs2<-cbind(km$surv[221:279],km$time[221:279]),
                                           ifelse((x3==2|x3=='3a')&gp=='Low',survs2<-cbind(km$surv[280:338],km$time[280:338]),
                                                  ifelse((x3=='3b')&gp=='Low',survs2<-cbind(km$surv[339:398],km$time[339:398]),
                                                         ifelse((x3=='4a')&gp=='Low',survs2<-cbind(km$surv[399:454],km$time[399:454])))))))))
        
        survs<-as.data.frame(survs)
        colnames(survs)<-c('Surv','Time')
        survs$Type<-paste('TNM Stage Mean','(',stg,')',sep='')
        survs2<-as.data.frame(survs2)
        colnames(survs2)<-c('Surv','Time')
        survs2$Type<-paste('Risk Group Mean','(',fingp,')',sep='')
        preds<-as.data.frame(t(predscalib))
        colnames(preds)<-'Surv'
        preds$Time<-1:60
        preds$Type<-'Individual Prediction'
        sum<-rbind(preds,survs,survs2)
        sum
      }
      {
        survplot<-ggplot(sum,aes(Time,Surv))+
          geom_smooth(se=FALSE,aes(col=Type),method='loess',span=0.05)+
          ylim(0,1)+
          scale_x_continuous(breaks=c(0,12,24,36,48,60),expand=c(0,0),limits=c(3,60))+
          xlab('Months After Surgery')+
          ylab('Probability of Survival')+
          ###ggtitle(paste('Stage --',stg, '/ Risk Group --', gp))+theme(legend.title = element_blank())
        ggtitle(paste('Stage --',stg))+theme(legend.title = element_blank())
        survplot
        ggs<-ggplotly(survplot)
        ggs$x$data[[3]]$visible<-'legendonly'
        ggs$x$data[[2]]$visible<-FALSE
        ggs
      }
    })
    output$plot<-renderPlotly({
      print(plotInput())
    })
    output$person1yr<-renderPlot({
      predscalib = inputdata()
      pred1yr<-(predscalib[,12])
      alive1yr<-round(pred1yr,digits=2)
      dead1yr<-round(1-pred1yr,digits=2)
      data1yr<-list(Alive=alive1yr,Dead=dead1yr)
      personograph(data1yr,colors=list(Alive='limegreen',Dead='orangered3'),fig.title="1 year",icon.style=11)
    })
    output$person2yr<-renderPlot({
      predscalib = inputdata()
      
      pred2yr<-predscalib[,24]
      alive2yr<-round(pred2yr,digits=2)
      dead2yr<-round(1-pred2yr,digits=2)
      data2yr<-list(Alive=alive2yr,Dead=dead2yr)
      personograph(data2yr,colors=list(Alive='limegreen',Dead='orangered3'),fig.title="2 years",icon.style=11)
    })
    output$person5yr<-renderPlot({
      predscalib= inputdata()
      
      pred5yr<-predscalib[,60]
      alive5yr<-round(pred5yr,digits=2)
      dead5yr<-round(1-pred5yr,digits=2)
      data5yr<-list(Alive=alive5yr,Dead=dead5yr)
      personograph(data5yr,colors=list(Alive='limegreen',Dead='orangered3'),fig.title="5 years",icon.style=11)
    })
    output$fullvar<-renderTable({
      as.data.frame(selvar)
    })
    output$usevar<-renderTable({
      as.data.frame(finalvar)
    })
    output$table1<-renderUI({
      print(dfsumm,method='render',bootstrap.css=FALSE,varnumbers=FALSE,valid.col=FALSE)
    })
    output$Completion<-renderUI({
      if (input$Neoadjuvant.Treatment!='None')
      selectInput("Completion.of.Neoadjuvant.Treatment", "Completion", choices = c("Completed", "Not Completed"),selected = 'Completed') else
        selectInput("Completion.of.Neoadjuvant.Treatment", "Completion", choices = c("Not Applicable"),selected = 'Not Applicable')
    })
    
    
    output$plotly1<-renderPlotly({
      axx <- list(
        title = "Time (months)"
      )
      axy <- list(
        title = "Positive Lymph Nodes"
      )
      axy1 <- list(
        title = "Age"
      )
      axz <- list(
        title = "Survival (%)",
        range=c(0,100)
      )
      p<-plot_ly(x= s$y,y=s$x,z=s$z) %>% add_surface(color = ~s$x, colorscale = 'Jet', showscale = TRUE) %>%
        layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
      p
    })
    output$plotly2<-renderPlotly({
      axx <- list(
        title = "Time (months)"
      )
      axy <- list(
        title = "Positive Lymph Nodes"
      )
      axy1 <- list(
        title = "Age"
      )
      axz <- list(
        title = "Survival (%)",
        range=c(0,100)
      )
      p1<-plot_ly(x= s1$y,y=s1$x,z=s1$z) %>% add_surface(color = ~s1$x, colorscale = 'Jet', showscale = TRUE) %>%
        layout(scene = list(xaxis=axx,yaxis=axy1,zaxis=axz))
      p1
    })
    dateinput<-reactive({
      dates<-c(input$dateapp,input$datesurg)
      dates
    })
    output$condsurv<-renderPlotly({
      predscalib<-inputdata()
      dates<-dateinput()
      time<-ceiling(as.numeric(((dates[1]-dates[2])/365)*12))
      {
        z<-c(time:60)
        z2<-as.data.frame(1:nrow(predscalib))
        for (i in time:(max(z)-1)){
          z2[(i-time+1)]<-predscalib[,i+1]/predscalib[,time]
        }
        colnames(z2)<-c((time+1):(max(z)))
        z2[z2>1]<-1
        z2[is.na(z2)]<-1
        z3<-select(z2,1:((60-time)))
        z3[,(60-time)+1]<-1
        colnames(z3)<-c(colnames(z3)[1:(60-time)],paste(time))
        z3<-z3[,c(ncol(z3),1:(60-time))]
      }
      {
        gg<-(seq(Sys.Date(),Sys.Date()+(((60-time)/12)*365),by=30))
        survplt<-as.data.frame(gg)
        survplt$survival<-(t(z3))
        names(survplt)<-c('time','survival')
        survplt
        survplot<-ggplot(survplt,aes(time,survival))+geom_smooth(colour='red',method='loess',se=FALSE)+
          ylim(0,1)+
          scale_x_date(breaks=c(seq(Sys.Date(),Sys.Date()+((60-time)/12)*365,by=round(0.5*365))),date_labels = "%b %y")+
          xlab('Date')+
          ylab('Probability of Survival')+
          ggtitle(c(paste('Predicted survival conditional on',time,'months prior survival')))+theme_bw()+
          geom_vline(xintercept=Sys.Date()+(((60-time)/12)*365))+
          annotate(geom="text", label=c(paste('5 years','post surgery',sep="\n")), x=Sys.Date()+((60-time)/12)*365, y=0, vjust=0,hjust=1)
        survplot2<-ggplotly(survplot,tooltip=c('survival'))
        survplot2
      }
      
  
    })
    output$condcalib<-renderPlot({
      predscalib<-inputdata()
      dates<-dateinput()
      time<-ceiling(as.numeric(((dates[1]-dates[2])/365)*12))
      plot<-condplot[[time]]
      plot
    })
      output$varimp1yr<-renderTable({
        as.data.frame(year1varimp)
      },rownames=TRUE)
      output$varimp2yr<-renderTable({
        as.data.frame(year2varimp)
      },rownames=TRUE)
      output$varimpcox<-renderTable({
        coxvarimp
      },rownames=TRUE,digits=3)
  })





