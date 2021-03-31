
packages<-c('dplyr','mice','vtreat','Amelia','caret','doParallel','foreach',
          'ggplot2','rms','parallel','pec','matrixStats','prodlim','qs','ranger','survival','timeROC')
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}
matrixStats::rowSds()
lapply(packages, require, character.only = TRUE)

####Prepare and clean dataset
survdata
####Impute for missing data, 10 imputations, 10 iterations
mice<- mice(survdata,m=10,maxit=10,cluster.seed=500)
qsave(mice,'mice.q')


####Train model on full dataset and tune hyperparameters in tgrid
{
  mice<-qread('mice.q')
  tgrid<-expand.grid(
    num.trees=c(200,250,300),
    .mtry=5:10,
    .splitrule=c("logrank"),
    .min.node.size=c(20,25,30,35,40)
  )
  
  oob<-1
  for (i in 1:length(mice)){
    ###For each imputed dataset train random forest model
  data<-mice::complete(mice,i)
  cl <- makeForkCluster(3)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(ranger),library(rms)))
  system.time(list<-foreach(j=1:nrow(tgrid)) %dopar%{
    rsf<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=tgrid$num.trees[j],mtry=tgrid$.mtry[j],min.node.size=tgrid$.min.node.size[j],splitrule=tgrid$.splitrule)
    ###Calculate out of bag error (equivalent to 1-cindex in out of bag samples as measure of performance) for each combination of tuning parameters
    oob[j]<-rsf$prediction.error
    rm(rsf)
    return(list(oob))
  })
  stopCluster(cl)
  registerDoSEQ()
  for (i in 1:nrow(tgrid)){
    oob[i]<-(list[[i]][[1]][[i]])
  }
  ####select combination of hyperparameters that gives the lowest prediction error and train final model
  rsf<-ranger(Surv(monthssurv,ons_death)~.,data,num.trees=tgrid[which.min(oob),]$num.trees,mtry=tgrid[which.min(oob),]$.mtry,min.node.size=tgrid[which.min(oob),]$.min.node.size,splitrule='logrank',importance='permutation')
  rsfintlist[i]<-list(rsf)
  }

  qsave(rsfintlist,'rsfintlist.q')
}

####Variable selection by bootstrap
{
  mice<-qread('survivalmicetrain.q')
  rsflist<-qread('rsfintlist.q')
  
  ###Calculate Raw VIMP for 1st imputation dataset
  {
  folds<-1:nrow(mice::complete(mice,1))
  finalvimp<-1
  cl <- makeForkCluster(3)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(ranger),library(rms),library(vtreat),library(plotrix),library(mice),library(Amelia)))
  system.time(list<-foreach(t=1:1)%dopar%{
    vimp<-as.data.frame(1:(ncol(mice::complete(mice,1))-2))
    for (n in 1:mice$m){
      datax<-mice::complete(mice,n)
      datax$monthssurv<-mice::complete(mice,n)$monthssurv
      datax$monthssurv<-ceiling(datax$monthssurv)
      data<-datax
      rsf1<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=rsflist[[n]]$num.trees,mtry=rsflist[[n]]$mtry,min.node.size=rsflist[[n]]$min.node.size,splitrule='logrank',importance='permutation')
      vimp[,n]<-rsf1$variable.importance
    }
    
    finalvimp[t]<-list(vimp)
    return(list(finalvimp))
  })
  stopCluster(cl)
  registerDoSEQ()
  
  vimp<-as.data.frame(list[[1]])
  rownames(vimp)<-colnames(mice::complete(mice,1)[,1:(ncol(mice::complete(mice,1))-2)])
  qsave(vimp,'mivimp.q')
}
  ####Create bootstrap resampling 
  {
  data<-mice::complete(mice,1)
  folds<-caret::createResample(1:nrow(mice::complete(mice,1)),times=1000)
  qsave(folds,'folds.q')
  folds<-qread('folds.q')
}
  
  ####Calculate boostrap vimp. Uses hyperparameters from full sample for computational reasons
  {
    finalvimp<-1
    cl <- makeForkCluster(3)
    registerDoParallel(cl)
    system.time(list<-foreach(t=1:1000)%dopar%{
      vimp<-as.data.frame(1:(ncol(mice::complete(mice,1))-2))
      for (n in 1:mice$m){
        
        datax<-mice::complete(mice,n)
        datax$monthssurv<-mice::complete(mice,n)$monthssurv
        datax$monthssurv<-ceiling(datax$monthssurv)
        data<-datax[folds[[t]],]
        
        
        rsf1<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=rsflist[[n]]$num.trees,mtry=rsflist[[n]]$mtry,min.node.size=rsflist[[n]]$min.node.size,splitrule='logrank',importance='permutation')
        vimp[,n]<-rsf1$variable.importance
        
      }
      
      finalvimp[t]<-list(vimp)
      return(list(finalvimp))
    })
    stopCluster(cl)
    registerDoSEQ()
    df<-as.data.frame(1:(ncol(mice::complete(mice,1))-2))
    rownames(df)<-colnames(mice::complete(mice,1)[,1:(ncol(mice::complete(mice,1))-2)])
    
    vimpdf<-array(unlist(df),dim=c((ncol(mice::complete(mice,1))-2),mice$m,length(folds)))
    list
    for (i in 1:1000){
      for (j in 1:mice$m){
        
        vimpdf[,j,i]<-(list[[i]][[1]][[i]][,j])
      }  
    }
    vimpdf
 
  qsave(vimpdf,'vimpdf.q')
  }
  
  
  ###Calculate standard error of vimp from bootstrap samples
  {
  df<-vimpdf
  stderr<-as.data.frame(1:(ncol(mice::complete(mice,1))-2))

  
  
  for (j in 1:mice$m){
    for (k in 1:(ncol(mice::complete(mice,1))-2)){
      stderr[k,j]<-rowSds(as.matrix(df[k,j,]))
    }
  }
  rownames(stderr)<-colnames(mice::complete(mice,1)[,c(1:40,43)])
  stderr
  }

  ###Combine across imputation datasets using Ruben's rules
  {
  mivimp<-qread('mivimp.q')

  rownames(mivimp)<-rownames(stderr)
  vimp<-mi.meld(q=mivimp,se=stderr,byrow=FALSE)
  vimpa<-as.data.frame(t(vimp$q.mi))
  vimpse<-as.data.frame(t(vimp$se.mi))
  vimpf<-cbind2(vimpa,vimpse)
  colnames(vimpf)<-c('VIMP','SE')
  }
  
  ####Select only variables with 99% LCI of >0
  
  vimpf$UCI<-(vimpf$VIMP+(2.576*vimpf$SE))
  vimpf$LCI<-(vimpf$VIMP-(2.576*vimpf$SE))
  vimpf2<-vimpf[order(-vimpf$VIMP),,drop=FALSE]
  vimpf3<-subset(vimpf2,vimpf2$LCI>0)
  vimpf3
  finalvar<-rownames(vimpf3)
  qsave(finalvar,'finalvar.q')

}

####create dummy variables
{
  mice<-qread('survivalmicetrain.q')
  finalvar<-qread('finalvar.q')
  ####Create dummy coding scheme with mortality indicator specified
  OSct<-mkCrossFrameNExperiment(mice::complete(mice,1),finalvar,'ons_death')
  dummies<-OSct$treatments
  qsave(dummies,'dummiesOS.q')
  rm(OSct,dummies)
}

###Train full Random Forest models
{
  mice<-qread('survivalmicetrain.q')
  finalvar<-qread('finalvar.q')
  dummies<-qread('dummiesOS.q')
  tgrid<-expand.grid(
    num.trees=c(200,300,400),
    .mtry=10:20,
    .splitrule=c("logrank"),
    .min.node.size=c(10,20,30,40)
  )
####For each mice imputation, train model, train hyperparameters and put into list
for (i in 1:length(mice)){
  data<-vtreat::prepare(dummies,mice::complete(mice,i),pruneSig=c())
  data$monthssurv<-mice::complete(mice,i)$monthssurv
  data$monthssurv<-ceiling(data$monthssurv)
  oob<-1
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(ranger),library(rms)))
  system.time(list<-foreach(j=1:nrow(tgrid)) %dopar%{
    rsf<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=tgrid$num.trees[j],mtry=tgrid$.mtry[j],min.node.size=tgrid$.min.node.size[j],splitrule=tgrid$.splitrule)
    oob[j]<-rsf$prediction.error
    rm(rsf)
    return(list(oob))
  })
  stopCluster(cl)
  registerDoSEQ()
  for (i in 1:nrow(tgrid)){
    oob[i]<-(list[[i]][[1]][[i]])
  }
  oob1<-oob
  rsf<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=tgrid[which.min(oob),]$num.trees,mtry=tgrid[which.min(oob),]$.mtry,min.node.size=tgrid[which.min(oob),]$.min.node.size,splitrule='logrank',importance='none')
  rsflist[i]<-list(rsf)
}
  qsave(rsflist,'rsflist.q')  
}

####Variable importance in final models by bootstrap (similar to original VIMP calculation)
{
  mice<-qread('survivalmicetrain.q')
  rsflist<-qread('rsflist.q')
  finalvar<-qread('finalvaros.q')
  finalvimp<-1
  cl <- makePSOCKcluster(8)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(ranger),library(rms),library(vtreat),library(plotrix),library(mice),library(Amelia),library(dplyr)))
  n<-1
  system.time(list<-foreach(t=1:1)%dopar%{
    vimp<-as.data.frame(1:(ncol(dplyr::select(mice::complete(mice,n),finalvar))))
    
    for (n in 1:mice$m){
      
      datax<-dplyr::select(mice::complete(mice,n),finalvar,ons_death)
      datax$monthssurv<-mice::complete(mice,n)$monthssurv
      datax$monthssurv<-ceiling(datax$monthssurv)
      data<-datax
      rsf1<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=rsflist[[n]]$num.trees,mtry=ifelse(rsflist[[n]]$mtry<15,rsflist[[n]]$mtry,14),min.node.size=rsflist[[n]]$min.node.size,splitrule='logrank',importance='permutation')
      vimp[,n]<-rsf1$variable.importance
      
    }
    
    finalvimp[t]<-list(vimp)
    return(list(finalvimp))
  })
  stopCluster(cl)
  registerDoSEQ()
  
  vimp<-as.data.frame(list[[1]])
  rownames(vimp)<-finalvar
  vimp
  list
  qsave(vimp,'mivimp2.q')
  
  
  data<-mice::complete(mice,1)
  folds<-qread('folds.q')
  n<-1
  t<-1
  
  {
    
    finalvimp<-1
    cl <- makePSOCKcluster(8)
    registerDoParallel(cl)
    clusterEvalQ(cl,c(library(ranger),library(rms),library(vtreat),library(plotrix),library(mice),library(Amelia),library(dplyr)))
    
    system.time(list<-foreach(t=1:1000)%dopar%{
      vimp<-as.data.frame(1:(ncol(dplyr::select(mice::complete(mice,n),finalvar))))
      
      for (n in 1:mice$m){
        
        datax<-dplyr::select(mice::complete(mice,n),finalvar,ons_death)
        datax$monthssurv<-mice::complete(mice,n)$monthssurv
        datax$monthssurv<-ceiling(datax$monthssurv)
        data<-datax[folds[[t]],]
        
        
        rsf1<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=rsflist[[n]]$num.trees,mtry=ifelse(rsflist[[n]]$mtry<15,rsflist[[n]]$mtry,14),min.node.size=rsflist[[n]]$min.node.size,splitrule='logrank',importance='permutation')
        vimp[,n]<-rsf1$variable.importance
        
      }
      
      finalvimp[t]<-list(vimp)
      return(list(finalvimp))
    })
    stopCluster(cl)
    registerDoSEQ()
    i<-100
    list[[i]][[1]][[i]][,j]
    
    df<-as.data.frame(1:14)
    rownames(df)<-finalvar

    dfa<-array(unlist(df),dim=c(14,mice$m,200))
    for (i in 1:200){
      for (j in 1:mice$m){
        dfa[,j,i]<-(list[[i]][[1]][[i]][,j])
      }  
    }
    dfa
  }
  stderr<-as.data.frame(1:14)
  
  
  for (j in 1:mice$m){
    for (k in 1:14){
      stderr[k,j]<-rowSds(as.matrix(dfa[k,j,]))
    }
  }
  rownames(stderr)<-finalvar
  stderr
  
  mivimp<-qread('mivimp2.q')
  
  vimp<-mi.meld(q=mivimp,se=stderr,byrow=FALSE)
  vimpa<-as.data.frame(t(vimp$q.mi))
  vimpse<-as.data.frame(t(vimp$se.mi))
  vimpf<-cbind2(vimpa,vimpse)
  colnames(vimpf)<-c('VIMP','SE')
  vimpf$UCI<-(vimpf$VIMP+(1.96*vimpf$SE))
  vimpf$LCI<-(vimpf$VIMP-(1.96*vimpf$SE))
  vimpf2<-vimpf[order(-vimpf$VIMP),,drop=FALSE]
  vimpf2

  qsave(vimpf2,'fullvimpfin.q')
  
}

####Calculate censoring times
{
rsflist<-qread('rsflist.q')
mice<-qread('survivalmicetrain.q')
finalvar<-qread('finalvarOS.q')
dummies<-qread('dummiesOS.q')

data<-vtreat::prepare(dummies,mice::complete(mice,1),pruneSig=c())
data$monthssurv<-mice::complete(mice,1)$monthssurv
rsf1<-rsflist[1]
deathtimes<-rsflist[[1]][[3]]
censortimes<-as.data.frame(1)
system.time(for (j in 1:length(deathtimes)){
  for (i in 1:nrow(data)){
    ifelse((data[i,"ons_death"]==1 & data[i,'monthssurv']<=deathtimes[j]), censortimes[i,j]<-1 , (ifelse(data[i,'monthssurv']<=deathtimes[j],censortimes[i,j]<-NA,censortimes[i,j]<-0)))
  }
})

qsave(censortimes,'censortimesfin.q')
}

####Internal validation - Generate bootstrap predictions
{
  rsflist<-qread('rsflist.q')
  mice<-qread('survivalmicetrain.q')
  finalvar<-qread('finalvarOS.q')
  dummies<-qread('dummiesOS.q')
  censortimes<-qread('censortimesfin.q')
  id<-as.data.frame(1:nrow(mice::complete(mice,1)))
  folds<-caret::createResample(1:nrow(mice::complete(mice,1)),times=1000)

  
  cl <- makePSOCKcluster(3)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(ranger),library(rms),library(vtreat),library(plotrix),library(mice),library(Amelia),library(CORElearn),library(miceadds),library(splines)))
  system.time(list<-foreach(t=1:1000)%dopar%{
    estimates<-1
    stander<-1
    estimates2<-1
    stander2<-1
    for (n in 1:length(mice)){
      datax<-vtreat::prepare(dummies,mice::complete(mice,n),pruneSig=c())
      datax$monthssurv<-mice::complete(mice,n)$monthssurv
      datax$monthssurv<-ceiling(datax$monthssurv)
      data<-datax[folds[[t]],]
      rsf1<-ranger(Surv(data$monthssurv,data$ons_death)~.,data,num.trees=rsflist[[n]]$num.trees,mtry=rsflist[[n]]$mtry,min.node.size=rsflist[[n]]$min.node.size,splitrule='logrank',importance='none')
      deathtimes<-predict(rsf1,datax[-folds[[t]],],type='response')$unique.death.times
      
      ###Generate log predictions and standard errors on out of training samples
      estimates[n]<-list(log(predict(rsf1,datax[-folds[[t]],],type='response')$survival))
      survfull<-log(predict(rsf1,datax[-folds[[t]],],type='response',predict.all=TRUE)$survival)
      standerx<-as.data.frame(1:nrow(datax[-folds[[t]],]))
      for (p in 1:length(deathtimes)){
        for (q in 1:nrow(datax[-folds[[t]],])){
          standerx[q,p]<-rowSds(as.matrix(survfull[q,p,c(1:rsf1$num.trees)]))
        }  
      }
      stander[n]<-list(standerx)
      deathtimes2<-predict(rsf1,datax[folds[[t]],],type='response')$unique.death.times
      
      ###Generate log predictions and standard errors on training samples
      estimates2[n]<-list(log(predict(rsf1,datax[folds[[t]],],type='response')$survival))
      survfull2<-log(predict(rsf1,datax[folds[[t]],],type='response',predict.all=TRUE)$survival)
      standerx2<-as.data.frame(1:nrow(datax[-folds[[t]],]))
      for (p in 1:length(deathtimes2)){
        for (q in 1:nrow(datax[folds[[t]],])){
          standerx2[q,p]<-rowSds(as.matrix(survfull2[q,p,c(1:rsf1$num.trees)]))
        }  
      }
      stander2[n]<-list(standerx2)
      
    }
    ####Combine predictions from imputed datasets for final out of training and training predictions
    mipreds<-as.data.frame(1:nrow(datax[-folds[[t]],]))
    for (m in 1:(length(deathtimes)-2)){
      mipreds[,m]<-as.data.frame(t(as.data.frame(mi.meld(q=cbind(estimates[[1]][,m],estimates[[2]][,m],estimates[[3]][,m],estimates[[5]][,m],
                                                                 estimates[[6]][,m],estimates[[7]][,m],estimates[[8]][,m],estimates[[9]][,m],estimates[[10]][,m]),
                                                         se=cbind(stander[[1]][,m],stander[[2]][,m],stander[[3]][,m],stander[[5]][,m],
                                                                  stander[[6]][,m],stander[[7]][,m],stander[[8]][,m],stander[[9]][,m],stander[[10]][,m]),byrow=FALSE)[1])))
    }
    

    mipredsa<-as.data.frame(1:nrow(datax[folds[[t]],]))
    for (m in 1:(length(deathtimes2)-2)){
      mipredsa[,m]<-as.data.frame(t(as.data.frame(mi.meld(q=cbind(estimates2[[1]][,m],estimates2[[2]][,m],estimates2[[3]][,m],estimates2[[5]][,m],
                                                                  estimates2[[6]][,m],estimates2[[7]][,m],estimates2[[8]][,m],estimates2[[9]][,m],estimates2[[10]][,m]),
                                                          se=cbind(stander2[[1]][,m],stander2[[2]][,m],stander2[[3]][,m],stander2[[5]][,m],
                                                                   stander2[[6]][,m],stander2[[7]][,m],stander2[[8]][,m],stander2[[9]][,m],stander2[[10]][,m]),byrow=FALSE)[1])))
    }
    testingpreds<-exp(mipreds)
    rm(mipreds)
    trainingpreds<-exp(mipredsa)
    rm(mipredsa)


    testingpreds$id<-id[-folds[[t]],]
    trainingpreds$id<-id[folds[[t]],]
    output<-list(testingpreds,trainingpreds)
    return(list(output)) 
  })
  stopCluster(cl)
  registerDoSEQ()
  qsave(list,'finallist.q')

}
  
####Internal validation - evaluation metrics
{
  ###Create id/survival dataframe
  ids<-select(survdata,monthssurv,onsdeath)
  ids$id<-1:nrow(ids)
  simplebootstrap<-1
  boot.632<-1
  calibplot<-1
  calibplot632<-1
  fullbt<-1
  full632<-1
  
  ###Separate out dataframes of predictions
  {
    df<-1
    for (i in 1:length(finlist)){
      df[i]<-list(finlist[[i]][[1]][[1]])
    }  
    for (i in 1:length(df)){
      df[i]<-list(merge(df[[i]],ids,by='id'))
    }
    df2<-1
    for (i in 1:length(finlist)){
      df2[i]<-list(dplyr::distinct(finlist[[i]][[1]][[2]]))
    }  
    for (i in 1:length(df2)){
      df2[i]<-list(merge(df2[[i]],ids,by='id'))
    }
    
    ###Recreate dataframe of predictions on original dataset for each bootstrap resample
    {
      df3<-1
      for (i in 1:length(df)){
        x<-df[[i]][,-2]
        y<-df2[[i]][,-2]
        x1<-rbind(x,y)
        df3[i]<-list(dplyr::distinct(x1))
      }
    }
    ###Simple Bootstrap validation
    {
      ###Calculate tROC, c-index and ibrier for bootstrap sample
      troc3<-1
      cid<-1
      ibrier<-1
      for (i in 1:length(df3)){
        ###set to 59 as all patients censored at 60
        troc3[i]<-timeROC(df3[[i]]$monthssurv,df3[[i]]$ons_death,1-df3[[i]][,59],times=59,cause=1,weighting='marginal',iid=FALSE)$AUC[2]
        surv.obj<-with(df3[[i]],Surv(df3[[i]]$monthssurv,df3[[i]]$ons_death))
        cid[i]<-rcorr.cens(x=df3[[i]][,59],S=surv.obj)[[1]]
        dat<-na.omit(df3[[i]])
        ####may need to edit column numbers in dat[,c(2:62)]
        ibrier[i]<-crps(pec(list(calib=as.matrix(dat[,c(2:62)])),Hist(monthssurv,ons_death)~1,data=dat))[2]
      }
      
      ###Generate calibration chart for each bootstrap resample
      plots<-1
      for (k in 1:length(df3)){
        {
          calibrsf<-df3[[k]]
          calib<-calibrsf
          calib$id<-1:nrow(calib)
          censordata<-as.data.frame(df3[[k]]$ons_death)
          censordata$monthssurv<-df3[[k]]$monthssurv
          colnames(censordata)<-c('ons_death','monthssurv')
          for (i in 1:nrow(censordata)){
            ifelse(censordata$monthssurv[i]>60, censordata$ons_death[i]<-0,censordata$ons_death[i]<-censordata$ons_death[i])
          }
          for (i in 1:nrow(censordata)){
            ifelse(censordata$monthssurv[i]>60, censordata$monthssurv[i]<-60,censordata$monthssurv[i]<-censordata$monthssurv[i])
          }
          
          
          calib$death<-censordata$ons_death
          calib$months<-censordata$monthssurv
          colnames(calib)
          
          calib$decile<-with(calibrsf,cut(calibrsf[,59],
                                          breaks = quantile(calibrsf[,59],probs=seq(0,1,by=0.2)),
                                          include.lowest = TRUE))
          
          levels(calib$decile)<-c('0-20','20-40','40-60','60-80','80-100')
          ###calib$decile<-with(calib,cut(V1.50,breaks = quantile(V1.50,probs=seq(0,1,by=0.25)),include.lowest = TRUE))
          ###levels(calib$decile)<-c('0-10','10-20','20-30','30-40')
          mts<-as.data.frame(cbind(censordata$monthssurv,censordata$ons_death))
          names(mts)<-c('monthssurv','ons_death')
          
          predscalib2<-calibrsf
          calib1<-calibrsf
          calib1$ons_death<-censordata$ons_death
          calib1$monthssurv<-censordata$monthssurv
          
          
          {
            predscalib2$decile<-calib$decile
            predscalib2$death<-calib1$ons_death
            predscalib2$months<-calib1$monthssurv
            dec1<-subset(predscalib2,predscalib2$decile=='0-20')
            dec1a<-select(dec1,-decile)
            estimatesdec1<-as.data.frame(1:63)
            estimatesdec1$Time<-estimatesdec1[,1]
            estimatesdec1$survival<-colMeans(dec1a)
            OSKMa <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)   
            times<-as.data.frame(1:60)
            colnames(times)<-'time'
            OSKM2a<-as.data.frame(cbind(OSKMa$time,OSKMa$surv))
            names(OSKM2a)<-c('time','surv') 
            if (nrow(OSKM2a)!=60) OSKM2a<-merge(times,OSKM2a,by='time',all=TRUE)
            if (is.na(OSKM2a[1,2])) OSKM2a[1,2]<-1
            for (i in 1:60){
              if (is.na(OSKM2a[i,2]))
                OSKM2a[i,2]<-OSKM2a[i-1,2]
            }
            predscalib2$decile<-calib$decile
            predscalib2$death<-calib1$ons_death
            predscalib2$months<-calib1$monthssurv
            dec1<-subset(predscalib2,predscalib2$decile=='20-40')
            dec1a<-select(dec1,-decile)
            estimatesdec1b<-as.data.frame(1:63)
            estimatesdec1b$Time<-estimatesdec1b[,1]
            estimatesdec1b$survival<-colMeans(dec1a)
            OSKMb <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
            times<-as.data.frame(1:60)
            colnames(times)<-'time'
            OSKM2b<-as.data.frame(cbind(OSKMb$time,OSKMb$surv))
            names(OSKM2b)<-c('time','surv') 
            if (nrow(OSKM2b)!=60) OSKM2b<-merge(times,OSKM2b,by='time',all=TRUE)
            if (is.na(OSKM2b[1,2])) OSKM2b[1,2]<-1
            for (i in 1:60){
              if (is.na(OSKM2b[i,2]))
                OSKM2b[i,2]<-OSKM2b[i-1,2]
            }
            predscalib2$decile<-calib$decile
            predscalib2$death<-calib1$ons_death
            predscalib2$months<-calib1$monthssurv
            dec1<-subset(predscalib2,predscalib2$decile=='40-60')
            dec1a<-select(dec1,-decile)
            estimatesdec1c<-as.data.frame(1:63)
            estimatesdec1c$Time<-estimatesdec1c[,1]
            estimatesdec1c$survival<-colMeans(dec1a)
            OSKMc <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
            times<-as.data.frame(1:60)
            colnames(times)<-'time'
            OSKM2c<-as.data.frame(cbind(OSKMc$time,OSKMc$surv))
            names(OSKM2c)<-c('time','surv') 
            if (nrow(OSKM2c)!=60) OSKM2c<-merge(times,OSKM2c,by='time',all=TRUE)
            if (is.na(OSKM2c[1,2])) OSKM2c[1,2]<-1
            for (i in 1:60){
              if (is.na(OSKM2c[i,2]))
                OSKM2c[i,2]<-OSKM2c[i-1,2]
            }
            predscalib2$decile<-calib$decile
            predscalib2$death<-calib1$ons_death
            predscalib2$months<-calib1$monthssurv
            dec1<-subset(predscalib2,predscalib2$decile=='60-80')
            dec1a<-select(dec1,-decile)
            estimatesdec1d<-as.data.frame(1:63)
            estimatesdec1d$Time<-estimatesdec1d[,1]
            estimatesdec1d$survival<-colMeans(dec1a)
            OSKMd <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
            times<-as.data.frame(1:60)
            colnames(times)<-'time'
            OSKM2d<-as.data.frame(cbind(OSKMd$time,OSKMd$surv))
            names(OSKM2d)<-c('time','surv') 
            if (nrow(OSKM2d)!=60) OSKM2d<-merge(times,OSKM2d,by='time',all=TRUE)
            if (is.na(OSKM2d[1,2])) OSKM2d[1,2]<-1
            for (i in 1:60){
              if (is.na(OSKM2d[i,2]))
                OSKM2d[i,2]<-OSKM2d[i-1,2]
            }
            predscalib2$decile<-calib$decile
            predscalib2$death<-calib1$ons_death
            predscalib2$months<-calib1$monthssurv
            dec1<-subset(predscalib2,predscalib2$decile=='80-100')
            dec1a<-select(dec1,-decile)
            estimatesdec1e<-as.data.frame(1:63)
            estimatesdec1e$Time<-estimatesdec1e[,1]
            estimatesdec1e$survival<-colMeans(dec1a)
            OSKMe <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
            times<-as.data.frame(1:60)
            colnames(times)<-'time'
            OSKM2e<-as.data.frame(cbind(OSKMe$time,OSKMe$surv))
            names(OSKM2e)<-c('time','surv') 
            if (nrow(OSKM2e)!=60) OSKM2e<-merge(times,OSKM2e,by='time',all=TRUE)
            if (is.na(OSKM2e[1,2])) OSKM2e[1,2]<-1
            for (i in 1:60){
              if (is.na(OSKM2e[i,2]))
                OSKM2e[i,2]<-OSKM2e[i-1,2]
            }
            
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
            
            full
            
            plot<-ggplot(full,aes(Time,survival))+
              geom_line(data=full,mapping=aes(colour=col,linetype=lt))+
              xlim(0,60)+ylim(0,1)+
              scale_color_discrete(name='')+
              scale_linetype_discrete(name='Probability Quintile',
                                      labels=c('0-20','20-40','40-60','60-80','80-100'))+
              theme_bw()
            
          }
        }
        plots[k]<-list(plot)
      }
      
      ###Average calibration chart for final simple bootstrap calibration chart
      surv<-as.data.frame(1:615)
      for (i in 1:length(plots)){
        surv[,i]<-plots[[i]]$data$survival
      }
      survs<-rowMeans(surv)
      
      full$survival<-survs
      
      fullx<-full[full$Time==1,]
      fullx$Time<-0
      fullx$survival<-1
      full2<-rbind(full,fullx)
      plot2<-ggplot(full2,aes(Time,survival))+
        geom_line(data=full2,mapping=aes(colour=col,linetype=lt))+
        xlim(0,60)+ylim(0,1)+
        scale_color_discrete(name='')+
        scale_linetype_discrete(name='Probability Quintile',
                                labels=c('0-20','20-40','40-60','60-80','80-100'))+
        theme_bw()
      bootstrapcalibration<-plot2
      
      
      ###Combine validation metrics across bootstrap resamples
      out<-as.data.frame(t(as.data.frame(c(mean(troc3),quantile(troc3,probs=c(0.025,0.975))))))
      out<-rbind(out,c(mean(cid),quantile(cid,probs=c(0.025,0.975))))
      out<-rbind(out,c(mean(ibrier),quantile(ibrier,probs=c(0.025,0.975))))
      colnames(out)<-c('mean','2.5%','97.5%')
      rownames(out)<-c('tROC','CiD','iBrier')
      simplebootstrap<-list(out)
      rm(out)
      
    }
    simplebootstrap
    bootstrapcalibration
    
    ###0.632 bootstrap validation
    {  
      ###Calculate tROC, c-index and ibrier for Testing Samples
      {
        trocTE<-1
        cidTE<-1
        ibrierTE<-1
        for (i in 1:length(df)){
          trocTE[i]<-timeROC(df[[i]]$monthssurv,df[[i]]$ons_death,1-df[[i]][,59],times=59,cause=1,weighting='marginal',iid=FALSE)$AUC[2]
          surv.obj<-with(df[[i]],Surv(df[[i]]$monthssurv,df[[i]]$ons_death))
          cidTE[i]<-rcorr.cens(x=df[[i]][,59],S=surv.obj)[[1]]
          dat<-na.omit(df[[i]])
          ibrierTE[i]<-crps(pec(list(calib=as.matrix(dat[,c(2:62)])),Hist(monthssurv,ons_death)~1,data=dat))[2]
        }
      }
      ###Calculate tROC, c-index and ibrier for Training Samples  
      {
        trocTR<-1
        cidTR<-1
        ibrierTR<-1
        for (i in 1:length(df2)){
          trocTR[i]<-timeROC(df2[[i]]$monthssurv,df2[[i]]$ons_death,1-df2[[i]][,59],times=59,cause=1,weighting='marginal',iid=FALSE)$AUC[2]
          surv.obj<-with(df2[[i]],Surv(df2[[i]]$monthssurv,df2[[i]]$ons_death))
          cidTR[i]<-rcorr.cens(x=df2[[i]][,59],S=surv.obj)[[1]]
          dat<-na.omit(df2[[i]])
          ibrierTR[i]<-crps(pec(list(calib=as.matrix(dat[,c(2:62)])),Hist(monthssurv,ons_death)~1,data=dat))[2]
        }
      }
      ###Combine testing and training in 0.632/0.368 ratio
      {
        troc632<-((trocTR*0.368)+(trocTE*0.632))
        cid632<-((cidTR*0.368)+(cidTE*0.632))
        ibrier632<-((ibrierTR*0.368)+(ibrierTE*0.632))
        out1<-as.data.frame(t(as.data.frame(c(mean(troc632),quantile(troc632,probs=c(0.025,0.975))))))
        out1<-rbind(out1,c(mean(cid632),quantile(cid632,probs=c(0.025,0.975))))
        out1<-rbind(out1,c(mean(ibrier632),quantile(ibrier632,probs=c(0.025,0.975))))
        colnames(out1)<-c('mean','2.5%','97.5%')
        rownames(out1)<-c('tROC','CiD','iBrier')
        boot.632<-list(out1)
        rm(out1)
      }
      
      ###Quintile calibration plots
      plots<-1
      {
        for (k in 1:length(df)){
          {
            ###Generate calibration plots for testing cases  
            calibrsf<-df[[k]]
            calibrsf
            calib<-calibrsf
            calib$id<-1:nrow(calib)
            censordata<-as.data.frame(df[[k]]$ons_death)
            censordata$monthssurv<-df[[k]]$monthssurv
            colnames(censordata)<-c('ons_death','monthssurv')
            for (i in 1:nrow(censordata)){
              ifelse(censordata$monthssurv[i]>60, censordata$ons_death[i]<-0,censordata$ons_death[i]<-censordata$ons_death[i])
            }
            for (i in 1:nrow(censordata)){
              ifelse(censordata$monthssurv[i]>60, censordata$monthssurv[i]<-60,censordata$monthssurv[i]<-censordata$monthssurv[i])
            }
            
            
            calib$death<-censordata$ons_death
            calib$months<-censordata$monthssurv
            
            colnames(calib)
            calib$decile<-with(calibrsf,cut(calibrsf[,59],
                                            breaks = quantile(calibrsf[,59],probs=seq(0,1,by=0.2)),
                                            include.lowest = TRUE))
            levels(calib$decile)<-c('0-20','20-40','40-60','60-80','80-100')
            ###calib$decile<-with(calib,cut(V1.50,breaks = quantile(V1.50,probs=seq(0,1,by=0.25)),include.lowest = TRUE))
            ###levels(calib$decile)<-c('0-10','10-20','20-30','30-40')
            mts<-as.data.frame(cbind(censordata$monthssurv,censordata$ons_death))
            names(mts)<-c('monthssurv','ons_death')
            
            predscalib2<-calibrsf
            calib1<-calibrsf
            calib1$ons_death<-censordata$ons_death
            calib1$monthssurv<-censordata$monthssurv
            
            
            {
              predscalib2$decile<-calib$decile
              predscalib2$death<-calib1$ons_death
              predscalib2$months<-calib1$monthssurv
              dec1<-subset(predscalib2,predscalib2$decile=='0-20')
              dec1a<-select(dec1,-decile)
              estimatesdec1<-as.data.frame(1:63)
              estimatesdec1$Time<-estimatesdec1[,1]
              estimatesdec1$survival<-colMeans(dec1a)
              OSKMa <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)   
              times<-as.data.frame(1:60)
              colnames(times)<-'time'
              OSKM2a<-as.data.frame(cbind(OSKMa$time,OSKMa$surv))
              names(OSKM2a)<-c('time','surv') 
              if (nrow(OSKM2a)!=60) OSKM2a<-merge(times,OSKM2a,by='time',all=TRUE)
              if (is.na(OSKM2a[1,2])) OSKM2a[1,2]<-1
              for (i in 1:60){
                if (is.na(OSKM2a[i,2]))
                  OSKM2a[i,2]<-OSKM2a[i-1,2]
              }
              predscalib2$decile<-calib$decile
              predscalib2$death<-calib1$ons_death
              predscalib2$months<-calib1$monthssurv
              dec1<-subset(predscalib2,predscalib2$decile=='20-40')
              dec1a<-select(dec1,-decile)
              estimatesdec1b<-as.data.frame(1:63)
              estimatesdec1b$Time<-estimatesdec1b[,1]
              estimatesdec1b$survival<-colMeans(dec1a)
              OSKMb <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
              times<-as.data.frame(1:60)
              colnames(times)<-'time'
              OSKM2b<-as.data.frame(cbind(OSKMb$time,OSKMb$surv))
              names(OSKM2b)<-c('time','surv') 
              if (nrow(OSKM2b)!=60) OSKM2b<-merge(times,OSKM2b,by='time',all=TRUE)
              if (is.na(OSKM2b[1,2])) OSKM2b[1,2]<-1
              for (i in 1:60){
                if (is.na(OSKM2b[i,2]))
                  OSKM2b[i,2]<-OSKM2b[i-1,2]
              }
              predscalib2$decile<-calib$decile
              predscalib2$death<-calib1$ons_death
              predscalib2$months<-calib1$monthssurv
              dec1<-subset(predscalib2,predscalib2$decile=='40-60')
              dec1a<-select(dec1,-decile)
              estimatesdec1c<-as.data.frame(1:63)
              estimatesdec1c$Time<-estimatesdec1c[,1]
              estimatesdec1c$survival<-colMeans(dec1a)
              OSKMc <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
              times<-as.data.frame(1:60)
              colnames(times)<-'time'
              OSKM2c<-as.data.frame(cbind(OSKMc$time,OSKMc$surv))
              names(OSKM2c)<-c('time','surv') 
              if (nrow(OSKM2c)!=60) OSKM2c<-merge(times,OSKM2c,by='time',all=TRUE)
              if (is.na(OSKM2c[1,2])) OSKM2c[1,2]<-1
              for (i in 1:60){
                if (is.na(OSKM2c[i,2]))
                  OSKM2c[i,2]<-OSKM2c[i-1,2]
              }
              predscalib2$decile<-calib$decile
              predscalib2$death<-calib1$ons_death
              predscalib2$months<-calib1$monthssurv
              dec1<-subset(predscalib2,predscalib2$decile=='60-80')
              dec1a<-select(dec1,-decile)
              estimatesdec1d<-as.data.frame(1:63)
              estimatesdec1d$Time<-estimatesdec1d[,1]
              estimatesdec1d$survival<-colMeans(dec1a)
              OSKMd <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
              times<-as.data.frame(1:60)
              colnames(times)<-'time'
              OSKM2d<-as.data.frame(cbind(OSKMd$time,OSKMd$surv))
              names(OSKM2d)<-c('time','surv') 
              if (nrow(OSKM2d)!=60) OSKM2d<-merge(times,OSKM2d,by='time',all=TRUE)
              if (is.na(OSKM2d[1,2])) OSKM2d[1,2]<-1
              for (i in 1:60){
                if (is.na(OSKM2d[i,2]))
                  OSKM2d[i,2]<-OSKM2d[i-1,2]
              }
              predscalib2$decile<-calib$decile
              predscalib2$death<-calib1$ons_death
              predscalib2$months<-calib1$monthssurv
              dec1<-subset(predscalib2,predscalib2$decile=='80-100')
              dec1a<-select(dec1,-decile)
              estimatesdec1e<-as.data.frame(1:63)
              estimatesdec1e$Time<-estimatesdec1e[,1]
              estimatesdec1e$survival<-colMeans(dec1a)
              OSKMe <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
              times<-as.data.frame(1:60)
              colnames(times)<-'time'
              OSKM2e<-as.data.frame(cbind(OSKMe$time,OSKMe$surv))
              names(OSKM2e)<-c('time','surv') 
              if (nrow(OSKM2e)!=60) OSKM2e<-merge(times,OSKM2e,by='time',all=TRUE)
              if (is.na(OSKM2e[1,2])) OSKM2e[1,2]<-1
              for (i in 1:60){
                if (is.na(OSKM2e[i,2]))
                  OSKM2e[i,2]<-OSKM2e[i-1,2]
              }
              
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
              
              full
              
              plot<-ggplot(full,aes(Time,survival))+
                geom_line(data=full,mapping=aes(colour=col,linetype=lt))+
                xlim(0,60)+ylim(0,1)+
                scale_color_discrete(name='')+
                scale_linetype_discrete(name='Probability Quintile',
                                        labels=c('0-20','20-40','40-60','60-80','80-100'))+
                theme_bw()
              
            }
          }
          ####0.632 Calibration plots in each resample
          {
            {###Generate calibration plots for training cases  
              calibrsf<-df2[[k]]
              calibrsf<-calibrsf[,-c(1,2,81,82)]
              calibrsf<-calibrsf[,1:61]
              calibrsf
              calib<-calibrsf
              calib$id<-1:nrow(calib)
              censordata<-as.data.frame(df2[[k]]$ons_death)
              censordata$monthssurv<-df2[[k]]$monthssurv
              colnames(censordata)<-c('ons_death','monthssurv')
              for (i in 1:nrow(censordata)){
                ifelse(censordata$monthssurv[i]>60, censordata$ons_death[i]<-0,censordata$ons_death[i]<-censordata$ons_death[i])
              }
              for (i in 1:nrow(censordata)){
                ifelse(censordata$monthssurv[i]>60, censordata$monthssurv[i]<-60,censordata$monthssurv[i]<-censordata$monthssurv[i])
              }
              
              
              calib$death<-censordata$ons_death
              calib$months<-censordata$monthssurv
              
              colnames(calib)
              calib$decile<-with(calibrsf,cut(calibrsf[,37],
                                              breaks = quantile(calibrsf[,37],probs=seq(0,1,by=0.2)),
                                              include.lowest = TRUE))
              levels(calib$decile)<-c('0-20','20-40','40-60','60-80','80-100')
              ###calib$decile<-with(calib,cut(V1.50,breaks = quantile(V1.50,probs=seq(0,1,by=0.25)),include.lowest = TRUE))
              ###levels(calib$decile)<-c('0-10','10-20','20-30','30-40')
              mts<-as.data.frame(cbind(censordata$monthssurv,censordata$ons_death))
              names(mts)<-c('monthssurv','ons_death')
              
              predscalib2<-calibrsf
              calib1<-calibrsf
              calib1$ons_death<-censordata$ons_death
              calib1$monthssurv<-censordata$monthssurv
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile=='0-20')
                dec1a<-select(dec1,-decile)
                estimatesdec1<-as.data.frame(1:63)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                OSKMa <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                times<-as.data.frame(1:60)
                colnames(times)<-'time'
                OSKM2a<-as.data.frame(cbind(OSKMa$time,OSKMa$surv))
                names(OSKM2a)<-c('time','surv') 
                if (nrow(OSKM2a)!=60) OSKM2a<-merge(times,OSKM2a,by='time',all=TRUE)
                if (is.na(OSKM2a[1,2])) OSKM2a[1,2]<-1
                for (i in 1:60){
                  if (is.na(OSKM2a[i,2]))
                    OSKM2a[i,2]<-OSKM2a[i-1,2]
                }
                
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile=='20-40')
                dec1a<-select(dec1,-decile)
                estimatesdec1b<-as.data.frame(1:63)
                estimatesdec1b$Time<-estimatesdec1b[,1]
                estimatesdec1b$survival<-colMeans(dec1a)
                OSKMb <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                times<-as.data.frame(1:60)
                colnames(times)<-'time'
                OSKM2b<-as.data.frame(cbind(OSKMb$time,OSKMb$surv))
                names(OSKM2b)<-c('time','surv') 
                if (nrow(OSKM2b)!=60) OSKM2b<-merge(times,OSKM2b,by='time',all=TRUE)
                if (is.na(OSKM2b[1,2])) OSKM2b[1,2]<-1
                for (i in 1:60){
                  if (is.na(OSKM2b[i,2]))
                    OSKM2b[i,2]<-OSKM2b[i-1,2]
                }
                
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile=='40-60')
                dec1a<-select(dec1,-decile)
                estimatesdec1c<-as.data.frame(1:63)
                estimatesdec1c$Time<-estimatesdec1c[,1]
                estimatesdec1c$survival<-colMeans(dec1a)
                OSKMc <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                times<-as.data.frame(1:60)
                colnames(times)<-'time'
                OSKM2c<-as.data.frame(cbind(OSKMc$time,OSKMc$surv))
                names(OSKM2c)<-c('time','surv') 
                if (nrow(OSKM2c)!=60) OSKM2c<-merge(times,OSKM2c,by='time',all=TRUE)
                if (is.na(OSKM2c[1,2])) OSKM2c[1,2]<-1
                for (i in 1:60){
                  if (is.na(OSKM2c[i,2]))
                    OSKM2c[i,2]<-OSKM2c[i-1,2]
                }
                
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile=='60-80')
                dec1a<-select(dec1,-decile)
                estimatesdec1d<-as.data.frame(1:63)
                estimatesdec1d$Time<-estimatesdec1d[,1]
                estimatesdec1d$survival<-colMeans(dec1a)
                OSKMd <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                times<-as.data.frame(1:60)
                colnames(times)<-'time'
                OSKM2d<-as.data.frame(cbind(OSKMd$time,OSKMd$surv))
                names(OSKM2d)<-c('time','surv') 
                if (nrow(OSKM2d)!=60) OSKM2d<-merge(times,OSKM2d,by='time',all=TRUE)
                if (is.na(OSKM2d[1,2])) OSKM2d[1,2]<-1
                for (i in 1:60){
                  if (is.na(OSKM2d[i,2]))
                    OSKM2d[i,2]<-OSKM2d[i-1,2]
                }
                
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile=='80-100')
                dec1a<-select(dec1,-decile)
                estimatesdec1e<-as.data.frame(1:63)
                estimatesdec1e$Time<-estimatesdec1e[,1]
                estimatesdec1e$survival<-colMeans(dec1a)
                OSKMe <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                times<-as.data.frame(1:60)
                colnames(times)<-'time'
                OSKM2e<-as.data.frame(cbind(OSKMe$time,OSKMe$surv))
                names(OSKM2e)<-c('time','surv') 
                if (nrow(OSKM2e)!=60) OSKM2e<-merge(times,OSKM2e,by='time',all=TRUE)
                if (is.na(OSKM2e[1,2])) OSKM2e[1,2]<-1
                for (i in 1:60){
                  if (is.na(OSKM2e[i,2]))
                    OSKM2e[i,2]<-OSKM2e[i-1,2]
                }
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
                plot1<-ggplot(full,aes(Time,survival))+
                  geom_line(data=full,mapping=aes(colour=col,linetype=lt))+
                  xlim(0,60)+ylim(0,1)+
                  scale_color_discrete(name='')+
                  scale_linetype_discrete(name='Probability Quintile',
                                          labels=c('0-20','20-40','40-60','60-80','80-100'))+
                  theme_bw()
                
              }
            }
            ####0.632 combination of testing and training cases
            plot3<-plot
            plot3$data$survival<-(0.632*(plot$data$survival)+0.368*(plot1$data$survival))
            plots[k]<-list(plot3)
          }
        }
        
        ####Average calibration plots across all bootstrap resamples
        {
          surv<-as.data.frame(1:615)
          for (i in 1:length(plots)){
            surv[,i]<-plots[[i]]$data$survival
          }
          survs<-rowMeans(surv)
          
          full$survival<-survs
          
          plot3<-ggplot(full,aes(Time,survival))+
            geom_line(data=full,mapping=aes(colour=col,linetype=lt))+
            xlim(0,60)+ylim(0,1)+
            scale_color_discrete(name='')+
            scale_linetype_discrete(name='Probability Quintile',
                                    labels=c('0-20','20-40','40-60','60-80','80-100'))+
            theme_bw()
          
          calibplot632<-list(plot3)
        }
        
      }
    }
    boot.632
    calibplot632
    
    ###decile plots RSF/oCPH 12s for 10reps
    system.time({
      decdata632<-1
      oskmdata632<-1
      modelcalprob<-1
      modelkmprob<-1
      for (m in c(2,4)){
        
        z<-m*2
        df<-1
        for (i in 1:length(finlist)){
          df[i]<-list(finlist[[i]][[1]][[z]])
        }  
        for (i in 1:length(df)){
          df[i]<-list(merge(df[[i]],ids,by='id'))
        }
        df2<-1
        for (i in 1:length(finlist)){
          df2[i]<-list(dplyr::distinct(finlist[[i]][[1]][[z+1]]))
        }  
        for (i in 1:length(df2)){
          df2[i]<-list(merge(df2[[i]],ids,by='id'))
        }
        for (k in 1:length(finlist)){
          {
            calibrsf<-df[[k]]
            calibrsf<-calibrsf[,-c(1,2,81,82)]
            calibrsf<-calibrsf[,1:61]
            calib<-calibrsf
            calib$id<-1:nrow(calib)
            censordata<-as.data.frame(df[[k]]$ons_death)
            censordata$monthssurv<-df[[k]]$monthssurv
            colnames(censordata)<-c('ons_death','monthssurv')
            for (i in 1:nrow(censordata)){
              ifelse(censordata$monthssurv[i]>60, censordata$ons_death[i]<-0,censordata$ons_death[i]<-censordata$ons_death[i])
            }
            for (i in 1:nrow(censordata)){
              ifelse(censordata$monthssurv[i]>60, censordata$monthssurv[i]<-60,censordata$monthssurv[i]<-censordata$monthssurv[i])
            }
            
            
            calib$death<-censordata$ons_death
            calib$months<-censordata$monthssurv
            
            colnames(calib)
            calib$decile<-with(calibrsf,cut(calibrsf[,37],
                                            breaks = quantile(calibrsf[,37],probs=seq(0,1,by=0.1)),
                                            include.lowest = TRUE))
            levels(calib$decile)<-c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')
            ###calib$decile<-with(calib,cut(V1.50,breaks = quantile(V1.50,probs=seq(0,1,by=0.25)),include.lowest = TRUE))
            ###levels(calib$decile)<-c('0-10','10-20','20-30','30-40')
            mts<-as.data.frame(cbind(censordata$monthssurv,censordata$ons_death))
            names(mts)<-c('monthssurv','ons_death')
            
            predscalib2<-calibrsf
            calib1<-calibrsf
            calib1$ons_death<-censordata$ons_death
            calib1$monthssurv<-censordata$monthssurv
            
            {
              group<-'0-10'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                a<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                aa<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                ab<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'10-20'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                b<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                ba<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                bb<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'20-30'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                c<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                ca<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                cb<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'30-40'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                d<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                da<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                db<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'40-50'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                e<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                ea<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                eb<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'50-60'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                f<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                fa<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                fb<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'60-70'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                g<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                ga<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                gb<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'70-80'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                h<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                ha<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                hb<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'80-90'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                ii<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                ia<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                ib<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'90-100'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                j<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                ja<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                jb<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
            }
          }
          {
            calibrsf<-df2[[k]]
            calibrsf<-calibrsf[,-c(1,2,81,82)]
            calibrsf<-calibrsf[,1:61]
            calib<-calibrsf
            calib$id<-1:nrow(calib)
            censordata<-as.data.frame(df2[[k]]$ons_death)
            censordata$monthssurv<-df2[[k]]$monthssurv
            colnames(censordata)<-c('ons_death','monthssurv')
            for (i in 1:nrow(censordata)){
              ifelse(censordata$monthssurv[i]>60, censordata$ons_death[i]<-0,censordata$ons_death[i]<-censordata$ons_death[i])
            }
            for (i in 1:nrow(censordata)){
              ifelse(censordata$monthssurv[i]>60, censordata$monthssurv[i]<-60,censordata$monthssurv[i]<-censordata$monthssurv[i])
            }
            
            
            calib$death<-censordata$ons_death
            calib$months<-censordata$monthssurv
            
            colnames(calib)
            calib$decile<-with(calibrsf,cut(calibrsf[,37],
                                            breaks = quantile(calibrsf[,37],probs=seq(0,1,by=0.1)),
                                            include.lowest = TRUE))
            levels(calib$decile)<-c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')
            ###calib$decile<-with(calib,cut(V1.50,breaks = quantile(V1.50,probs=seq(0,1,by=0.25)),include.lowest = TRUE))
            ###levels(calib$decile)<-c('0-10','10-20','20-30','30-40')
            mts<-as.data.frame(cbind(censordata$monthssurv,censordata$ons_death))
            names(mts)<-c('monthssurv','ons_death')
            
            predscalib2<-calibrsf
            calib1<-calibrsf
            calib1$ons_death<-censordata$ons_death
            calib1$monthssurv<-censordata$monthssurv
            
            {
              group<-'0-10'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                a2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                a2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                a2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'10-20'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                b2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                b2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                b2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'20-30'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                c2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                c2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                c2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'30-40'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                d2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                d2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                d2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'40-50'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                e2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                e2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                e2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'50-60'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                f2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                f2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                f2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'60-70'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                g2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                g2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                g2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'70-80'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                h2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                h2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                h2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'80-90'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                i2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                i2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                i2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
              group<-'90-100'
              {
                predscalib2$decile<-calib$decile
                predscalib2$death<-calib1$ons_death
                predscalib2$months<-calib1$monthssurv
                dec1<-subset(predscalib2,predscalib2$decile==group)
                dec1a<-select(dec1,-decile,-months,-death)
                dec2<-subset(calib,calib$decile==group)
                dec2a<-select(dec2,-decile,-months,-death,-id)
                estimatesdec1<-as.data.frame(1:61)
                estimatesdec1$Time<-estimatesdec1[,1]
                estimatesdec1$survival<-colMeans(dec1a)
                estimatesdec2<-as.data.frame(1:61)
                estimatesdec2$Time<-estimatesdec2[,1]
                estimatesdec2$survival<-colMeans(dec2a)
                OSKM <- survfit(Surv(dec1$months, dec1$death)~1, data=dec1)                                 
                OSKM2<-as.data.frame(cbind(OSKM$time,OSKM$surv))
                names(OSKM2)<-c('time','surv')  
                j2<-ggplot(estimatesdec1,aes(Time,survival))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)
                j2a<-ggplot(OSKM2,aes(time,surv))+geom_line(colour='red')+xlim(0,60)+ylim(0,1)
                j2b<-ggplot(estimatesdec2,aes(Time,survival))+geom_line(colour='green')+xlim(0,60)+ylim(0,1)
              }
            }
          }
          blackplotte<-list(a,b,c,d,e,f,g,h,ii,j)
          blackplottr<-list(a2,b2,c2,d2,e2,f2,g2,h2,i2,j2)
          redplotte<-list(aa,ba,ca,da,ea,fa,ga,ha,ia,ja)
          redplottr<-list(a2a,b2a,c2a,d2a,e2a,f2a,g2a,h2a,i2a,j2a)
          
          blackplot632<-1
          for (qq in 1:10){
            newdf<-as.data.frame((0.632*approx(blackplotte[[qq]]$data[,-1],n=100)$y)+(0.368*approx(blackplottr[[qq]]$data[,-1],n=100)$y))
            newdf$time<-((0.632*approx(blackplotte[[qq]]$data,n=100)$x)+(0.368*approx(blackplottr[[qq]]$data,n=100)$x))
            colnames(newdf)<-c('surv','time')
            blackplot632[qq]<-list(newdf)
          }
          redplot632<-1
          for (qqq in 1:10){
            newdf<-as.data.frame((0.632*approx(redplotte[[qqq]]$data,n=100)$y)+(0.368*approx(redplottr[[qqq]]$data,n=100)$y))
            newdf$time<-((0.632*approx(redplotte[[qqq]]$data,n=100)$x)+(0.368*approx(redplottr[[qqq]]$data,n=100)$x))
            colnames(newdf)<-c('surv','time')
            redplot632[qqq]<-list(newdf)
          }
          
          
          decdata632[k]<-list(blackplot632)
          oskmdata632[k]<-list(redplot632)
        }
        
        calprobdec<-1
        probdec<-1
        
        for (t in 1:10){
          probdec<-decdata632[[1]][[t]]
          surv<-as.data.frame(1:100)
          for (i in 1:length(finlist)){
            surv[,i]<-decdata632[[i]][[t]]$surv
          }
          probdec$surv<-rowMeans(surv)
          calprobdec[t]<-list(probdec)
        }
        
        kmprobdec<-1
        probdec<-1
        
        for (t in 1:10){
          probdec<-oskmdata632[[1]][[t]]
          surv<-as.data.frame(1:100)
          for (i in 1:length(finlist)){
            surv[,i]<-oskmdata632[[i]][[t]]$surv
          }
          probdec$surv<-rowMeans(surv)
          kmprobdec[t]<-list(probdec)
        }
        modelcalprob[m]<-list(calprobdec)
        modelkmprob[m]<-list(kmprobdec)
      }
      plots<-1
      for (i in 1:10){
        plots[i]<-list(ggplot(modelcalprob[[2]][[i]],aes(time,surv))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)+
                         geom_line(data=modelkmprob[[2]][[i]],aes(time,surv),colour='red')
        )
      }
      
      
      plots2<-1
      for (i in 1:10){
        plots2[i]<-list(ggplot(modelcalprob[[4]][[i]],aes(time,surv))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)+
                          geom_line(data=modelkmprob[[4]][[i]],aes(time,surv),colour='red'))
      }
      
      plots3<-1
      for (i in 1:10){
        plots3[i]<-list(ggplot(modelcalprob[[2]][[i]],aes(time,surv))+geom_line(colour='black')+xlim(0,60)+ylim(0,1)+
                          geom_line(data=modelkmprob[[2]][[i]],aes(time,surv),colour='red')+
                          geom_line(data=modelcalprob[[4]][[i]],aes(time,surv),colour='green'))
      }
      
    })
    qsave(modelcalprob,'modelcalprob.q')
    qsave(modelkmprob,'modelkmprob.q')
    qsave(plots,'deccalplotrsf.q')
    qsave(plots2,'deccalplotocph.q')
    qsave(plots3,'deccalplotorsfcph.q')
    qsave(simplebootstrap,'simplebootstrap.q')
    qsave(boot.632,'boot.632.q')
    qsave(calibplot,'calibplot1.q')
    qsave(calibplot632,'calibplot632.q')
    qsave(fullbt,'fullbt.q')
    qsave(full632,'full632.q')
    qsave(trocplots,'trocplots.q')
    
  }
  
}