# provide residuals for an irreversible illness-death model(state1, state2, state3)
# install.packages(c('survival', 'mstate', 'stats'))

Multistate_residuals<-function(time2,time3,status2,status3,dataset,residual
                               # time2: time of transition to state2
                               # time3: time of transition to state3
                               # status2: status for transition to state2 (if the transition is observed equle 1, otherwise 0)
                               # status3: status for transition to state3 (if the transition is observed equle 1, otherwise 0)
                               # residual=c("martingale.residual","deviance.resisdual","deviance") 
                               # the structure of initial dataset should be the same as for mstate package
){
  res<-NULL
  library(survival)
  library(mstate)
  library(stats)
  tmat <- matrix(NA, 3, 3)
  tmat[1, 2:3] <- 1:2
  tmat[2, 3] <- 3
  
  time2<-deparse(substitute(time2)) 
  time3<-deparse(substitute(time3)) 
  status2<-deparse(substitute(status2)) 
  status3<-deparse(substitute(status3)) 
  resid<-c("martingale.residual","deviance.resisdual","deviance")
  
  dimnames(tmat) <- list(from =c("state1","state2","state3") , to =c("state1","state2","state3"))
  msdata <- msprep(time = c(NA, time2, time3), status = c(NA, status2, status3), data = dataset, trans = tmat)
  
  datatrans12<-msdata[ which(msdata$trans==1 ),]
  datatrans13<-msdata[ which(msdata$trans==2 ),]
  datatrans23<-msdata[ which(msdata$trans==3 ),]
  
  # fit null model
  fittrans12<- coxph(Surv(time,status)~1, data=datatrans12, method = "breslow")
  fittrans13<- coxph(Surv(time,status)~1, data=datatrans13, method = "breslow")
  fittrans23<- coxph(Surv(time,status)~1, data=datatrans23, method = "breslow")
  
  #calculate residuals
  if(residual==resid[1]){
    Martingale.residual12<-resid(fittrans12, type=c("martingale"))
    Martingale.residual13<-resid(fittrans13, type=c("martingale"))
    Martingale.residual23<-resid(fittrans23, type=c("martingale"))
    res<-list(Martingale.residual12=Martingale.residual12,Martingale.residual13=Martingale.residual13,Martingale.residual23=Martingale.residual23)
  }
  
  if(residual==resid[2]){
    Deviance.residual12<-resid(fittrans12, type=c("deviance"))
    Deviance.residual13<-resid(fittrans13, type=c("deviance"))
    Deviance.residual23<-resid(fittrans23, type=c("deviance"))
    res<-list(Deviance.residual12=Deviance.residual12,Deviance.residual13=Deviance.residual13,Deviance.residual23=Deviance.residual23)
  }
  
  if(residual==resid[3]){
    Martingale.residual12<-resid(fittrans12, type=c("martingale"))
    Martingale.residual13<-resid(fittrans13, type=c("martingale"))
    Martingale.residual23<-resid(fittrans23, type=c("martingale"))
    Deviance.residual12<-resid(fittrans12, type=c("deviance"))
    Deviance.residual13<-resid(fittrans13, type=c("deviance"))
    Deviance.residual23<-resid(fittrans23, type=c("deviance"))
    Deviance12<-(sign(Martingale.residual12)*(Deviance.residual12))^2
    Deviance13<-(sign(Martingale.residual13)*(Deviance.residual13))^2
    Deviance23<-(sign(Martingale.residual23)*(Deviance.residual23))^2
    res<-list(Deviance12=Deviance12,Deviance13=Deviance13,Deviance23=Deviance23)
  }
  
  return(res)
}