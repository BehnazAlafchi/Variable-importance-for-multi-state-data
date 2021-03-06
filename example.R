
library(haven)
multistate_data <- as.data.frame(read_sav("J:/multistate_data.sav"))
covs<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9")

res.martingale.resid<-Multistate_residuals(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,dataset=multistate_data,residual="martingale.residual")
res.deviance<-Multistate_residuals(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,dataset=multistate_data,residual="deviance.resisdual")
res.deviance.resid<-Multistate_residuals(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,dataset=multistate_data,residual="deviance")

GBM.varimp.martingale.resid<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.martingale.resid,covs,multistate_data,
                                           method="GBM")

GBM.varimp.martingale.resid<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.martingale.resid,covs,multistate_data,
                                           method="GBM", n.trees=200, shrinkage = 0.5, cv.folds = 5)

RF.varimp.martingale.resid<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.martingale.resid,covs,multistate_data,
                                          method="RF", ntree=100, mtry = 1)

ANN.varimp.martingale.resid<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.martingale.resid,covs,multistate_data,
                                           method="ANN", size=10)

GBM.varimp.deviance<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.deviance,covs,multistate_data,
                                   method="GBM")

RF.varimp.deviance<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.deviance,covs,multistate_data,
                                  method="RF", ntree=1000)

ANN.varimp.deviance<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.deviance,covs,multistate_data,
                                   method="ANN", size=10)

GBM.varimp.deviance.resid<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.deviance.resid,covs,multistate_data,
                                         method="GBM")

RF.varimp.deviance.resid<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.deviance.resid,covs,multistate_data,
                                        method="RF", ntree=1000)

ANN.varimp.deviance.resid<-varImportance(time2=t.trans2,time3=t.trans3,status2=s.trans2,status3=s.trans3,residual=res.deviance.resid,covs,multistate_data,
                                         method="ANN", size=10)
