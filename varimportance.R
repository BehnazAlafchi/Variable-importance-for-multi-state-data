# require("gbm")
# require("randomForest")
# require("nnet")
# require("NeuralNetTools")


varImportance<-function(time2,time3,status2,status3,residual,covs,dataset,method,
                        #GBM
                        distribution = "gaussian", n.trees=500 , shrinkage=0.1 , cv.folds=1, GBM.weights=NULL,
                        var.monotone = NULL, interaction.depth = 1, train.fraction = 1, keep.data = FALSE,
                        verbose = FALSE, class.stratify.cv = NULL, n.cores = NULL,
                        #ANN
                        size=2, ANN.weights=NULL, contrasts=NULL, linout = FALSE, entropy = FALSE, 
                        softmax = FALSE, censored = FALSE, skip = FALSE, rang = 0.7, decay = 0,
                        maxit = 100, Hess = FALSE, trace = TRUE, MaxNWts = 1000, abstol = 1.0e-4, reltol = 1.0e-8, 
                        #RF
                        ntree=500, mtry=1, replace=TRUE, strata=NULL, nodesize = 5,
                        maxnodes = NULL, importance=TRUE, localImp=FALSE, nPerm=1,
                        #
                        p=1, na.action=NULL, subset=NULL,
                        ...
                        # residual is an output of Multistate_residuals function
                        # p is the probability used to construct train data (defult is 1 (use all data))
                        # covs is the vector of covariate's name
                        # method=c("GBM","RF","ANN") # Gradient Boosting(GBM) & Random Forest (RF) & Neural Network (ANN)
){
  
  varimp<-NULL
  meth<-c("GBM","RF","ANN")
  
  library(survival)
  library(mstate)
  tmat <- matrix(NA, 3, 3)
  tmat[1, 2:3] <- 1:2
  tmat[2, 3] <- 3
  
  time2<-deparse(substitute(time2)) 
  time3<-deparse(substitute(time3)) 
  status2<-deparse(substitute(status2)) 
  status3<-deparse(substitute(status3)) 
  
  dimnames(tmat) <- list(from =c("state1","state2","state3") , to =c("state1","state2","state3"))
  msdata <- msprep(time = c(NA, time2, time3), status = c(NA, status2, status3), data = dataset, trans = tmat,  keep = covs)
  msdata <- expand.covs(msdata, covs, append = TRUE, longnames = FALSE)
  
  datatrans12<-msdata[ which(msdata$trans==1 ),]
  datatrans13<-msdata[ which(msdata$trans==2 ),]
  datatrans23<-msdata[ which(msdata$trans==3 ),]
  
  a=8+length(covs)+1
  b=NULL
  for (i in 1:(length(covs)-1)){
    a=a+3
    b<-c(b,a)
  }
  b
  data12<-datatrans12[c(8+length(covs)+1,b)]
  data13<-datatrans13[c(8+length(covs)+2,b+1)]
  data23<-datatrans23[c(8+length(covs)+3,b+2)]
  
  data_train12 <- sample(1:nrow(data12), size = round(p * nrow(data12)),replace = FALSE)
  data_train13 <- sample(1:nrow(data13), size = round(p * nrow(data13)),replace = FALSE)
  data_train23 <- sample(1:nrow(data23), size = round(p * nrow(data23)),replace = FALSE)
  
  
  # Gradient boosting
  if(method==meth[1]){
    library("gbm")
    b12<-NULL
    b13<-NULL
    b23<-NULL
    boosttrans12 <- gbm(
      residual[[1]]~ ., data =data12,
      n.trees = n.trees,
      shrinkage = shrinkage,
      cv.folds = cv.folds,
      weights=GBM.weights,
      var.monotone = var.monotone,
      interaction.depth = interaction.depth,
      train.fraction = train.fraction,
      keep.data = keep.data,
      verbose = verbose,
      class.stratify.cv = class.stratify.cv,
      n.cores = n.cores)
    
    boosttrans13 <- gbm(
      residual[[2]]~ ., data =data13,
      n.trees = n.trees,
      shrinkage = shrinkage,
      cv.folds = cv.folds,
      weights=GBM.weights,
      var.monotone = var.monotone,
      interaction.depth = interaction.depth,
      train.fraction = train.fraction,
      keep.data = keep.data,
      verbose = verbose,
      class.stratify.cv = class.stratify.cv,
      n.cores = n.cores)
    
    boosttrans23 <- gbm(
      residual[[3]]~ ., data =data23,
      n.trees = n.trees,
      shrinkage = shrinkage,
      cv.folds = cv.folds,
      weights=GBM.weights,
      var.monotone = var.monotone,
      interaction.depth = interaction.depth,
      train.fraction = train.fraction,
      keep.data = keep.data,
      verbose = verbose,
      class.stratify.cv = class.stratify.cv,
      n.cores = n.cores)
    
    b12<-summary(boosttrans12)
    b13<-summary(boosttrans13)
    b23<-summary(boosttrans23)
    varimp<-list(boost.residual12=b12, boost.residual13=b13, boost.residual23=b23)
  }
  
  # RF
  if(method==meth[2]){
    library("randomForest")
    a12<-NULL
    a13<-NULL
    a23<-NULL
    RFtrans12<-randomForest(residual[[1]]~ ., data=data12, 
                            ntree=ntree,
                            na.action=na.action,
                            mtry=mtry,
                            replace=replace, 
                            strata=strata,
                            nodesize = nodesize,
                            maxnodes = maxnodes,
                            importance=TRUE, 
                            localImp=localImp, 
                            nPerm=nPerm,...)
    RFtrans13<-randomForest(residual[[2]]~ ., data=data13, 
                            ntree=ntree,
                            na.action=na.action,
                            mtry=mtry,
                            replace=replace, 
                            strata=strata,
                            nodesize = nodesize,
                            maxnodes = maxnodes,
                            importance=TRUE, 
                            localImp=localImp, 
                            nPerm=nPerm,...)
    RFtrans23<-randomForest(residual[[3]]~ ., data=data23, 
                            ntree=ntree,
                            na.action=na.action,
                            mtry=mtry,
                            replace=replace, 
                            strata=strata,
                            nodesize = nodesize,
                            maxnodes = maxnodes,
                            importance=TRUE, 
                            localImp=localImp, 
                            nPerm=nPerm,...)
    a12<-as.matrix(sort(RFtrans12$importance[,2], decreasing = T))
    a13<-as.matrix(sort(RFtrans13$importance[,2], decreasing = T))
    a23<-as.matrix(sort(RFtrans23$importance[,2], decreasing = T))
    
    varimp<-list(RF.varimp12(IncNodePurity)=a12,RF.varimp13(IncNodePurity)=a13,RF.varimp23(IncNodePurity)=a23)
  }
  
  # Neural network
  if(method==meth[3]){
    library("nnet")
    library("NeuralNetTools")
    
    AN12 <- nnet(residual[[1]] ~. , data = data12, 
                 size=size, 
                 weights=ANN.weights,
                 contrasts=contrasts,
                 linout = linout, 
                 entropy = entropy, 
                 softmax = softmax,
                 censored = censored, 
                 skip = skip, 
                 rang = rang, 
                 decay = decay,
                 maxit = maxit, 
                 Hess = Hess, 
                 trace = trace, 
                 MaxNWts = MaxNWts,
                 abstol = abstol, 
                 reltol = reltol, ...,
                 na.action=na.action,
                 subset=subset
                 ,...)
    ANN12<-garson(AN12,bar_plot = FALSE)
    ANNO12<-olden(AN12,bar_plot = FALSE)
    
    AN13 <- nnet(residual[[2]] ~. , data = data13, 
                 size=size, 
                 weights=ANN.weights,
                 contrasts=contrasts,
                 linout = linout, 
                 entropy = entropy, 
                 softmax = softmax,
                 censored = censored, 
                 skip = skip, 
                 rang = rang, 
                 decay = decay,
                 maxit = maxit, 
                 Hess = Hess, 
                 trace = trace, 
                 MaxNWts = MaxNWts,
                 abstol = abstol, 
                 reltol = reltol, ...,
                 na.action=na.action,
                 subset=subset
                 ,...)
    ANN13<-garson(AN13,bar_plot = FALSE)
    ANNO13<-olden(AN13,bar_plot = FALSE)
    
    AN23 <- nnet(residual[[3]] ~. , data = data23,  
                 size=size, 
                 weights=ANN.weights,
                 contrasts=contrasts,
                 linout = linout, 
                 entropy = entropy, 
                 softmax = softmax,
                 censored = censored, 
                 skip = skip, 
                 rang = rang, 
                 decay = decay,
                 maxit = maxit, 
                 Hess = Hess, 
                 trace = trace, 
                 MaxNWts = MaxNWts,
                 abstol = abstol, 
                 reltol = reltol, ...,
                 na.action=na.action,
                 subset=subset
                 ,...)
    ANN23<-garson(AN23,bar_plot = FALSE)
    ANNO23<-olden(AN23,bar_plot = FALSE)
    
    varimp<-list(ANN.residual12=ANN12, ANN.residual13=ANN13, ANN.residual23=ANN23)
  }
  
  return(varimp)
}