
eval_acc=function(y,y1){
  fg=y1[which(y==1)];bg= y1[which(y==0)];roc_obj <- roc(y, y1,quiet=TRUE)
  roc<-roc.curve(scores.class0 = fg, scores.class1 = bg)
  pr<-pr.curve(scores.class0 = fg, scores.class1 = bg)
  eval=c(auc(roc_obj),roc$auc,pr$auc.integral); names(eval)<- c("auc","auc2","aupr")
  return(eval)
}

plot_auc_roc=function(y,y1){#y1 contains multiple methods 
  palet=c("red","blue","gray","brown","yellow")
  for (i in 1:ncol(y1)){
    roc_obj <- roc(y, y1[,i],quiet=TRUE,curve=TRUE)
    roc_rose <- plot(roc, color = palet[i],add=(i>1),auc.main=(i>1),rand.plot=(i==1))
  }
  return(roc_rose)  
}

plot_auc_pproc=function(y,y1){#y1 contains multiple methods 
  fg=y1[which(y==1),];bg= y1[which(y==0),];
  palet=c("red","blue","gray","brown","yellow")
  for (i in 1:ncol(y1)){
    roc<- roc.curve(scores.class0 = fg[,i], scores.class1 = bg[,i],curve=TRUE,max.compute=TRUE)
    roc_rose <- plot(roc, color = palet[i],add=(i>1),auc.main=(i>1),rand.plot=(i==1))
  }
  return(roc_rose)  
}



train_classifier=function(x,y,train_ids,val_ids,method,standardize=TRUE,low_l=-Inf,balanced=FALSE){
  if (method=="rf")return(train_classifier_rf(x,y,train_ids,val_ids,standardize=standardize,low_l=low_l,balanced=balanced))
  if (method=="lda")return(train_classifier_lda(x,y,train_ids,val_ids,standardize=standardize,low_l=low_l,balanced=balanced))
  if (method=="qda")return(train_classifier_qda(x,y,train_ids,val_ids,standardize=standardize,low_l=low_l,balanced=balanced))
  if (method=="pfda")return(train_classifier_pfda(x,y,train_ids,val_ids,standardize=standardize,low_l=low_l,balanced=balanced))
  #glmnet by default
  return(train_classifier_glmnet(x,y,train_ids,val_ids,standardize=standardize,low_l=low_l,balanced=balanced))
}  

train_classifier_glmnet=function(x,y,train_ids,val_ids,standardize=TRUE,low_l=-Inf,balanced=FALSE){
  train_ids2=train_ids;if (balanced){train_ids2=c(train_ids[which(y[train_ids]==1)],sample(train_ids[which(y[train_ids]==0)],length(which(y[train_ids]==1))))}
  if (low_l==-Inf)res=glmnet(t(x[,train_ids2]),y[train_ids2],family="binomial",standardize=standardize)#Was lower.limits=lower.limits, but raises error
  if (low_l==0)res=glmnet(t(x[,train_ids2]),y[train_ids2],family="binomial",standardize=standardize,lower.limits=low_l);
  y1=predict(res,t(x),type="response")#predict for every lambda
  evals=apply(y1[val_ids,],2,function(x)eval_acc(y[val_ids],x))
  bestlambda=which.max(evals[3,])
  return(list(res=res,y1=y1,evals=evals,bestlambda=bestlambda,bestmodel=y1[,bestlambda]))
}

train_classifier_rf=function(x,y,train_ids,val_ids,standardize=TRUE,low_l=-Inf,balanced=FALSE){
  train_ids2=train_ids;if (balanced){train_ids2=c(train_ids[which(y[train_ids]==1)],sample(train_ids[which(y[train_ids]==0)],length(which(y[train_ids]==1))))}
  res=randomForest(t(x[,train_ids2]),as.factor(y[train_ids2]));
  y1=predict(res,t(x),type="prob")[,2]#predict for every lambda
  evals=eval_acc(y[val_ids],y1[val_ids])
  #bestlambda=which.max(evals[3,])
  return(list(res=res,y1=y1,evals=evals,bestmodel=y1))#,evals=evals,bestlambda=bestlambda,bestmodel=y1[,bestlambda]))
}

train_classifier_lda=function(x,y,train_ids,val_ids,standardize=TRUE,low_l=-Inf,balanced=FALSE){
  train_ids2=train_ids;if (balanced){train_ids2=c(train_ids[which(y[train_ids]==1)],sample(train_ids[which(y[train_ids]==0)],length(which(y[train_ids]==1))))}
  res=lda(t(x[,train_ids2]),as.factor(y[train_ids2]));
  y1=predict(res,t(x))$posterior[,2]#predict for every lambda
  evals=eval_acc(y[val_ids],y1[val_ids])
  #bestlambda=which.max(evals[3,])
  return(list(res=res,y1=y1,evals=evals,bestmodel=y1))#,evals=evals,bestlambda=bestlambda,bestmodel=y1[,bestlambda]))
}

train_classifier_qda=function(x,y,train_ids,val_ids,standardize=TRUE,low_l=-Inf,balanced=FALSE){
  train_ids2=train_ids;if (balanced){train_ids2=c(train_ids[which(y[train_ids]==1)],sample(train_ids[which(y[train_ids]==0)],length(which(y[train_ids]==1))))}
  res=qda(t(x[,train_ids2]),as.factor(y[train_ids2]));
  y1=predict(res,t(x))$posterior[,2]#predict for every lambda
  evals=eval_acc(y[val_ids],y1[val_ids])
  #bestlambda=which.max(evals[3,])
  return(list(res=res,y1=y1,evals=evals,bestmodel=y1))#,evals=evals,bestlambda=bestlambda,bestmodel=y1[,bestlambda]))
}

train_classifier_pfda=function(x,y,train_ids,val_ids,standardize=TRUE,low_l=-Inf,balanced=FALSE){
  train_ids2=train_ids;if (balanced){train_ids2=c(train_ids[which(y[train_ids]==1)],sample(train_ids[which(y[train_ids]==0)],length(which(y[train_ids]==1))))}
  res=pfda(t(x[,train_ids2]),as.factor(y[train_ids2]));
  y1=predict(res,t(x))$P[,2]#predict for every lambda
  evals=eval_acc(y[val_ids],y1[val_ids])
  #bestlambda=which.max(evals[3,])
  return(list(res=res,y1=y1,evals=evals,bestmodel=y1))#,evals=evals,bestlambda=bestlambda,bestmodel=y1[,bestlambda]))
}
