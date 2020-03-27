rmax=0.3;dataname="ovarian_mirna100";topfs=c(10,20,50,100,200,300,500);
method=""
dirx="Downloads/aziz_test/"
load(paste0("Desktop/testing_expression/prevres/",dataname,".RData2"))
library(glmnet)
library(caret)
library(pROC)
library(PRROC)
library(randomForest)
library(probFDA)
source(paste0(dirx,"Rproject/functions.r"))
source(paste0(dirx,"Rproject/classification/classification_functions.r"))

#Variable used: finalexpr,finalpheno,w0,genenames,names

candidates_breastcancer=c("miR-1246", "miR-1307-3p", "miR-4634", "miR-6861-5p", "miR-6875-5p");
candidates_ovcancer=c("miR-320a", "miR-665", "miR-3184-5p", "miR-6717-5p", "miR-4459", "miR-6076", "miR-3195", "miR-1275", "miR-3185", "miR-4640-5p");

indpaper=which(substr(genenames,5,100) %in% c("miR-320a", "miR-665", "miR-3184-5p", "miR-6717-5p", "miR-4459", "miR-6076", "miR-3195", "miR-1275", "miR-3185", "miR-4640-5p"));#Ovariancancer
#indpaper=which(substr(genenames,5,100) %in% c("miR-1246", "miR-1307-3p", "miR-4634", "miR-6861-5p", "miR-6875-5p"))

binarize=TRUE;makepositive=F;low_limits=0;
binarize=F;makepositive=T;low_limits=0;
binarize=F;makepositive=F;low_limits=-Inf;


nrep=1000000;kfold=5;accuracythresh=1/nrep;
testprop=0.15;
test_n=round(testprop*length(which(finalpheno==1)))
test_ids=c(sample(which(finalpheno==1),test_n),sample(which(finalpheno==0),test_n))
cv_ids=(1:length(finalpheno))[-test_ids] #1:ncol(finalexpr)#Cross validation ids

#Doing Cross-Validation
flds <- createFolds(cv_ids, k = kfold, list = TRUE, returnTrain = FALSE)
resfold=list(); length(resfold)<- kfold; pred=array(0,dim=c(length(finalpheno),6,length(topfs),kfold));valpred=pred;

for (k in 1:length(flds)){
  
  val_ids=cv_ids[flds[[k]]];
  train_ids=cv_ids[unlist(flds[-k])];
  #Running our test and Limma
  finalpheno_t=finalpheno[train_ids];finalexpr_t=finalexpr[,train_ids];
  ptm= proc.time()
  pval_wilc=apply(finalexpr_t,1,function(e)wilcox.test(e[which(finalpheno_t==1)],e[which(finalpheno_t==0)])$p.value);
  design <- model.matrix(~ as.factor(finalpheno_t));#colnames(design) <- c("cases", "controls");
  fit <- lmFit(finalexpr_t, design=design); fit2 <- eBayes(fit)
  results2=topTable(fit2,coef=2,n=Inf,sort="none")
  pval_limma=results2[,4];fdr_limma=results2[,5];
  res_esa=list();length(res_esa)<- nrow(finalexpr_t);
  for (i in 1:nrow(finalexpr_t)){e=finalexpr_t[i,];res_esa[[i]]<- aziz.test(finalpheno_t,e,w0,rep=0,ignoremax=0);}
  res_ref=reformat_results(res_esa);
  calibration=calibrate_test(finalpheno_t,w0,rep=nrep,doall=TRUE,ignoremax=0);
  esa_p=get_calibrated_pvalues(calibration,res_ref$es);
  esa_cor=correct_zero_pvalues(calibration,log((nrep:1)/nrep),res_ref$es,esa_p,accuracythresh,maxpvaltrain=0.05)
  fdr_esa=p.adjust(esa_cor, method="fdr");fdr_w=p.adjust(pval_wilc, method="fdr");
  print(proc.time()-ptm)
  show_res2=data.frame(genenames,names,fdr_esa,fdr_limma,res_ref,pval2=esa_p,approx=esa_cor,limma=pval_limma,fdr_w,wilc=pval_wilc);
  partition=data.frame(train=rep(0,length(finalpheno)),val=0,test=0);partition[train_ids,1]=1;partition[val_ids,2]=1;partition[test_ids,3]=1;
  write.table(show_res2,paste0("~/Downloads/classification/CV_",dataname,"_res_",k,".txt"))
  write.table(partition,paste0("~/Downloads/classification/CV_",dataname,"_par_",k,".txt"))
  
#}  
for (k in 1:length(flds)){
  show_res2=read.table(paste0("~/Downloads/classification/CV_",dataname,"_res_",k,".txt"),header=T)
  partition=read.table(paste0("~/Downloads/classification/CV_",dataname,"_par_",k,".txt"),header=T); train_ids=which(partition$train==1);val_ids=which(partition$val==1);test_ids=which(partition$test==1)

#  finalexprold=finalexpr#Very weird adding this line preserves finalexpr , instead of it taking the value of finalexpr2 later (only happens if I am loading multiple in for loop) 
  #feature Selection based on previous testing
  n_factors=1;addedfactors=c();if(n_factors>0)addedfactors=length(rg)+c(1:n_factors)
  for (j in 1:length(topfs)){
    topf=topfs[j]; 
    if (length(which(show_res2$r[rg]<rmax))>=topf){
    if (topf==0){ind=c();ind0=c();ind1=c();}else{
    ord=order(show_res2$fdr_esa[rg]);ord2=order(show_res2$fdr_limma[rg]);
    ind=ord[which(show_res2$r[ord]<rmax)][1:topf];
    ind1=ord[1:topf];ind0=ord2[1:topf];
    }
  #transforming the data with the selected features (based on our test)
    finalexpr2=finalexpr;
    for (i in ind){
      indab=which_aberrant(finalexpr[i,],finalexpr[i,train_ids],show_res2$es_sign[i],show_res2$es_indmax[i])
      finalexpr2[i,-indab]=0;if (binarize)finalexpr2[i,indab]=1; if(makepositive)finalexpr2[i,indab]=abs(finalexpr2[i,indab]);
    }
    
    #Build classifiers (only uses train_ids to train)
    yn=train_classifier(finalexpr2[c(ind,addedfactors),], finalpheno,train_ids,val_ids,method,low_l=low_limits)
    y=train_classifier(finalexpr[c(ind,addedfactors),], finalpheno,train_ids,val_ids,method)
    y1=train_classifier(finalexpr[c(ind1,addedfactors),], finalpheno,train_ids,val_ids,method)
    y0=train_classifier(finalexpr[c(ind0,addedfactors),], finalpheno,train_ids,val_ids,method)
    if (j==1)yall=train_classifier(finalexpr[c(rg,addedfactors),], finalpheno,train_ids,val_ids,method)#all features, does not depend on j
    if (j==1)ypaper=train_classifier(finalexpr[c(indpaper,addedfactors),], finalpheno,train_ids,val_ids,method)#all features, does not depend on j
    resfold[[k]]<- list(yn,y,y1,y0,yall,ypaper);
    pred[,,j,k]=cbind(yn$bestmodel,y$bestmodel,y1$bestmodel,y0$bestmodel,yall$bestmodel,ypaper$bestmodel);
    valpred[val_ids,,j,k]=pred[val_ids,,j,k];
  }else{print(paste0(j, "failed"))}}
}

finalpred=apply(pred,c(1,2,3),mean);finalvalpred=apply(valpred,c(1,2,3),sum);
apply(finalpred[test_ids,,],c(2,3),function(x)eval_acc(finalpheno[test_ids],x))
apply(finalvalpred[-test_ids,,],c(2,3),function(x)eval_acc(finalpheno[-test_ids],x))

save.image(paste0("~/Downloads/classification/all_",dataname,"_",rmax,"_",method,n_factors,".RData"))

j=3;
plot_auc_pproc(finalpheno[test_ids],finalpred[test_ids,,j])




#ovarian cancer data titles processing
title=sapply(as.character(gse$title),function(x){t=strsplit(x," ")[[1]];return(paste(t[-length(t)],collapse=" "))})
ov_cancer_id=which(title=="Ovarian Cancer");ctrl_id=which(title=="non-Cancer");
ov_tumor_other_id=which(title %in% c("Borderline Ovarian Tumor","OV_others"));ov_dis_id=which(title=="Benign Ovarian Disease");
other_cancer_id=(1:length(title))[-c(ctrl_id,ov_cancer_id,ov_tumor_other_id,ov_dis_id)];
ctrl_choice=ctrl_id
finalexpr=finalexpr[,c(ov_cancer_id,ctrl_choice)];finalpheno=c(rep(1,length(ov_cancer_id)),rep(0,length(ctrl_choice)))

#observations:
Breast cancer: 
I can get to auc of 0.96-0.97 with my method. while all others are around auc random, previous probes are at auc 0.55 aupr 0.54. 
Adding just one Peer factor bring all other methods to 0.94 auc,0.96 aupr. My method is still better at auc/aupr 0.98
Using 10 PEER factors, even alone, leads to full classification, 5 factors gets all methods to 0.996 auc/aupr (factors alone).

Ovarian cancer: vs ctrl 
I can get to 0.94-0.95 auc/aupr, other methods are around 0.59 auc/ 0.63-0.69 aupr, previous probes are random
Adding just one Peer factor bring all other methods to 0.92-0.93 auc,0.90 aupr. My method is still better at auc/aupr 0.95-0.96. With 3 factors we are at 0.98 and 0.96, 0.96 explained by the factors alone
25 factors give full classification on their own.

Ovarian cancer vs other cancers:
  Best we can do is auc 0.64 and aupr 0.61 . same with limma features. no transformation. Only 5 miRNA had FDR<0.1
  with 25 factors we can get auc of 0.94-0.96
