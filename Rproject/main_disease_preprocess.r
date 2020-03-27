#Requires  dir, fname, k, k2, confounders, pheno, indcase, indctrl, genes, genes2, expr, preprocess_param, patids. Obtained from loading a GEO dataset and assigning the right variables to the right name
#To run on a previous R version; Clean the data, use Peer and save Everything on a name,k,.Rdata file
#library(peer)
library("peer", lib.loc="/PHShome/mnm92/")
source(paste(dir,"Rproject/functions.r",sep=""))

genenames=as.character(genes);names=genes2;
indmissing=c();
if (length(confounders)){indmissing=which(apply(confounders,1,function(x)any(is.na(x)))); print(paste("Missing values for confounders removed: ",length(indmissing)))}

pheno0=rep(-1,length(pheno));pheno0[indcase]=1;pheno0[indctrl]=0; if(length(indmissing))pheno0[indmissing]=-1;
selected=which(pheno0==0 | pheno0==1);print(paste("After correcting for missingness we have",length(which(pheno0==1)),"cases and ",length(which(pheno0==0)),"controls"))

exprselected=expr[,selected];rownames(exprselected)<- rownames(expr);colnames(exprselected)<- colnames(expr)[selected]
#Remove NAs and infs
nalines=which(apply(exprselected,1,function(x)any(is.na(x) | x==Inf| x==-Inf)));
if (length(nalines)){exprselected=exprselected[-nalines,];genenames=genenames[-nalines];names=names[-nalines];print(paste("Genes NA removed: ",length(nalines)));}
#remove lowly expressed and unvarying. Removing genes that are lowly expressed in more than 90% of participants
notexpressed=which(apply(exprselected,1,function(x)mean(x<preprocess_param[[1]])>preprocess_param[[2]]) | apply(exprselected,1,var)==0);
if (length(notexpressed)){exprselected=exprselected[-notexpressed,];genenames=genenames[-notexpressed];names=names[-notexpressed];print(paste("Low expression genes removed: ",length(notexpressed)));}
if (preprocess_param[[3]]){exprselected=log(exprselected+preprocess_param[[3]]);}
exprscaled=t(scale(t(exprselected)));rownames(exprscaled)<- rownames(exprselected);colnames(exprscaled)<- colnames(exprselected);
exprselected=NULL;

confounders_final=confounders[selected,];
if (k2) confounders_final=data.frame(confounders_final,prcomp(t(exprscaled))$x[,1:k2]);#ADD top k2 PCs as confounders


model = PEER()
PEER_setPhenoMean(model,t(as.matrix(exprscaled)))
PEER_setNk(model,k)
if (length(confounders)){PEER_setCovariates(model, as.matrix(confounders_final))}
PEER_setAdd_mean(model,TRUE)
#PEER_setNmax_iterations(model,500)
PEER_update(model)

residuals = PEER_getResiduals(model)
factors = PEER_getX(model)
weights = PEER_getW(model)


factornames=c();
if (ncol(factors)){
uselessfactors=which(apply(factors,2,var)==0);
if (length(uselessfactors)){factors=factors[,-uselessfactors]; weights=weights[,-uselessfactors];}#Remove bias term and unvarying factors
if (ncol(factors)){factornames=paste("P_factor",1:ncol(factors),sep="_");}
print(dim(factors))
}
finalexpr=t(cbind(residuals,factors)); genenames=c(genenames,factornames);names=c(names,factornames); #Add PEER factors
finalpheno=pheno0[selected];finalpatids=patids[selected];

#Showing gene correlations before and after pre-processing
subsample_feat=min(10000,nrow(finalexpr)-ncol(factors));
cor_analyze=min(1000000,subsample_feat^2);#0 if i don't want to do it
if (cor_analyze){#Analysing correlation
  ind_feat=sample(1:(nrow(finalexpr)-ncol(factors)),subsample_feat)  # Only 
  cor_post= sample(cor(t(finalexpr[ind_feat,])),cor_analyze);
  cor_pre= sample(cor(t(exprscaled[ind_feat,])),cor_analyze);
  cor_quant=rbind(get_cor_quartiles(cor_pre),get_cor_quartiles(cor_post))
  pdf(paste(dir,fname,"_cor_pre.pdf",sep=""));hist(cor_pre);dev.off();
  pdf(paste(dir,fname,"_cor_post.pdf",sep=""));hist(cor_post);dev.off();
  print(cor_quant)
}  

exprscaled=NULL;residuals = NULL;cor_pre=NULL;cor_post=NULL;ind_feat=NULL;#Freeing up space
save.image(paste(dir,fname,".RData",sep=""))
