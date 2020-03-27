#!/hpf/tools/centos6/R/3.2.3/bin/Rscript
if (length(commandArgs(TRUE))<4){print("Error: incorrect number of arguments");if(NA)print("Error");}
dir0=commandArgs(TRUE)[1];
name=commandArgs(TRUE)[2];
k0=as.numeric(commandArgs(TRUE)[3]);
exact0=as.numeric(commandArgs(TRUE)[4]);


local=TRUE;
if(local){
  dir0="/home/aziz/Desktop/aziz/testing_expression/";
  name="Alz_rnaseq"
  k0=30; exact0=0;
}

preprocessing=TRUE;
if (preprocessing){
load(paste(dir0,name,"/",name,"30",".RData",sep=""))
  dir0=commandArgs(TRUE)[1];
  name=commandArgs(TRUE)[2];
  k0=as.numeric(commandArgs(TRUE)[3]);
  exact0=as.numeric(commandArgs(TRUE)[4]);
k=k0;exact=exact0;dir=dir0;

library(GEOquery)
source(paste(dir,"functions.r",sep=""))

nrep=10000000;accuracythresh=1/nrep;#recommended <0.05; 0 only fixes 0 permutation pval
perm=FALSE;
dir.create(paste(dir,name,"/",sep=""))

library(peer)

#expr=exprs(gse)
genenames=as.character(genes);
#names=featureNames(gse)

indmissing=c();
if (length(confounders)){indmissing=which(apply(confounders,1,function(x)any(is.na(x)))); print(paste("Missing values for confounders removed: ",length(indmissing)))}

pheno0=rep(-1,length(pheno));pheno0[indcase]=1;pheno0[indctrl]=0; if(length(indmissing))pheno0[indmissing]=-1;
selected=which(pheno0==0 | pheno0==1);print(paste("After correcting for missingness we have",length(which(pheno0==1)),"cases and ",length(which(pheno0==0)),"controls"))

exprselected=expr[,selected];rownames(exprselected)<- rownames(expr);colnames(exprselected)<- colnames(expr)[selected]
#Remove NAs and infs
nalines=which(apply(exprselected,1,function(x)any(is.na(x) | x==Inf| x==-Inf)));
if (length(nalines)){exprselected=exprselected[-nalines,];genenames=genenames[-nalines];names=names[-nalines];print(paste("Genes NA removed: ",length(nalines)));}
#remove lowly expressed and unvarying. Removing genes that are lowly expressed in more than 90% of participants
minexpr=0.2;perc_minexpr=0.9
notexpressed=which(apply(exprselected,1,function(x)mean(x<minexpr)>.9) | apply(exprselected,1,var)==0);
if (length(notexpressed)){exprselected=exprselected[-notexpressed,];genenames=genenames[-notexpressed];names=names[-notexpressed];print(paste("Low expression genes removed: ",length(notexpressed)));}
exprscaled=t(scale(t(exprselected)));rownames(exprscaled)<- rownames(exprselected);colnames(exprscaled)<- colnames(exprselected);
exprselected=NULL;

model = PEER()
PEER_setPhenoMean(model,t(as.matrix(exprscaled)))
PEER_setNk(model,k)
if (length(confounders)){PEER_setCovariates(model, as.matrix(confounders[selected,]))}
PEER_setAdd_mean(model,TRUE)
#PEER_setNmax_iterations(model,500)
PEER_update(model)

residuals = PEER_getResiduals(model)
#factors = PEER_getX(model)
#weights = PEER_getW(model)

finalexpr=t(residuals)#only considers selected indiv
finalpheno=pheno0[selected]

#Showing gene correlations before and after pre-processing
subsample_feat=min(10000,nrow(finalexpr));
cor_analyze=min(1000000,subsample_feat^2);#0 if i don't want to do it
if (cor_analyze){#Analysing correlation
 ind_feat=sample(1:nrow(finalexpr),subsample_feat)  
 cor_post= sample(cor(t(finalexpr[ind_feat,])),cor_analyze);
 cor_pre= sample(cor(t(exprscaled[ind_feat,])),cor_analyze);
 cor_quant=rbind(get_cor_quartiles(cor_pre),get_cor_quartiles(cor_post))
 pdf(paste(dir,name,"/",name,"_cor_pre.pdf",sep=""));hist(cor_pre);dev.off();
 pdf(paste(dir,name,"/",name,k,"_cor_post.pdf",sep=""));hist(cor_post);dev.off();
 print(cor_quant)
}  

exprscaled=NULL;residuals = NULL;factors = NULL;weights =NULL;cor_pre=NULL;cor_post=NULL;ind_feat=NULL;#Freeing up space
save.image(paste(dir,name,"/",name,k,".RData",sep=""))

}else{
load(paste(dir0,name,"/",name,k0,".RData",sep=""))
k=k0;exact=exact0;dir=dir0;
source(paste(dir,"functions.r",sep=""))
}
#source(paste(dir,"functions.r",sep=""))
print("Running differential expression tests...")

pval_t=apply(finalexpr,1,function(e)t.test(e[which(finalpheno==1)],e[which(finalpheno==0)])$p.value);
fdr_t=p.adjust(pval_t, method="fdr");
library(limma)
design <- model.matrix(~ as.factor(finalpheno));#colnames(design) <- c("cases", "controls");
fit <- lmFit(finalexpr, design=design)
fit2 <- eBayes(fit)
results2=topTable(fit2,coef=2,n=Inf,sort="none")
fdr_limma=results2[,5]
pval_limma=results2[,4]

print("Computing enrichment scores...")

w0=NULL;
#closed form es normalization ,no permutations 1min
ptm= proc.time()
res_esa=list();length(res_esa)<- nrow(finalexpr);
for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa[[i]]<- psea(finalpheno,e,w0,rep=0);}
res_ref=reformat_results(res_esa);
print(proc.time()-ptm)

#First approximation: Estimate from 1 set of permutations
print("Computing first pvalues approximation using a single set of permutations for all genes...");ptm= proc.time()
ptm= proc.time()
calibration=calibrate_psea(finalpheno,w0,rep=nrep,doall=TRUE);
esa_m=get_calibrated_pvalues(calibration,res_ref$es);
print(proc.time()-ptm)

#Exact computation : Full permutation on every gene 2h
esa_p=esa_m;
if(exact){
print("Computing exact p-values with permutations for every gene...")
res_esa_p=list();length(res_esa_p)<- nrow(finalexpr); 
ptm= proc.time()
for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa_p[[i]]<- psea(finalpheno,e,w0,rep=nrep);}
print(proc.time()-ptm)
es_p=rep(0,nrow(finalexpr));esa_p=rep(1,nrow(finalexpr));for (i in 1:nrow(finalexpr)){es_p[i]=res_esa_p[[i]]$es;esa_p[i]=(res_esa_p[[i]])$pval;}
}

#modeling p-values to compute good approximations
esa_cor=correct_zero_pvalues(calibration,log((nrep:1)/nrep),res_ref$es,esa_p,accuracythresh,maxpvaltrain=0.05)
if (exact)esa_cor=correct_zero_pvalues(res_ref$es,log(esa_p),res_ref$es,esa_p,accuracythresh,maxpvaltrain=0.05)
#plot(log10(predict_pval(res_ref$es)),log10(esa_p));lines(c(0,-20),c(0,-20))#training
pdf(paste(dir,name,"/",name,k,"_pvalues_predicted_versus_estimated.pdf",sep=""));
plot(-log10(esa_p),-log10(predict_pval(res_ref$es)), xlab="Negative Log10 Estimated p-values",ylab="Negative Log10 predicted p-values",main="Predicted p-values quality");lines(c(0,20),c(0,20));dev.off();
pdf(paste(dir,name,"/",name,k,"_pvalues_function_of_scores.pdf",sep=""));
plot((1:1000)/100,-log10(predict_pval((1:1000)/100)), xlab="Max Enrichment Score ",ylab="Negative Log10 p-values",main="Predicted p-values function");dev.off();

fdr_esa=p.adjust(esa_cor, method="fdr");
#print(genenames[which(fdr_esa<0.2)])
#print(genenames[which(fdr_limma<0.2)])

signifthresh=0.05/length(genenames)
indg=which(esa_p<= signifthresh)
for (i in indg){e=finalexpr[i,];print(paste(genenames[i],res_esa[[i]]$summary,format(predict_pval(res_esa[[i]]$es),digits=3)));}
for (i in indg){e=finalexpr[i,];print(paste(genenames[i],res_esa[[i]]$details));}
#for (i in indg){print("");expr_plot(finalexpr[i,],finalpheno,paste(dir,name,"/",name,k,"_",genenames[i],sep=""))}
print("Number of genes found with FDR<0.2");printvenn(which(fdr_esa<0.2),which(fdr_limma<0.2));
print("Number of genes found (Bonferroni)");printvenn(which(esa_cor<signifthresh),which(pval_limma<signifthresh));

write.table(data.frame(gname=genenames,score=res_ref$es,esa_pval=esa_p,approx=esa_cor,limma=pval_limma,fdr=fdr_esa, fdrlimma=fdr_limma),paste(dir,name,"/",name,k,".txt",sep=""),row.names=FALSE,sep="\t")

expr=NULL;exprscaled=NULL;residuals = NULL;factors = NULL;weights =NULL;#Freeing up space
save.image(paste(dir,name,"/",name,k,".RData2",sep=""))

indg=which(fdr_esa<0.1 | fdr_limma<0.1);
col=c("blue","green");


plot_graph=function(indg,suff){
pdf(paste(dir,name,"/",name,k,suff,".pdf",sep=""))
plot(-log10(pval_limma[indg]), -log10(esa_cor[indg]), xlab="Limma Negative Log p-value",ylab="Aberration Test Negative Log p-value",main="Significant genes",col=col[res_ref$es_sign[indg]],
xlim=minrange(c(-0.55,3)+range(-log10(esa_cor[indg]))) )
text(-log10(pval_limma[indg]), -log10(esa_cor[indg]),labels=genenames[indg],cex=0.75,pos=3,col=col[res_ref$es_sign[indg]]);lines(c(-1,20),c(-1,20));
lines(c(-1,20),c(-log10(signifthresh),-log10(signifthresh)),col="red");lines(c(-log10(signifthresh),-log10(signifthresh)),c(-1,20),col="red");
legend("bottomright",c("Under-expressed","Over-expressed"),fill=col)
dev.off()
}

indg=which(fdr_esa<0.1 | fdr_limma<0.1);if(length(indg))plot_graph(indg,"_fdr01")
indg=which(fdr_esa<0.2 | fdr_limma<0.2);if(length(indg))plot_graph(indg,"_fdr02")
coord=read.table(paste(dir,"gene_coordinates.txt",sep=""));
library(qqman)
manhattan_plot(coord[,c(5,2,3)],genenames,data.frame(Aberration_test=esa_cor,Limma=pval_limma),paste(dir,name,"/",name,k,"_manhattan.pdf",sep=""),annotatePval=0.001)
pdf(paste(dir,name,"/",name,k,"_qq.pdf",sep=""));qq(esa_cor);dev.off()
plot_2g(finalexpr,res_esa, indg[1],indg[2],paste(dir,name,"/",name,k,sep=""),genenames)

pdf(paste(dir,name,"/","qq.pdf",sep=""));qq(esa_p);dev.off()

respca=prcomp(finalexpr,center = TRUE,scale. = TRUE)
pdf(paste("eigen_",name,".pdf",sep=""));plot(1:length(respca$sdev),respca$sdev);dev.off()

if (perm){
finalpheno2=sample(finalpheno,length(finalpheno))
ptm= proc.time()
res_esa2=list();length(res_esa2)<- nrow(finalexpr);
for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa2[[i]]<- psea(finalpheno2,e,w0,rep=0);}
print(proc.time()-ptm)
#es2=rep(0,nrow(finalexpr));for (i in 1:nrow(finalexpr)){es2[i]=(res_esa2[[i]])$es;}
res_ref2=reformat_results(res_esa2);
esa_m2=get_calibrated_pvalues(calibration,res_ref2$es);
#pdf(paste(dir,name,"/","qq_perm.pdf",sep=""));qq(esa_m2);dev.off()
esa_cor2=correct_zero_pvalues(calibration,log((nrep:1)/nrep),res_ref2$es,esa_m2,accuracythresh,maxpvaltrain=0.05)
fdr_esa2=p.adjust(esa_cor2, method="fdr");
}


compositeplot=FALSE; #All significant genes; One example gene (gsname); qq plot; qq plot of permuted
gsname="CRBN"#"UQCRH"
if (compositeplot){
#indg=which(fdr_esa<0.1 | fdr_limma<0.1);
pdf(paste(dir,name,"/","res_all_",gsname,".pdf",sep=""));
par(mfrow = c(2,2))
col=c("blue","green");
plot(-log10(pval_limma[indg]), -log10(esa_cor[indg]), xlab="Limma Negative Log p-value",ylab="Aberration Test Negative Log p-value",main="Significant genes",col=col[res_ref$es_sign[indg]],
     xlim=minrange(c(-0.55,3)+range(-log10(esa_cor[indg]))) )
text(-log10(pval_limma[indg]), -log10(esa_cor[indg]),labels=genenames[indg],cex=0.6,pch=8,pos=3,col=col[res_ref$es_sign[indg]]);lines(c(-1,20),c(-1,20));
lines(c(-1,20),c(-log10(signifthresh),-log10(signifthresh)),col="red");lines(c(-log10(signifthresh),-log10(signifthresh)),c(-1,20),col="red");
legend("topleft",c("Under-expressed","Over-expressed"),fill=col,cex=0.75)
col=c("black","red");e=finalexpr[which(genenames==gsname & ( (1:length(genenames)) %in% indg))[1],];
plot(density(e[which(finalpheno==0)]),col=col[1],main=paste(gsname,"expression levels"));lines(density(e[which(finalpheno==1)]),col=col[2]);
legend("topleft",c("controls","cases"),fill=col)
qq(esa_p);
qq(esa_m2,ylim=c(0,7.5));
dev.off()
}

#sort -g -k4 CD/CD30.txt | awk '{if (NR>1)print substr($1,2,length($1)-2)}' - | head
#awk '{t=0.000002;if (NR>1 && ($6<0.1 && $7>0.1 || $4<t && $5>t)&& $4/($5+0.0000000001*t) <0.1 )print substr($1,2,length($1)-2)}' Alzheimer_ctrl/Alzheimer_ctrl30.txt > Alzheimer_ctrl/novel_genes.txt
x=esa_cor/pval_limma
y=order(fdr_esa)[1:length(which(fdr_esa<0.1))]
genenames[y[ which(x[y]<0.1 & pval_limma[y]>signifthresh) ]]#genenames[y[ which(x[y]<0.1 & fdr_limma[y]>0.1) ]]

#Alzheimer "ATP6V1D" "CRBN"    "POMP"    "PPCDC"   "FBP1"    "FDXR"        "GPER"   "DISC1"   "LRP3"    "TLR2"    "UIMC1"   "TMEM188"
#ctrl "EIF1AY"   "GNL2"     "LY96"     "CARD14"   "ITGB3BP"  "PSMA3"   "SLC39A8"  "RPS27L"   "HLA-A"    "EIF1AY"   "CD69"     "DCP2"    
#... "MAPKAPK3" "HSPE1"    "HAT1"     "PDCD10"   "SNORD33"  "ZNF559"  
#... "C6orf66"  "GRK5"     "FBXW8" 

