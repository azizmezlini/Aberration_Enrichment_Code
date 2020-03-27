
#combine Alzheimer data:
library(GEOquery)
dir="/home/aziz/Desktop/aziz/testing_expression/"
source("/home/aziz/Desktop/aziz/testing_expression/functions.r")

library(peer)
#options('download.file.method.GEOquery'='auto')
gse=getGEO(filename="/home/aziz/Desktop/aziz/testing_expression/expression-data/Alzheimer-GSE63063-GPL10558_series_matrix.txt.gz")
gse2=getGEO(filename="/home/aziz/Desktop/aziz/testing_expression/expression-data/Alzheimer-GSE63063-GPL6947_series_matrix.txt.gz")

interprobes=intersect(fData(gse)[,1],fData(gse2)[,1]);
match1=match(interprobes,fData(gse)[,1]);match2=match(interprobes,fData(gse2)[,1])
genes=fData(gse)[match1,13];genenames=as.character(genes);
expr=cbind(exprs(gse)[match1,],exprs(gse2)[match2,])

pheno=c(as.character(gse$characteristics_ch1),as.character(gse2$characteristics_ch1))
getethnicity=function(x){if(grepl("*British*",x) |grepl("*Western European*",x))return(1); if(grepl("*Other Caucasian*",x))return(0);return(NA)}
getage=function(x){spl=unlist(strsplit(x," ")); return(as.numeric(spl[2])); }
confounders=data.frame(cohort=as.factor(c(rep(1,ncol(exprs(gse))),rep(2,ncol(exprs(gse2)) ))),
ethnicity=as.factor(sapply(c(as.character(gse$characteristics_ch1.1),as.character(gse2$characteristics_ch1.1)),getethnicity)),age=sapply(c(as.character(gse$characteristics_ch1.2),as.character(gse2$characteristics_ch1.2)),getage),gender=as.factor(c(as.character(gse$characteristics_ch1.3),as.character(gse2$characteristics_ch1.3))))
indcase=which(pheno=="status: CTL");indctrl=which(pheno=="status: MCI");#| pheno=="status: MCI"pheno=="status: AD"|
name="ctrl_MCI"
source("main_disease.r")
#phenoData(gse)varLabels(gse)
#names=featureNames(gse)


pheno0=rep(-1,length(pheno));pheno0[indcase]=1;pheno0[indctrl]=0;
selected=which(pheno0==0 | pheno0==1)
#modcombat=model.matrix(~1,data=pheno0)
#combat_data=ComBat(dat=expr,batch=confounders[,1]-1,par.prior=TRUE)
exprscaled=t(scale(t(expr[,selected])));rownames(exprscaled)<- rownames(expr);colnames(exprscaled)<- colnames(expr)[selected]
k=30
model = PEER()
PEER_setPhenoMean(model,t(as.matrix(exprscaled)))
PEER_setNk(model,k)
PEER_setCovariates(model, as.matrix(confounders[selected,]))
PEER_setAdd_mean(model,TRUE)
PEER_update(model)

factors = PEER_getX(model)
weights = PEER_getW(model)
residuals = PEER_getResiduals(model)
#peerres=peer_analysis(expr,k);


#Testing expression differences
source("/home/aziz/Desktop/aziz/testing_expression/functions.r")
finalexpr=t(residuals)#[,selected]
finalpheno=pheno0[selected]
pval_t=apply(finalexpr,1,function(e)t.test(e[which(finalpheno==1)],e[which(finalpheno==0)])$p.value);
fdr_t=p.adjust(pval_t, method="fdr");
library(limma)
design <- model.matrix(~ as.factor(finalpheno));#colnames(design) <- c("cases", "controls");
fit <- lmFit(finalexpr, design=design)
fit2 <- eBayes(fit)
results2=topTable(fit2,coef=2,n=Inf,sort="none")
fdr_limma=results2[,5]
pval_limma=results2[,4]

early=0;w0=NULL;half=FALSE;flatten=0;
#closed form es normalization ,no permutations 1min
res_esa=list();length(res_esa)<- nrow(finalexpr);
ptm= proc.time()
for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa[[i]]<- psea(finalpheno,scale(e+rnorm(length(e),0,0.00000001)),w0,early=early,rep=0,half=half,flatten=flatten);}
proc.time()-ptm
es1=rep(0,nrow(finalexpr));for (i in 1:nrow(finalexpr)){es1[i]=(res_esa[[i]])$es;}

#closed form es normalization + permutations 2h
res_esa_p=list();length(res_esa_p)<- nrow(finalexpr); 
ptm= proc.time()
for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa_p[[i]]<- psea(finalpheno,scale(e+rnorm(length(e),0,0.00000001)),w0,early=early,rep=2000000,half=half,flatten=flatten);}
proc.time()-ptm
es_p=rep(0,nrow(finalexpr));esa_p=rep(1,nrow(finalexpr));for (i in 1:nrow(finalexpr)){es_p[i]=res_esa_p[[i]]$es;esa_p[i]=(res_esa_p[[i]])$pval;}

#randomized es distribution: 
times=1;
es_rand=matrix(0,nrow(finalexpr),times);
for (k in 1:times){
finalpheno_rand=sample(finalpheno,length(finalpheno));
es_rand[,k]=apply(finalexpr,1,function(e)psea(finalpheno_rand,scale(e+rnorm(length(e),0,0.00000001)),w0,early=early,rep=2000000,half=half,flatten=flatten)$es);
}


#Assume ES gaussian
#esa1_comp=2*(1-pnorm(es1,mean=median(es_rand),sd=mad(es_rand)))
#esa2_comp=2*(1-pnorm(es2,mean=mean(es_rand),sd=sd(es_rand)))
#z=density(es_rand); mode=z$x[z$y==max(z$y)]; sdp=sqrt(sum((es_rand[which(es_rand>0)]-mode)^2)/length(which(es_rand>0)) )
#Assume ES Gumbel distributed
#beta=sd(es_rand)*sqrt(6)/pi;mu=median(es_rand)+beta*log(log(2))
#or beta=qnorm(1- (1/(exp(1)*n)))- qnorm(1- (1/n)); mu=qnorm(1- (1/n)), for a given n
#pgumbel=function(x)exp(-exp(-(x-mu)/beta))
#If we consider extreme distribution: max over dependant standardized variables,
ind=which(esa_p>0.0000005)
res_lm=lm(log(esa_p[ind])~ es1[ind]+ (es1^2)[ind]+ (es1^3)[ind])
predict_pval=function(esx)exp(cbind(1,esx,esx^2,esx^3)%*%res_lm$coefficients)
#res=lm(log(esa_p[ind])~ es1[ind] + (es1^2)[ind])
#predict_pval=function(esx)exp(cbind(1,esx,esx^2)%*%res$coefficients)
#predict_pval2=function(esx) (1-pnorm(esx,mean=mean(es_rand),sd=sd(es_rand)))#Gaussian
#predict_pval3=function(esx,div=40) (1-pnorm(esx)^(length(finalpheno)/(2*div)))#power relaxed
#predict_pval4=function(esx,div) (1- exp(-exp(-(esx-qnorm(1- (1/div)))/(qnorm(1- (1/(exp(1)*div)))- qnorm(1- (1/div))))) )#Extreme gumbel

esa=predict_pval(es1)
esa_cor=esa_p;esa_cor[which(esa_p==0)]=esa[which(esa_p==0)];
plot(log10(predict_pval(es1)),log10(esa_p));lines(c(0,-20),c(0,-20))#reconstruction
plot(log10(predict_pval((1:1000)/100)),(1:1000)/100);lines(c(0,-20),c(0,-20))

#re-evaluate potential hits without approximating the normalization
#res_esa2=res_esa; esa2=esa;
#for (i in which(esa<0.001)){e=finalexpr[i,];res_esa2[[i]]<- psea(finalpheno,scale(e+rnorm(length(e),0,0.00000001)),early=0,rep=2000000,approx=0);esa2[i]=(res_esa2[[i]])$pval;}

fdr_esa=p.adjust(esa_p, method="fdr");
print(genenames[which(fdr_t<0.15)])
print(genenames[which(fdr_esa<0.15)])


indg=which(esa<=0.000002)
for (i in indg){e=finalexpr[i,];print(paste(genenames[i],res_esa[[i]]$summary,predict_pval2(res_esa[[i]]$es)));}

name="Alz_ctl";
dir="/home/aziz/Desktop/aziz/testing_expression/"
write.table(data.frame(genenames,es1,esa,es_p,esa_p,es_rand),paste(dir,name,".txt",sep=""),row.names=FALSE,col.names=FALSE,sep="\t")
save.image(paste(dir,name,".RData",sep=""))

pdf(paste(dir,name,".pdf",sep=""))
#indg=which(esa<0.00001 | pval_t<0.00001)
indg=which(fdr_esa<0.1 | fdr_limma<0.1)
plot(log10(esa_cor[indg]), log10(pval_limma[indg]), xlab="Log p-value according to abnormality test",ylab="Log p-value according to Limma",main="Significant genes",col="red",lty="solid")
text(log10(esa_cor[indg]), log10(pval_limma[indg]),labels=genenames[indg],cex=0.75,pos=3);lines(c(0,-20),c(0,-20));
lines(c(0,-20),c(log10(0.000002),log10(0.000002)),col="red");lines(c(log10(0.000002),log10(0.000002)),c(0,-20),col="red");
dev.off()
#significant (bonferroni) hits by esa ordered by decreasing t-test significance
#POMP   EIF3E     C5orf24  MRLC2  UQCRH LDHB
#in alzheimer-contrl CRBN mental retardation

test with chisquare

for (sel in indg){
e=finalexpr[sel,];print(res_esa[[sel]]$summary)
best_chisq(finalpheno,e)
}
plot(1:length(res_esa[[sel]]$escurve),res_esa[[sel]]$escurve)


#Calibrating under random
res_esa_c=list();length(res_esa_c)<- nrow(finalexpr); 
finalpheno_rand=sample(finalpheno,length(finalpheno));
ptm= proc.time()
for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa_c[[i]]<- psea(finalpheno_rand,scale(e+rnorm(length(e),0,0.00000001)),w0,early=early,rep=20000,half=half,flatten=flatten);}
proc.time()-ptm
es_c=rep(0,nrow(finalexpr));esa_c=rep(1,nrow(finalexpr));for (i in 1:nrow(finalexpr)){es_c[i]=res_esa_c[[i]]$es;esa_c[i]=(res_esa_c[[i]])$pval;}
ind=which( esa_c>0.0000005)
res=lm(log(esa_c[ind])~ es_c[ind] + (es_c^2)[ind]+ (es_c^3)[ind]+ (es_c^4)[ind] + (es_c^5)[ind])
predict_pval=function(esx)exp(cbind(1,esx,esx^2,esx^3,esx^4,esx^5)%*%res$coefficients)#res$coefficients=c(-1.4680327,1.2494153,-0.6588248)

esa_cp=predict_pval(es_c)

plot(log10(predict_pval(es_c)),log10(esa_c));lines(c(0,-20),c(0,-20))#reconstruction

plot(log10(predict_pval(es_c)),es_c);lines(c(0,-20),c(0,-20))#reconstruction
plot(log10(predict_pval((1:1000)/100)),(1:1000)/100);lines(c(0,-20),c(0,-20))#reconstruction

