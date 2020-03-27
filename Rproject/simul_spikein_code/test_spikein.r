#!/hpf/tools/centos6/R/3.2.3/bin/Rscript
if (length(commandArgs(TRUE))<7){print("Error: incorrect number of arguments");if(NA)print("Error");}
dir=commandArgs(TRUE)[1];
r=as.numeric(commandArgs(TRUE)[2]);
testing_imbalance_ctrl=as.numeric(commandArgs(TRUE)[3]);
testing_imbalance_case=as.numeric(commandArgs(TRUE)[4]);
ndef=as.numeric(commandArgs(TRUE)[5]);
multiple=as.numeric(commandArgs(TRUE)[6]);
k_ind=as.numeric(commandArgs(TRUE)[7]);
testing_adding_ncase=0;if (length(commandArgs(TRUE))==8){testing_adding_ncase=as.numeric(commandArgs(TRUE)[8])}

local=FALSE;
if(local){
dir="/home/aziz/Desktop/aziz/testing_expression/";
testing_imbalance_ctrl=0;testing_imbalance_case=0;#At most one of these is True
r=0.1; ndef=100;#only used if testing imbalance
multiple=100
k_ind=1;
}

print(paste("Parameters:",r,testing_imbalance_ctrl,testing_imbalance_case));
library(limma)
library(car)#levene
source(paste(dir,"functions.r",sep=""))
ctrlexpr=read.table(paste(dir,"data/healthy_expr.txt",sep=""))

ns=c(200,300,450,600)
if (ndef==500)ns=c(500,650,800,1000)
if (ndef==100)ns=c(100,200,300,400)
ctrlmultiplier=c(0.1,0.25,0.5,0.8,1,1.2,2,4,10);#1 balanced
d=3; g=10; s=0.01;exact=TRUE;
nrep=10000000;accuracythresh=1/nrep;
w0=NULL;flatten=0.5;

if (testing_imbalance_case | testing_imbalance_ctrl){var=ctrlmultiplier} else var=ns;
allstats=array(0,dim=c(6,3,length(var),multiple));allstatsbonf=allstats;

ptmall= proc.time();
#Inverse gaussian transformation is optional for t-test, Levene and combined test
invtrans=function(x,doinverse=FALSE){if (doinverse)qnorm( (rank(x)- 3/8)/ (length(x) - 6/8 +1)   ); return(x)};


for (u in 1:length(var)){
  if (!testing_imbalance_case && !testing_imbalance_ctrl){n=ns[u];n0=ns[u];}#n:cases, n0:controls
  if (testing_imbalance_case){n0=ndef;n=round(ctrlmultiplier[u]*ndef);}
  if (testing_imbalance_ctrl){n=ndef;n0=round(ctrlmultiplier[u]*ndef);}
  if (testing_imbalance_case && testing_imbalance_ctrl){n=round(ctrlmultiplier[u]*ndef);n0=round(ctrlmultiplier[u]*ndef);}#For testing the addition of equal amounts of controls and neutral cases

  finalpheno=c(rep(1,n),rep(0,n0));
  indcase=which(finalpheno==1);
  indctrl=which(finalpheno==0);
  ptm= proc.time();
  calibration=calibrate_test(finalpheno,w0,rep=nrep,doall=TRUE,flatten=flatten);
  print(proc.time()-ptm)
  print(paste("Starting multiple trials for n=",n,"n0=",n0))
  for (k in 1:multiple){
    #Sample with noise from data and spike in
    finalexpr=t(apply(ctrlexpr,1,function(x)(sample(x,length(finalpheno),replace=TRUE)+rnorm(length(finalpheno),s=s))))
    truegenes=sample(1:nrow(finalexpr),g);
    for (i in truegenes){
    affected=sample(1:n,r*n);if (testing_adding_ncase)affected=sample(1:n,r*ndef);
    sign=runif(1)-0.5;
    finalexpr[i,affected]=finalexpr[i,affected]+ (sign/abs(sign))*d*sd(finalexpr[i,])
    }
    
    #t-test and wilcoxon
    pval_t=apply(finalexpr,1,function(e){e=invtrans(e);t.test(e[indcase],e[indctrl])$p.value});
    fdr_t=p.adjust(pval_t, method="fdr");
    pval_wilc=apply(finalexpr,1,function(e)wilcox.test(e[indcase],e[indctrl])$p.value);
    fdr_w=p.adjust(pval_wilc, method="fdr");

    #closed form es normalization ,no permutations 1min
    res_esa=list();length(res_esa)<- nrow(finalexpr);
    for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa[[i]]<- aziz.test(finalpheno,e,w0,rep=0,flatten=flatten);}
    res_ref=reformat_results(res_esa);es1=res_ref$es;

    #First approximation: Estimate from 1 set of permutations
    esa_m=get_calibrated_pvalues(calibration,es1);

    esa_p=esa_m;
    fdr_esa=p.adjust(esa_p, method="fdr");
    #limma
    design <- model.matrix(~ as.factor(finalpheno));#colnames(design) <- c("cases", "controls");
    fit <- lmFit(finalexpr, design=design)
    fit2 <- eBayes(fit)
    results2=topTable(fit2,coef=2,n=Inf,sort="none")

    #Levene test
    levres=apply(finalexpr,1,function(x){x=invtrans(x);leveneTest(x,group=as.factor(finalpheno))$"Pr(>F)"[1]})
    wf=-2*(log(pval_t)+log(levres));
    combres=1-pchisq(wf,df=4);
    fdr_comb=p.adjust(combres, method="fdr");fdr_lev=p.adjust(levres, method="fdr")

    bft=0.05/nrow(finalexpr);
    
    stats=function(x){if (length(x)==0)return(c(0,0,0));i1=length(which(x %in% truegenes));i2=length(which(!(x %in% truegenes)));
          return(c(i1/g,i2/length(x),i2/(nrow(finalexpr)-g)));};# Power, FDR, FPR
    #allstats[1,,u,k]=stats(which(fdr_esa<0.1));allstats[2,,u,k]=stats(which(results2[,5]<0.1));
    #allstats[3,,u,k]=stats(which(fdr_t<0.1));allstats[4,,u,k]=stats(which(fdr_lev<0.1));
    #allstats[5,,u,k]=stats(which(fdr_comb<0.1));allstats[6,,u,k]=stats(which(fdr_w<0.1))
    resind=list(which(fdr_esa<0.1),which(results2[,5]<0.1),which(fdr_t<0.1),which(fdr_lev<0.1),which(fdr_comb<0.1),which(fdr_w<0.1))
    resindbonf=list(which(esa_p<bft),which(results2[,4]<bft),which(pval_t<bft),which(levres<bft),which(combres<bft),which(pval_wilc<bft))
    allstats[,,u,k]=t(mapply(stats,resind));
    allstatsbonf[,,u,k]=t(mapply(stats,resindbonf));
  }
  print(apply(allstats[,,u,],c(1,2),mean));
}
print(proc.time()-ptmall)
resname=paste(dir,"spikein/spikein_",r,'-',n,'-',n0,"-",k_ind,sep="")

finalexpr=NULL;ctrlexpr=NULL; res_esa=NULL;calibration=NULL;fit=NULL;fit2=NULL; results2=NULL;#Free space
pval_t=NULL;pval_wilc=NULL;levres=NULL;fdr_t=NULL;fdr_lev=NULL;esa_p=NULL;esa_m=NULL;fdr_esa=NULL;res_ref=NULL;fdr_w=NULL;fdr_comb=NULL;combres=NULL;wf=NULL;#Free more space
save.image(paste(resname,".RData",sep=""))

library(ggplot2)
#plotdistr=data.frame(cond=as.factor(finalpheno),val=finalexpr[truegenes[1],])
#pdf("example025.pdf");ggplot(plotdistr, aes(x=val,fill=cond))+geom_density(alpha=0.3) ;dev.off()

if (length(var)>1){
selmeth=c(1,2,3,4,5,6)
col=c("black","blue","cyan","red","orange","yellow")
meth=c("Aberration test","Limma","t-test","Levene test","Combined test","Wilcoxon")
mat=apply(allstats[selmeth,1,,],c(1,2),mean);mat2=apply(allstats[selmeth,2,,],c(1,2),mean);

pdf(paste(resname,".pdf",sep=""));
plot(as.factor(var),c(rep(0,length(var)-1),1),pch=NA,ylab="Power (solid) / FDR (Dashed)",border="white");
for (i in 1:length(selmeth))print(lines(as.factor(var),mat[i,],col=col[i],lty="solid"));
for (i in 1:length(selmeth))print(lines(as.factor(var),mat2[i,],col=col[i],lty="dashed"));
legend("topleft",meth,fill=col)
dev.off()
}


