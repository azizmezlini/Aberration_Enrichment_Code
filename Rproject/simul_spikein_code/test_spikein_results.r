

dir="/hpf/largeprojects/agoldenb/aziz/testing_expression/";#dir="/home/aziz/Desktop/aziz/testing_expression/";
r=0.05;n=1000; n0=1000;
#r=0.1;n=400;n0=400
#r=0.1;n=5000;n0=500
#r=0.05;n=750;n0=7500
#r=0.05;n0=750;n=7500;
njobs=10;
resname0=paste(dir,"spikein/spikein_",r,'-',n,'-',n0,sep="");

resname=paste(resname0,"-",1,sep="");
load(paste(resname,".RData",sep=""));
allstatsk=array(0,dim=c(1,1,1,njobs)*dim(allstats));allstatsbonfk=allstatsk;
allstatsk[,,,1:multiple]=allstats;allstatsbonfk[,,,1:multiple]=allstatsbonf;

if (njobs>1){
for (k_ind in 2:njobs){
resname=paste(resname0,"-",k_ind,sep="");
load(paste(resname,".RData",sep=""));
allstatsk[,,,(1:multiple)+((k_ind-1)*multiple)]=allstats;
allstatsbonfk[,,,(1:multiple)+((k_ind-1)*multiple)]=allstatsbonf;
}
}

for (u in 1:length(var))print(apply(allstatsk[,,u,],c(1,2),mean));
for (u in 1:length(var))print(apply(allstatsbonfk[,,u,],c(1,2),mean));

t1="# of cases";
t0="# of controls";
if (n0>n){axe=t0;title=paste(t1,"=",n,". Varying ",t0)};
if (n>n0){axe=t1;title=paste(t0,"=",n0,". Varying ",t1)};
if(n0==n){axe="Sample size"; title="Power analysis";}

library(ggplot2)
if (length(var)>1){
  selmeth=c(1,2,3,4,5,6)
  col=c("red","blue","cyan","orange","black","yellow")
  meth=c("Aberration test","Limma","t-test","Levene test","Combined test","Wilcoxon")
  mat=apply(allstatsk[selmeth,1,,],c(1,2),mean);mat2=apply(allstatsk[selmeth,2,,],c(1,2),mean);
  
  pdf(paste(resname0,".pdf",sep=""));
  plot(as.factor(var),c(rep(0,length(var)-1),1),pch=NA,xlab=axe,ylab="Power (solid) / FDR (Dashed)",border="white");
  for (i in 1:length(selmeth))print(lines(as.factor(var),mat[i,],col=col[i],lty="solid"));
  for (i in 1:length(selmeth))print(lines(as.factor(var),mat2[i,],col=col[i],lty="dashed"));
  legend("topleft",meth,fill=col)
  title(title)
  dev.off()
}



