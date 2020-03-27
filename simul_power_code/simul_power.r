dir="/home/aziz/Desktop/aziz/aziz_test/";
source(paste(dir,"functions.r",sep=""))
library(car)
multiple=200
dirname="robustness_s/";
ns=c(100,200,300,400,500,600,800,1200,2000)
rs=c(0.02,0.035,0.05,0.075,0.1,0.15,0.2)
ds=c(1.5,2,2.5,3,4)
#dirname="robustness_precise/";
#ns=c(800,1200,1800,2700,4000,6000,9000,14000)
#rs=c(0.005,0.01,0.02,0.03)
#ds=c(1.5,2,2.5,3,4)
#dirname="robustness_broad2/";
#ns=c(20,30,50,75,100,150,200,300)
#rs=c(0.2,0.3,0.5,0.75,1)
#ds=c(1,1.5,2,2.5,3)
#dirname="robustness_low_d/";
#ns=c(200,300,500,750,1000,1500,2500,5000)
#rs=c(0.025,0.05,0.075,0.1,0.2)
#ds=c(0.15, 0.3, 0.5, 0.7, 1.)
#dirname="robustness_low_d2/";
#ns=c(50,100,200,300,500,750,1000)
#rs=c(0.1,0.2,0.3,0.5,0.75,1)
#ds=c(0.15, 0.3, 0.5, 0.7, 1.)

vs=c(0.2)#vs=c(0.2,0.3,0.5,1)
res=array(0,dim=c(length(ns),length(rs),length(ds),length(vs),multiple,8))
w0=NULL;nrep=10000000;accuracythresh=1/nrep;

resdirname=paste(dir,dirname,sep="");dir.create(resdirname);

for (i in 1:length(ns)){
n=ns[i];pheno=c(rep(1,n),rep(0,n));
#create an e matrix for each ns and test all together to gain time
expr=matrix(0,multiple*length(vs)*length(ds)*length(rs),2*n);


for (j in 1:length(rs)){
for (u in 1:length(ds)){
for (v in 1:length(vs)){
for (k in 1:multiple){

r=rs[j];m=100;
s=vs[v]*m;
dis=ds[u]*s;

e=rnorm(n*2,m,s);#cases+controls 
# some cases are shifted
ind=sample(1:n,floor(r*n)); 
sign=sample(c(-1,1),1);
if(length(ind))e[ind]<- e[ind]+sign* dis;
#pheno=c(rep(1,n),rep(0,n));expr_plot(e,pheno,paste("example_",r,"_",n,".pdf",sep=""));pheno2=as.factor(c(rep("cases",n),rep("controls",n)))
#p1 <- ggplot(data.frame(expression=e,phenotype=pheno2), aes(x=expression, fill=phenotype)) + geom_density(alpha=.3);
#p2<- ggplot(data.frame(expression=e,phenotype=pheno2), aes(x=expression, fill=phenotype,xmin=-6)) + geom_density(alpha=.3);
#pdf(paste("nicerexample_",r,"_",n,"2.pdf",sep=""));
#grid.arrange(p1, p2,ncol=1)
#dev.off()
expr[k+(v-1)*multiple+(u-1)*multiple*length(vs)+(j-1)*multiple*length(vs)*length(ds),]=e;

}}}}

pval_t=apply(expr,1,function(x)t.test(x[1:n],x[(n+1):(2*n)])$p.value);
pval_lev=apply(expr,1,function(x)leveneTest(x,group=as.factor(pheno))$"Pr(>F)"[1])#Levene test
pval_ks=apply(expr,1,function(x)ks.test(x[1:n],x[(n+1):(2*n)])$p.value);
pvals_kw=apply(expr,1,function(x)kruskal.test(x,as.factor(pheno))$p.value)
pval_wilc=apply(expr,1,function(x)wilcox.test(x[1:n],x[(n+1):(2*n)])$p.value);

#closed form es normalization ,no permutations 1min
res_esa=list();length(res_esa)<- nrow(expr);
for (p in 1:nrow(expr)){x=expr[p,];res_esa[[p]]<- aziz.test(pheno,x,w0,rep=0);}
res_ref=reformat_results(res_esa);es1=res_ref$es;


#e=rnorm(2*n);res_m <- aziz.test(pheno,e,w0,rep=nrep,doall=TRUE);
#sortperm=sort(res_m$perm); esa_m=1- findInterval(es1,sortperm,left.open=TRUE)/nrep; 
calibration=calibrate_test(pheno,w0,rep=nrep,doall=TRUE);
esa_m=get_calibrated_pvalues(calibration,es1);
#correct if 0
esa_cor=correct_zero_pvalues(calibration,log((nrep:1)/nrep),es1,esa_m,accuracythresh)


#print(paste(res_esa$pval,res_esa$summary))
#plot(1:length(res_esa$escurve),res_esa$escurve)
#best_chisq(pheno,expr)

#testbinom=sapply(1:n, function(x)dbinom(sum(pheno[ord[1:x]]),x,0.5));
#hyper=testhypergeom(pheno,expr)
#mindan=mindaniel(pheno,expr)
#meandan=meandaniel(pheno,expr)

pvals=list(pval_t,pval_lev,pval_ks,pvals_kw,pval_wilc,esa_m,esa_cor,es1)
reshape_revdepth=function(x)aperm(array(x,dim=c(multiple,length(vs),length(ds),length(rs))));
print(i);
for (p in 1:length(pvals))res[i,,,,,p]=reshape_revdepth(pvals[[p]]);#filling the array in reverse depth
#res[i,j,u,v,k,]=c(pval_t,pval_lev,pval_ks,pvals_kw,pval_wilc,res_esa$pval,res_esa$es)
}

#resold=res;deg=which(dim(res)==1);library(abind);#Degenerate dimension quick fix; Useful if I am not changing vs
#res=abind(resold,resold,along=deg)

vsel=1;#results for a specific value of vs
resv=res[,,,vsel,,];

library(lattice)
res2=round(apply(resv<=0.05,c(1,2,3,5),sum)/(multiple))
res2hard=round(apply(resv<=0.05/25000,c(1,2,3,5),sum)/(multiple))
res3t=apply(log10(resv[,,,,6]+accuracythresh)-log10(resv[,,,,1]+accuracythresh),c(1,2,3),mean)
res3lev=apply(log10(resv[,,,,6]+accuracythresh)-log10(resv[,,,,2]+accuracythresh),c(1,2,3),mean)
res3wil=apply(log10(resv[,,,,6]+accuracythresh)-log10(resv[,,,,5]+accuracythresh),c(1,2,3),mean)
chosencol=function(x){tcol=x[,,6]+2*x[,,2]+4*x[,,1]+8*x[,,5];rownames(tcol)<- ns;colnames(tcol)<- rs; return(tcol)}
chosencol3=function(x){tcol=x;tcol[tcol>3.99] <- 3.99;tcol[tcol< -3.99] <- -3.99;rownames(tcol)<- ns;colnames(tcol)<- rs; return(tcol)}

for (dsel in 1:length(ds)){
pres2=res2[,,dsel,];pres2hard=res2hard[,,dsel,];
color=chosencol(pres2);colorhard=chosencol(pres2hard);
logdiffs=list(-res3t[,,dsel],-res3lev[,,dsel],-res3wil[,,dsel])
leastlogdiff=pmin(logdiffs[[1]],logdiffs[[2]],logdiffs[[3]])
color3=lapply(c(list(leastlogdiff),logdiffs),chosencol3)
suffix=c("All","t","l","w")

pdf(paste(resdirname,"soft_",dsel,"_",vsel,".pdf",sep=""))
print(levelplot(color,xlab="N",ylab="R",col.regions=c("black",rainbow(14),"white"), at=seq(0,16,0.99)))
dev.off()
pdf(paste(resdirname,"hard_",dsel,"_",vsel,".pdf",sep=""))
print(levelplot(colorhard,xlab="N",ylab="R",col.regions=c("black",rainbow(14),"white"), at=seq(0,16,0.99)))
dev.off()
for (ldiff in 1:length(color3)){
pdf(paste(resdirname,"log10_",suffix[ldiff],"_",dsel,"_",vsel,".pdf",sep=""))
print(levelplot(color3[[ldiff]],xlab="N",ylab="R",col.regions=colorRampPalette(c("darkred","white","darkblue"))(16), at=seq(-4,4,0.4999999)))
dev.off()
}
}

save.image(paste(resdirname,"allres.RData",sep=""))




#Calibrating ES-max  score and pvalue
ns=c(100,200,300,400,500,600,800,1200)
r=0;
m=100;
s=25;
multiple=1000000
res_c=array(0,dim=c(multiple,length(ns),2))
for (i in 1:length(ns)){
n=ns[i];
for (k in 1:multiple){
pheno=c(rep(1,n),rep(0,n))
expr=rnorm(n*2,m,s);#cases+controls 
res_esa=aziz.test(pheno,expr,w0,rep=2000000)
res_c[k,i,]=c(res_esa$pval,res_esa$es)
}
}

#Separate regression
coef=matrix(0,length(ns),6)
for (i in 1:length(ns)){
ind=which( esa_c>0.0000005)
res_lm=lm(log(esa_c[ind])~ res_c[ind,i,2] + (res_c^2)[ind,i,2]+ (res_c^3)[ind,i,2]+ (res_c^4)[ind,i,2] + (res_c^5)[ind,i,2])
coef[i,]=res_lm$coefficients
}
predict_pval=function(esx)exp(cbind(1,esx,esx^2,esx^3,esx^4,esx^5)%*%coef[4,])
i=8;plot(log10(predict_pval(res_c[,i,2])),log10(res_c[,i,1]));lines(c(0,-20),c(0,-20))#reconstruction

#Together
es_c=as.numeric(res_c[,,2]);esa_c=as.numeric(res_c[,,1]);
lognmat=matrix(log(ns),length(ns),multiple);logn=as.numeric(t(lognmat));
ind=which( esa_c>0.0000005)
res_lm=lm(log(esa_c[ind])~ logn[ind]+ es_c[ind]+ (es_c*logn)[ind] + (es_c^2)[ind]+ (es_c^3)[ind]+ (es_c^4)[ind] + (es_c^5)[ind])

predict_pval_all=function(esx,logn=5)exp(cbind(1,logn,esx,esx*logn,esx^2,esx^3,esx^4,esx^5)%*%res_lm$coefficients)

i=4;plot( log10( predict_pval_all(res_c[,i,2],logn=log(ns[i]) ) ),log10(res_c[,i,1]));lines(c(0,-20),c(0,-20))#reconstruction



#esa_cp=predict_pval(es_c)
esa_c=res_c[,i,1];es_c=res_c[,i,2];
plot(log10(predict_pval(es_c)),log10(esa_c));lines(c(0,-20),c(0,-20))#reconstruction

plot(log10(predict_pval(es_c)),es_c);#reconstruction


