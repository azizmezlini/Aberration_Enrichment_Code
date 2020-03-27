load(paste(dir0,fname2,".RData",sep="")) #ncol(factors) is used to only adjust for genes p-vals
dir=dir0; print(k==k0); #Resetting the new directory and verifying we are using the same k
#nrep=10000000; exact=FALSE;accuracythresh=0.05/nrow(finalexpr)
source(paste(dir,"Rproject/functions.r",sep=""))

print("Running differential expression tests...")

#t-test and Limma

pval_wilc=apply(finalexpr,1,function(e)wilcox.test(e[which(finalpheno==1)],e[which(finalpheno==0)])$p.value);
pval_t=apply(finalexpr,1,function(e)t.test(e[which(finalpheno==1)],e[which(finalpheno==0)])$p.value);
library(limma)
design <- model.matrix(~ as.factor(finalpheno));#colnames(design) <- c("cases", "controls");
fit <- lmFit(finalexpr, design=design)
fit2 <- eBayes(fit)
results2=topTable(fit2,coef=2,n=Inf,sort="none")
pval_limma=results2[,4];
fdr_limma=results2[,5];
fdr_t=p.adjust(pval_t, method="fdr");fdr_w=p.adjust(pval_wilc, method="fdr");
print("Computing enrichment scores...");

w0=NULL;
#closed form es normalization ,no permutations 1min
ptm= proc.time()
res_esa=list();length(res_esa)<- nrow(finalexpr);
for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa[[i]]<- aziz.test(finalpheno,e,w0,rep=0,ignoremax=0);}
res_ref=reformat_results(res_esa);
print(proc.time()-ptm)

#First approximation: Estimate from 1 set of permutations
print("Computing first pvalues approximation using a single set of permutations for all genes...");ptm= proc.time()
ptm= proc.time()
calibration=calibrate_test(finalpheno,w0,rep=nrep,doall=TRUE,ignoremax=0);
esa_m=get_calibrated_pvalues(calibration,res_ref$es);
print(proc.time()-ptm)

#Exact computation : Full permutation on every gene 2h
esa_p=esa_m;
if(exact){
  print("Computing exact p-values with permutations for every gene...")
  res_esa_p=list();length(res_esa_p)<- nrow(finalexpr); 
  ptm= proc.time()
  for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa_p[[i]]<- aziz.test(finalpheno,e,w0,rep=nrep);}
  print(proc.time()-ptm)
  es_p=rep(0,nrow(finalexpr));esa_p=rep(1,nrow(finalexpr));for (i in 1:nrow(finalexpr)){es_p[i]=res_esa_p[[i]]$es;esa_p[i]=(res_esa_p[[i]])$pval;}
}

#modeling p-values to compute good approximations
esa_cor=correct_zero_pvalues(calibration,log((nrep:1)/nrep),res_ref$es,esa_p,accuracythresh,maxpvaltrain=0.05)
if (exact)esa_cor=correct_zero_pvalues(res_ref$es,log(esa_p),res_ref$es,esa_p,accuracythresh,maxpvaltrain=0.05)
#plot(log10(predict_pval(res_ref$es)),log10(esa_p));lines(c(0,-20),c(0,-20))#training

fdr_esa=p.adjust(esa_cor, method="fdr");

#recalculate FDR so we only take gene pvals and not factors pvals into account.
rg=1:(length(genenames)-ncol(factors)); #rg : real genes
fdr_esa[rg]=p.adjust(esa_cor[rg], method="fdr");
fdr_limma[rg]=p.adjust(pval_limma[rg], method="fdr");fdr_t[rg]=p.adjust(pval_t[rg], method="fdr");fdr_w[rg]=p.adjust(pval_wilc[rg], method="fdr");

#names=genenames;rg=1:(length(genenames))# for outdated versions
show_res=data.frame(genenames,names,fdr_esa,fdr_limma,res_ref,pval2=esa_p,approx=esa_cor,limma=pval_limma,fdr_w,wilc=pval_wilc);
#show_res=data.frame(genenames,names,fdr_esa,fdr_w,res_ref,pval2=esa_p,approx=esa_cor,wilc=pval_wilc);

signifthresh=0.05/length(genenames)
print(format( show_res[esa_p<signifthresh,],digits=3));
print(format( show_res[fdr_esa<0.1 & res_ref$r<0.3,],digits=3));
#indplotexpr=which(fdr_esa<0.1 & res_ref$r<0.3);for (i in indplotexpr){gnamei=genenames[i];if(gnamei=="")gnamei=names[i];  pdf(paste(gnamei,"_gg.pdf",sep=""));expr_plot_gg(finalexpr[i,],finalpheno,gnamei,res_ref[i,]);dev.off();}
print("Number of genes found with FDR<0.1");printvenn(which(fdr_esa[rg]<0.1),which(fdr_limma[rg]<0.1));
print("Number of genes found (Bonferroni)");printvenn(which(esa_cor[rg]<signifthresh),which(pval_limma[rg]<signifthresh));
print("With wilcoxon: fdr/Bonf:");printvenn(which(fdr_esa[rg]<0.1),which(fdr_w[rg]<0.1));printvenn(which(esa_cor[rg]<signifthresh),which(pval_wilc[rg]<signifthresh));

write.table(show_res,paste(dir,fname,".txt",sep=""),row.names=FALSE,sep="\t")

expr=NULL;exprscaled=NULL;residuals = NULL;factors = NULL;weights =NULL;#Freeing up space
save.image(paste(dir,fname,".RData2",sep=""))


plot_graph=function(indg, cex=0.75,pos=3,pch=NULL,lcex=1){
  plot(-log10(pval_limma[indg]), -log10(esa_cor[indg]), xlab="Limma Negative Log p-value",ylab="Aberration Test Negative Log p-value",main="Significant genes",col=col[res_ref$es_sign[indg]],
       xlim=minrange(c(-0.55,3)+range(-log10(esa_cor[indg]))) )
  text(-log10(pval_limma[indg]), -log10(esa_cor[indg]),labels=genenames[indg],cex=cex,pos=pos,pch=pch,col=col[res_ref$es_sign[indg]]);lines(c(-1,20),c(-1,20));
  lines(c(-1,20),c(-log10(signifthresh),-log10(signifthresh)),col="red");lines(c(-log10(signifthresh),-log10(signifthresh)),c(-1,20),col="red");
  legend("bottomright",c("Under-expressed","Over-expressed"),fill=col,cex=lcex)
}

col=c("blue","green");
indg=which(fdr_esa[rg]<0.1 | fdr_limma[rg]<0.1);pdf(paste(dir,fname,"_fdr01.pdf",sep=""));if(length(indg))plot_graph(indg);dev.off();
indg=which(fdr_esa[rg]<0.2 | fdr_limma[rg]<0.2);pdf(paste(dir,fname,"_fdr02.pdf",sep=""));if(length(indg))plot_graph(indg);dev.off();

library(qqman)
coord=read.table(paste(dir,"data/gene_coordinates.txt",sep=""));
manhattan_plot(coord[,c(5,2,3)],genenames[rg],data.frame(Aberration_test=esa_cor[rg],Limma=pval_limma[rg]),paste(dir,fname,"_manhattan.pdf",sep=""),annotatePval=0.001)
pdf(paste(dir,fname,"_qq.pdf",sep=""));qq(esa_cor[rg]);dev.off(); #esa_cor instead of esa_p

if (perm){
  finalpheno2=sample(finalpheno,length(finalpheno))
  ptm= proc.time()
  res_esa2=list();length(res_esa2)<- nrow(finalexpr);
  for (i in 1:nrow(finalexpr)){e=finalexpr[i,];res_esa2[[i]]<- aziz.test(finalpheno2,e,w0,rep=0);}
  print(proc.time()-ptm)
  res_ref2=reformat_results(res_esa2);
  esa_m2=get_calibrated_pvalues(calibration,res_ref2$es);
  esa_cor2=correct_zero_pvalues(calibration,log((nrep:1)/nrep),res_ref2$es,esa_m2,accuracythresh,maxpvaltrain=0.05)
  fdr_esa2=p.adjust(esa_cor2, method="fdr");fdr_esa2[rg]=p.adjust(esa_cor2[rg], method="fdr");
}


compositeplot=TRUE; candidates=which(fdr_esa[rg]<0.1); #All significant genes; One example gene (gsname); qq plot; qq plot of permuted
if (perm & compositeplot &length(candidates)){
  order_candidates=order(res_ref$es_indmax[candidates]);
  gsind=candidates[order_candidates][1];gsname=genenames[gsind];if (is.na(gsname)| gsname=="")gsname=names[gsind];#"CRBN"#"UQCRH"
  pdf(paste(dir,fname,"_res_all_",gsname,".pdf",sep=""));
  par(mfrow = c(2,2))
  col=c("blue","green");  indg=which(fdr_esa[rg]<0.1 | fdr_limma[rg]<0.1);
  if(length(indg))p1=plot_graph(indg,cex=0.6,pos=3,pch=8,lcex=0.75);
  #plot(-log10(pval_limma[indg]), -log10(esa_cor[indg]), xlab="Limma Negative Log p-value",ylab="Aberration Test Negative Log p-value",main="Significant genes",col=col[res_ref$es_sign[indg]],
  #     xlim=minrange(c(-0.55,3)+range(-log10(esa_cor[indg]))) )
  #text(-log10(pval_limma[indg]), -log10(esa_cor[indg]),labels=genenames[indg],cex=0.6,pch=8,pos=3,col=col[res_ref$es_sign[indg]]);lines(c(-1,20),c(-1,20));
  #lines(c(-1,20),c(-log10(signifthresh),-log10(signifthresh)),col="red");lines(c(-log10(signifthresh),-log10(signifthresh)),c(-1,20),col="red");
  #legend("topleft",c("Under-expressed","Over-expressed"),fill=col,cex=0.75)
  e=finalexpr[gsind,];p2=expr_plot2(e,finalpheno,gsname);
  #col=c("black","red");e=finalexpr[gsind,];
  #plot(density(e[which(finalpheno==0)]),col=col[1],main=paste(gsname,"expression levels"));lines(density(e[which(finalpheno==1)]),col=col[2]);
  #legend("topleft",c("controls","cases"),fill=col)
  p3=qq(esa_p[rg],ylim=c(0,7.5));p4=qq(esa_m2[rg],ylim=c(0,7.5));#Maybe I should use esa_cor and esa_cor2
  dev.off()
}

#sort -g -k4 CD/CD30.txt | awk '{if (NR>1)print substr($1,2,length($1)-2)}' - | head
#awk '{t=0.000002;if (NR>1 && ($6<0.1 && $7>0.1 || $4<t && $5>t)&& $4/($5+0.0000000001*t) <0.1 )print substr($1,2,length($1)-2)}' Alzheimer_ctrl/Alzheimer_ctrl30.txt > Alzheimer_ctrl/novel_genes.txt
x=esa_cor/pval_limma
y=order(fdr_esa)[1:length(which(fdr_esa<0.1))]
genenames[y[ which(x[y]<0.1 & pval_limma[y]>signifthresh) ]]#genenames[y[ which(x[y]<0.1 & fdr_limma[y]>0.1) ]]

#Plot expression
#for (i in indg){print("");expr_plot(finalexpr[i,],finalpheno,paste(dir,name,"/",name,k,"_",genenames[i],sep=""))}

#Correlation analysis
#if (length(indg)>1)plot_2g(finalexpr,res_esa, indg[1],indg[2],paste(dir,fname,sep=""),genenames)

#PCA plot (found useless)
#respca=prcomp(finalexpr,center = TRUE,scale. = TRUE)
#pdf(paste("eigen_",name,".pdf",sep=""));plot(1:length(respca$sdev),respca$sdev);dev.off()

exactall=F;#verfify the approx pvalues are the same as the full perm ones
if(exactall){
  indrep=which(show_res$fdr_esa<0.1);indrep2=sample(indrep,50)
  print("Computing exact p-values with permutations for every gene...")
  res_esa_p=list();length(res_esa_p)<- nrow(finalexpr); 
  ptm= proc.time()
  for (i in indrep2){e=finalexpr[i,];res_esa_p[[i]]<- aziz.test(finalpheno,e,w0,rep=nrep,doall = T);}
  print(proc.time()-ptm)
  es_p=rep(0,nrow(finalexpr));esa_p=rep(1,nrow(finalexpr));for (i in indrep2){es_p[i]=res_esa_p[[i]][[1]]$es;esa_p[i]=(res_esa_p[[i]][[1]])$pval;}

pdf("pvalues_match.pdf");
plot(-log10(esa_p[indrep2]),-log10(esa_m[indrep2]),xlab="Full permutation p-values",ylab="Approximative p-values",xlim=c(3,8),ylim=c(3,8))
lines(c(3,8),c(3,8))
dev.off()

ggplot(data.frame(res_esa[[2344]]$escurve))+ geom_step(aes(x=seq_along(res_esa[[2344]]$escurve), y=res_esa[[2344]]$escurve)) +
  xlab("Individuals ordered by decreasing gene expression level") +ylab("Normalized enrichment score") + theme_bw()#CRBN curve plot
}

