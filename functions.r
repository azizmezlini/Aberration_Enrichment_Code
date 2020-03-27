
reformat_results=function(res_esa){# From results listed by gene first to results listed by attribute first
  es1=rep(0,length(res_esa));es_sign=rep(-1,length(res_esa));es_indmax=rep(-1,length(res_esa));
  r=rep(1,length(res_esa));or=r;ncas=r-1;pval=r;#initializing
  for (i in 1:length(res_esa)){
    es1[i]=(res_esa[[i]])$es;
    es_sign[i]=(res_esa[[i]])$direction;
    es_indmax[i]=(res_esa[[i]])$esind[es_sign[i]];
    r[i]=(res_esa[[i]])$oddcas;or[i]=(res_esa[[i]])$oddratio;
    ncas[i]=(res_esa[[i]])$ncas;pval[i]=(res_esa[[i]])$pval;
  }#WARNING the direction attribute was not previously returned
  return(data.frame(es=es1,es_sign=es_sign,es_indmax=es_indmax,r=r,or=or,ncases=ncas,pval=pval))  
}

calibrate_test=function(pheno,w0,rep=10000000,doall=TRUE,unidirectional=0,flatten=0.5,ignoremax=0){
  e=rnorm(length(pheno));res_m <- aziz.test(pheno,e,w0,rep=rep,doall=doall,unidirectional=unidirectional,flatten=flatten,ignoremax=ignoremax);
  return(sort(res_m$perm)); 
}

get_calibrated_pvalues=function(calibration,es1){
  return(1- findInterval(es1,calibration)/length(calibration)); 
  #return(1- findInterval(es1,calibration,left.open=TRUE)/length(calibration)); #Might be better but does not work on hpf
}

correct_zero_pvalues=function(trainx,trainy,es1,esa_m,accuracythresh,maxpvaltrain=0.05){#calibration,log((length(calibration):1)/length(calibration))
  regdata=data.frame(x=trainx,y=trainy);
  ind=which(regdata[,2]>= log(accuracythresh) & regdata[,2]< log(maxpvaltrain));
  res_lm=lm(y~ poly(x,3),data=regdata[ind,]);
  predict_pval=function(esx)exp(predict(res_lm,data.frame(x=esx)));
  esa=predict_pval(es1);
  esa_cor=esa_m;
  toreplace=which(esa_m<accuracythresh);#Was <accuracythresh instead of ==0
  if (length(toreplace))esa_cor[toreplace]=sapply(esa[toreplace],function(x)min(x,accuracythresh))
  #pdf(paste(dir,name,"/",name,k,"_pvalues_predicted_versus_estimated.pdf",sep=""));
  #plot(-log10(esa_p),-log10(predict_pval(res_ref$es)), xlab="Negative Log10 Estimated p-values",ylab="Negative Log10 predicted p-values",main="Predicted p-values quality");lines(c(0,20),c(0,20));dev.off();
  #pdf(paste(dir,name,"/",name,k,"_pvalues_function_of_scores.pdf",sep=""));
  #plot((1:1000)/100,-log10(predict_pval((1:1000)/100)), xlab="Max Enrichment Score ",ylab="Negative Log10 p-values",main="Predicted p-values function");dev.off();
  return(esa_cor);
}


aziz.test=function(pheno,pred,w=NULL,rep=100000,doall=FALSE,eps=0.000000001,unidirectional=0,flatten=0.5,ignoremax=0,approx=1,weighting=TRUE,trimmed=0.75){# negative pred are used as weights by default
 prep=preprocess_weights(pheno,pred,w=w,eps=eps,weighting=weighting,flatten=flatten)
 ord=prep$ord;signmult=prep$signmult;w=prep$w;
 
 esnorm=generate_es_normalizer_generic(w[ord],pheno,approx,weighting=weighting);
 esnormrev=generate_es_normalizer_generic(w[rev(ord)],pheno,approx,weighting=weighting);
 esa0=function(s,o,norm)estransformer(cumsum( (s*w)[o]),norm,trimmed=trimmed); #Enrichment score is computed and standardized here.
 esa_pos=function(s)(as.numeric(esa0(s,ord,esnorm)))
 esa_neg=function(s)(as.numeric(esa0(s,rev(ord),esnormrev)))
 esa_pos_m=function(s)max(esa_pos(s));esa_neg_m=function(s)max(esa_neg(s));
 if (ignoremax>0){esa_pos_m=function(s)max(esa_pos(s)[-(1:ignoremax)]);esa_neg_m=function(s)max(esa_neg(s)[-(1:ignoremax)]);}
 esa_2d=function(s)max(c(esa_pos_m(s),esa_neg_m(s)))
 #if (ignoremax>0)esa_2d=function(s)max(c(esa_pos(s)[-(1:ignoremax)],esa_neg(s)[-(1:ignoremax)])); #treating both directions
 esa_applied=esa_2d; if (unidirectional==1)esa_applied=esa_pos_m;if (unidirectional==-1)esa_applied=esa_neg_m;#both direction, positive direction (underexpression) or negative direction
 esperm=permute(rep,signmult,esa_applied,doall);
 dir1=esa_pos(signmult);dir2=esa_neg(signmult);
 esm=c(max(dir1),max(dir2)); esind=c(which.max(dir1),which.max(dir2));
 direction=1; escurve=dir1;if(esm[2]>esm[1]){direction=2;escurve=dir2;};# WARNING The ignoremax parameter could cause discrepancies here
 vip=1:esind[direction];if (direction==2)vip=(length(pheno)+1-esind[direction]):length(pheno);
 ncas=sum(pheno[ord][vip]);nctr=esind[direction]-ncas; 
 ncastotal=sum(pheno);nctrtotal=length(pheno)-ncastotal;
 oddcas=ncas/ncastotal;oddctr=nctr/nctrtotal;oddratio=oddcas/oddctr;
 #summary=paste("Score:",format(esperm$real,digits=3),"direction:",direction,"r:",format(oddcas,digits=3),"pval:",esperm$pval)
 #details=paste("indmax:",esind[direction],"cases:",ncas,"odd ratio:",format( oddratio ,digits=3))
 return(list(es=esperm$real,esm=esm,esind=esind,pval=esperm$pval,direction=direction,escurve=escurve,vip=vip,perm=esperm$perm,ncas=ncas,oddratio=oddratio,oddcas=oddcas))
}

preprocess_weights=function(pheno,pred,w=NULL,eps=0.000000001,weighting=TRUE,flatten=0){
  if (eps>0)pred=scale(pred+rnorm(length(pred),0,eps));#Add noise to avoid equal values biasing the ranking
  if (length(w)==0){w=abs(pred);if(flatten)w[which(w<flatten)]=flatten;}
  if (length(w)==1)w=rep(1,length(pred));
  ord=order(pred,decreasing=FALSE);
  ind1=which(pheno==1);ind0=which(pheno==0);
  signmult=pheno;signmult[ind0]=-1;
  if(weighting){signmult[ind1]=1/length(ind1);signmult[ind0]=-1/length(ind0);}
  return(list(ord=ord,w=w,signmult=signmult))
}

estransformer=function(x0,esnorm,trimmed=0.75){#only return the 3/4 by default, 
x=(x0/esnorm$sd)#Was x=((x0-esnorm$mean)/esnorm$sd)
x=x[1:round(trimmed*length(x))];
return(x)
}

permute=function(rep, phen,f,doall=FALSE){#rep assumed >100 or 0
 real=f(phen)
 if (rep){
   lrep=100;perm=sapply(1:lrep,function(x){f(sample(phen,length(phen)))});
   pval=sum(perm>=real)/length(perm);
   lrepall=c(1000,10000,100000,1000000,10000000); lrepall=lrepall[which(lrepall<rep)]
   for (lrep in lrepall){
    if (!doall & pval<(100/lrep)){#pval < 0.1 (100/1000)  -> Do lrep=1000
     perm=sapply(1:lrep,function(x){f(sample(phen,length(phen)))});
     pval=sum(perm>=real)/length(perm);
    }
   }
   if ((pval< 100/rep & rep>100)| doall){
    perm=sapply(1:rep,function(x){f(sample(phen,length(phen)))})
    pval=sum(perm>=real)/length(perm);
   }
 }else {pval=1;perm=c();}
 #for (i in 1:rep){ permphen=sample(phen,length(phen)); perm[i]=f(permphen); } 
 return(list(real=real,perm=perm,pval=pval))
}

generate_es_normalizer_generic=function(w,pheno,approx=1,weighting=TRUE){ #weighting corrects for case-control imbalance within the steps
if (approx==1 & weighting) return(generate_es_normalizer(w,pheno))
if (approx==1) return(generate_es_normalizer_old(w,pheno))
return(generate_es_normalizer2(w,pheno))
}

generate_es_normalizer=function(w,pheno){ #w is the ranked expression
p=length(which(pheno==1))/length(pheno)
#f(x)=(1/q1)(sum(positivecases%*%w_positive)-(1/q2)sum(controls%*%w_controls))=(1/q1 + 1/q2)(sum(positivecases%*%w_positive)-(1/q2)sum(everything%*%w)
#mean=mean(w)*(np/q1+ np/q2 - n/q2)= mean(w)*(n/N + nq1/Nq2 - n/q2)= 0;sd=replace the 2 from prevous sd by (1/q1 + 1/q2)
#ES score is 
cumsumw=cumsum(w)
cumsumwsq=cumsum(w^2);
m=0
s=( (1/length(which(pheno==1))) + (1/length(which(pheno==0))) )*sqrt(p*(1-p)/(length(pheno)-1))*sapply(1:length(w),function(n)sqrt( length(pheno)*cumsumwsq[n] - (cumsumw[n])^2 ))
#Can also be written   n (length(pheno)*cummeanwsq[n] - n (cummeanw[n])^2 ) =  n (length(pheno)*var(w)+ cummeanw[n])^2  *(length(pheno)-n) )
#So the only new term is the addition of  p*(1-p) n (length(pheno)*var(w)) to the sqrt
return(list(mean=m,sd=s))
}

generate_es_normalizer_old=function(w,pheno){ #w is the ranked expression
p=length(which(pheno==1))/length(pheno)
#x from hypergeometric. Corresponding es score is f(x)=(sum(positivecases%*%w_positive)-sum(controls%*%w_controls));  We know that sum(everyone%*%w)=sum(positivecases%*%w_positive)+sum(controls%*%w_controls)
#So we can rewrite f(x)=2*(sum(positivecases%*%w_positive)- sum(everyone%*%w)= 2*sum(positivecases%*%w_positive)- n*mean(w)
#The mean is obviously mean(w) (2mean(x)-n). The variance is harder: If we take X=(x1,x2,...,x_n) the draws, the variance is Var(w%*%X) =w %*% cov(X) %*% w, because w is a contant.
#Now we need cov(X), cov(x_i,x_j)= E(x_i*x_j)-E(x_i)E(x_j) = P(x_i*x_j =1) - p^2 = P(x_i=1)*P(x_j =1/x_i=1) -p^2 = p* ((p*N-1)/(N-1))- p^2= (N*p^2 -p- Np^2+ p^2)/(N-1)= -p(1-p)/(N-1)
#So cov(x) is a matrix with p(1-p) on the diagonal and -p(1-p)/(N-1) everywhere else.  w%*%cov(X) is a vector p(1-p) (w1-sum(w_other)/(N-1)), w2-sum(w_other)/(N-1)),   ..., wn-sum(w_other)/(N-1)) 
#=p(1-p)/(N-1) ((N-1)w1-sum(w_other), (N-1)w-sum(w_other) ,..., (N-1)wn-sum(w_other))= p(1-p)/(N-1) (N w1-sum(w), N w2-sum(w) ,..., N wn-sum(w))
#So w%*%cov(X)%*%w= p(1-p)/(N-1) (N sum(wi^2)-sum(w)^2)
cumsumw=cumsum(w)
cumsumwsq=cumsum(w^2);
m=cumsumw*(2*p -1)
s=2*sqrt(p*(1-p)/(length(pheno)-1))*sapply(1:length(w),function(n)sqrt( length(pheno)*cumsumwsq[n] - (cumsumw[n])^2 ))
#Can also be written   n (length(pheno)*cummeanwsq[n] - n (cummeanw[n])^2 ) =  n (length(pheno)*var(w)+ cummeanw[n])^2  *(length(pheno)-n) )
#So the only new term is the addition of  p*(1-p) n (length(pheno)*var(w)) to the sqrt
return(list(mean=m,sd=s))
}


generate_es_normalizer2=function(w,pheno,nperm=2000){ #estimate variance under the null by permutations
p=length(which(pheno==1))/length(pheno);
#x from hypergeometric. Corresponding es score is f(x)=mean(w)(sum(positivecases)-sum(controls))=mean(w)(2x-n). Mean=mean(w)(2*mean(hypergeom)-n). 
cummeanw=cumsum(w)/(1:length(w));
m=cummeanw*(2*(1:length(w))*p -(1:length(w)))
resp=mapply(function(k)cumsum(2*(sample(pheno,length(pheno))-0.5)*w),1:nperm)
s=rep(1,length(w));
s[1:(length(w)/2)]=apply(resp[1:(length(w)/2),],1,sd)
#s[1:(length(w)/2)]=sapply(1:(length(w)/2),function(i)sqrt(sum((rand_es(i)-m[i])^2)/nperm ))
m=apply(resp,1,mean)
return(list(mean=m,sd=s))
}

best_chisq=function(pheno,e){
tp=1;js=1;
for (j in 10:(length(pheno)-10)){
ca=sum(pheno[order(e)][1:j])
t=chisq.test(cbind(c(ca,length(which(pheno==1))-ca),c(j-ca,length(which(pheno==0))-j+ca)))$p.value
if (t<tp){tp=t; js=j;}
}
print(paste(js,tp))
}


print_details=function(x){
  print(paste("indmax:",x$esind[x$direction],"cases:",x$ncas,"odd ratio:",format( x$oddratio ,digits=3)))
}

print_summary=function(x){
  print(paste("Score:",format(x$esperm$real,digits=3),"direction:",x$direction,"r:",format(x$oddcas,digits=3),"pval:",esperm$pval))
}


manhattan_plot=function(coord,genenames,pvalues, manfilename,annotatePval=0.0001,mirror=FALSE){#coord is a table with [gene,chr,pos]
#pvalues can be one vector or a matrix with two columns
pvalues1=pvalues;pvalues2=NULL;if(length(dim(pvalues))){if (dim(pvalues)[2]==2){pvalues1=pvalues[,1];pvalues2=pvalues[,2];mirror=TRUE;}}  
res=data.frame(SNP=genenames,CHR=26,BP=sample(1:10000000,length(genenames)),P=pvalues1);
chr=sapply(as.character(coord[,2]),function(x)substring(strsplit(x,"_")[[1]][1],4));
chr[chr=="X"]=23;chr[chr=="Y"]=24;chr[chr=="M"]=25;chr[chr=="Un"]=26;
ind=match(genenames,coord[,1]);res[which(!is.na(ind)),2]=as.numeric(chr[ind[which(!is.na(ind))]]);
res[which(!is.na(ind)),3]=coord[ind[which(!is.na(ind))],3];
res2=NULL;if (length(pvalues2)){res2=res;res2[,4]=pvalues2;}
ylim=c(0,0.5+max(-log10(pvalues1)));#if (mirror)ylim=rev(ylim);
library(qqman)
pdf(manfilename,width=11,height=9);
par(mfrow = c(2,1))
manhattan(res, chr="CHR", bp="BP", snp="SNP", p="P" ,annotatePval=annotatePval,annotateTop=FALSE,
          suggestiveline=FALSE,genomewideline=-log10(0.05/length(genenames)),ylim=ylim)
if (length(res2))manhattan(res2, chr="CHR", bp="BP", snp="SNP", p="P" ,annotatePval=annotatePval,annotateTop=FALSE,
          suggestiveline=FALSE,genomewideline=-log10(0.05/length(genenames)),ylim=rev(ylim))
#,ylab=expression('Limma   -log'[10]*'(p)')
dev.off();  
}

expr_plot=function(expr,pheno,name){
  col=c("black","red")
  pdf(paste(name,"_dens"));
  plot(density(expr[which(pheno==0)]),col=col[1],main="Distribution of expression levels for cases and controls");lines(density(expr[which(pheno==1)]),col=col[2]);
  legend("topleft",c("controls","cases"),fill=col)
  dev.off();
}

expr_plot2=function(expr,pheno,gsname){#without writing to file
  col=c("black","red")
  plot(density(expr[which(pheno==0)]),xlab=NA,col=col[1],main=paste(gsname,"expression levels"));lines(density(expr[which(pheno==1)]),col=col[2]);
  legend("topleft",c("controls","cases"),fill=col)
}

expr_plot_gg=function(expr,pheno,name,res_ref=NULL){
  library(ggplot2);
  col=c("black","red")
  ph=factor(pheno, levels=c(1,0),labels=c("Cases","Controls"))
  gg=ggplot(data.frame(Pheno=ph,Expression=expr), aes(x=Expression,fill=Pheno,linetype=Pheno))+geom_density(alpha=0.3)+ 
          ggtitle(name)+theme(axis.text=element_text(size=14), axis.title=element_text(size=14),plot.title = element_text(size=30,hjust = 0.5),legend.text=element_text(size=16),legend.title=element_text(size=18,hjust = 0.5)) ;
  if (length(res_ref)){
    ordexpr=sort(expr);
    exprval=ordexpr[res_ref$es_indmax]; ratio=res_ref$es_indmax/(8*length(pheno));#ratio is just to position text
    if (res_ref$es_sign==2){exprval=ordexpr[length(expr)+1-res_ref$es_indmax];ratio=1-ratio}
    gg=gg+geom_vline(xintercept=exprval,  color="black", size=1)#Add line for indmax 
    gg=gg+geom_text(x=(res_ref$es_sign-1.5)*Inf,y=Inf,label=paste(" ",res_ref$ncases,"cases /",res_ref$es_indmax," \n  O.R = ", format(res_ref$or,digits=3),"  "),angle=0,vjust=10,hjust=res_ref$es_sign-1)#Add text for ncase/indmax and oddratio  #ordexpr[round(ratio*length(expr))]
  }  
  print(gg);#removed dev.off() recently
}

expr_plot_v=function(expr,pheno,name){
  library(ggplot2);
  col=c("black","red")
  pdf(paste(name,"_violin.pdf"));
  print(ggplot(data.frame(ex=expr,phen=as.factor(pheno)), aes(x=phen, y=ex)) +geom_violin() );
  dev.off();
}

which_aberrant_old=function(finalexpri,res_esai){
  if (res_esai$direction==1){
  ipoint=order(finalexpri)[res_esai$esind[res_esai$direction]];
  return(which(finalexpri<= finalexpri[ipoint]));
  }else{
    ipoint=order(finalexpri)[ res_esai$esind[res_esai$direction]];
    return(which(finalexpri>=finalexpri[ipoint]));
  }
}

aberrant_expression=function(finalexpri,di,esmi){#Expression level corresponding to max ES
  if (di==1)return(finalexpri[order(finalexpri)[esmi]]);
  return(finalexpri[order(finalexpri,decreasing = T)[esmi]])
}

which_aberrant=function(finalexprj,finalexpri,di,esmi){#Return the individuals with finalexprj in the aberrant interval (Learned on finalexpri)
  ipoint=aberrant_expression(finalexpri,di,esmi)
  if (di==1)return(which(finalexprj<= ipoint))#direct sense
  if (di==2)return(which(finalexprj>= ipoint));#antisense
  return(NULL)
}

plot_2g=function(finalexpr,res_esa, i,j,name,genenames){
  pti=which_aberrant(finalexpr[i,],res_esa[[i]]);ptj=which_aberrant(finalexpr[j,],res_esa[[j]]);
  ptinter=intersect(pti,ptj);
  print(pti); print(ptj);print(ptinter);
  pdf(paste(name,"2G_",genenames[i],"_",genenames[j],".pdf",sep=""));
  plot(finalexpr[i,],finalexpr[j,],xlab=genenames[i],ylab=genenames[j]);
  title("Expression correlation")
  points(finalexpr[i,pti],finalexpr[j,pti],col="red");
  points(finalexpr[i,ptj],finalexpr[j,ptj],col="blue");
  if (length(ptinter))points(finalexpr[i,ptinter],finalexpr[j,ptinter],col="purple");
  dev.off();
}

minrange=function(x,mn=0,mx=8)c(min(x[1],mn),max(x[2],mx))

printvenn=function(x,y)print(paste(length(x),length(y),length(intersect(x,y))))

get_cor_quartiles=function(x,quart=c(-1,-0.5,-0.3,-0.1,0,0.1,0.3,0.5,1)){
  perc=rep(0,length(quart)-1);
  for (i in 1: (length(quart)-1)){perc[i]=sum(x<= quart[i+1] & x>=quart[i]);}
  return(perc/length(x))
}

mindaniel=function(pheno,expr,rep=100000){
  ord=order(expr,decreasing=FALSE);
  mindaniel0=function(pheno,o) (-expr[o[which(pheno[o]==1)[1]]]+ expr[o[which(pheno[o]==0)[1]]]);
  mindaniel_2d=function(pheno) max(c(mindaniel0(pheno,ord),-mindaniel0(pheno,rev(ord))));
  minperm=permute(rep,pheno,mindaniel_2d)
  return(list(pval=minperm$pval))
}

meandaniel=function(pheno,expr,rep=100000,t=0.1){
  ord=order(expr,decreasing=FALSE);
  bound=1:round(length(pheno)*t/2);
  revbound=((length(pheno)/2 )- round(length(pheno)*t/2)+1):(length(pheno)/2 );
  meandaniel0=function(pheno) t.test(expr[ord[which(pheno[ord]==1)[bound]]], expr[ord[which(pheno[ord]==0)[bound]]])$p.val;
  meandaniel1=function(pheno) t.test(expr[ord[which(pheno[ord]==1)[revbound]]], expr[ord[which(pheno[ord]==0)[revbound]]])$p.val;
  meandaniel_2d=function(pheno) -min(c(meandaniel0(pheno),meandaniel1(pheno)));
  meanperm=permute(rep,pheno,meandaniel_2d)
  return(list(pval=meanperm$pval))
}

testhypergeom=function(pheno,expr,rep=100000){
  ord=order(expr,decreasing=FALSE);
  n=length(which(pheno>0.5));
  testhypergeom0=function(pheno,o)sapply(1:n, function(x)dhyper(sum(pheno[o[1:x]]),n,length(which(pheno<=0.5)),x))
  testhypergeom_2d=function(pheno)-min(c(testhypergeom0(pheno,ord),testhypergeom0(pheno,rev(ord))));#- because permute look the upper end
  h_ind=c(which.min(testhypergeom0(pheno,ord)),which.min(testhypergeom0(pheno,rev(ord))));
  hyperperm=permute(rep,pheno,testhypergeom_2d)
  return(list( hypermin=-hyperperm$real,pval=hyperperm$pval, h_ind=h_ind))
}

plot_es_all_old=function(finalpheno,finalexpr,indg,info=NULL){#indg=indcrbn[1]
  #oldres=res_esa[indg];
  res= aziz.test(finalpheno,finalexpr[indg,],w0,rep=0,doall = T,trimmed=1)#recalculate with trimmed close to 1
    ggplot(data.frame(res$escurve))+
     geom_step(aes(x=seq_along(res$escurve), y=res$escurve)) +
     xlab("Individuals ordered by decreasing gene expression level") +
     ylab("Normalized enrichment score") + theme_bw()
}

plot_es_all=function(finalpheno,finalexpr,indg,info=NULL){#indg=indcrbn[1]; Run as plot_es_all(finalpheno,finalexpr,indg,show_res[indg,])
  res= aziz.test(finalpheno,finalexpr[indg,],w0,rep=0,doall = T,trimmed=1);
  gsea.es.profile=res$escurve;enrichment.score.range=range(c(gsea.es.profile,0));gsea.hit.indices=1+(1:length(gsea.es.profile));gsea.gene.set="CRBN";
  gsea.enrichment.score=info$es;gsea.rnk=list(metric=sort(finalexpr[indg,],decreasing=TRUE));metric.range=range(gsea.rnk$metric)
  def.par <- par(no.readonly = TRUE)
  dev.new(width = 3, height = 3)
  gsea.layout <- layout(matrix(c(1, 2, 3)), heights = c(2, 0.3, 0.15))
  layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)+1),#last was length(gsea.rnk$metric)
       c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Normalized Enrichment Score",
       ylim = enrichment.score.range,cex.lab=2.5,cex.axis=2,
       main = list(gsea.gene.set, font = 1, cex = 2.5),
       panel.first = {
         abline(h = seq(round(enrichment.score.range[1], digits = 1),
                        enrichment.score.range[2], 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  abline(v = info$es_indmax+1, lwd = 0.75,col="red")
  plot.coordinates <- par("usr")
  if(length(info)) {
  #  text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
  #       paste("FDR:", info$fdr_esa,"\nOur p-value:",info$pval2,  
  #             "\nLimma p-value:",info$limma, "\nWilcoxon p-value:",info$wilc), adj = c(0, 0))
  #} else {
  text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
  paste("FDR:", format(info$fdr_esa,digits=2) ,"\nNominal p-value:",format(info$pval2,digits=2),  
               "\nLimma:",format(info$limma,digits=2), "\nWilcoxon:",format(info$wilc,digits=2)), adj = c(1, 1),cex=1.8)  #}
  }
  text(length(gsea.rnk$metric) * 0.065, 0.5,"Interval \nof \nInterest",cex=1.8,col="red")
  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk$metric)))#Was
  abline(v = which(c(0,gsea.es.profile)<c(gsea.es.profile,0)), lwd = 0.75)
  
  par(mar = c(0, 5, 0, 2))
  rank.colors <- gsea.rnk$metric - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(gsea.rnk$metric)))
  box()
  #text(length(gsea.rnk$metric) / 2, 0.7,labels =gsea.gene.set)
  text(length(gsea.rnk$metric) * 0.01, 0.6,cex=2, "Overexpression", adj = c(0, NA))
  text(length(gsea.rnk$metric) * 0.99, 0.6,cex=2, "Underexpression", adj = c(1, NA))
  
}
