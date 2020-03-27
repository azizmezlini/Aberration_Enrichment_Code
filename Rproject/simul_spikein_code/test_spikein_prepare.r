
#combine Alzheimer data:
library(GEOquery)
library(peer)
#options('download.file.method.GEOquery'='auto')
gse=getGEO(filename="/home/aziz/Desktop/aziz/testing_expression/expression-data/Alzheimer-GSE63063-GPL10558_series_matrix.txt.gz")
gse2=getGEO(filename="/home/aziz/Desktop/aziz/testing_expression/expression-data/Alzheimer-GSE63063-GPL6947_series_matrix.txt.gz")

interprobes=intersect(fData(gse)[,1],fData(gse2)[,1]);
match1=match(interprobes,fData(gse)[,1]);match2=match(interprobes,fData(gse2)[,1])
genes=fData(gse)[match1,13];genenames=as.character(genes);
expr=cbind(exprs(gse)[match1,],exprs(gse2)[match2,])
varLabels(gse)
pheno=c(as.character(gse$characteristics_ch1),as.character(gse2$characteristics_ch1))
confounders=data.frame(cohort=c(rep(1,ncol(exprs(gse))),rep(2,ncol(exprs(gse2)) )),
                       ethnicity=c(gse$characteristics_ch1.1,gse2$characteristics_ch1.1),age=c(gse$characteristics_ch1.2,gse2$characteristics_ch1.2),gender=c(gse$characteristics_ch1.3,gse2$characteristics_ch1.3))
phenoData(gse)
names=featureNames(gse)

indcase=which(pheno=="status: AD" | pheno=="status: MCI");indctrl=which(pheno=="status: CTL");pheno0=rep(-1,length(pheno));pheno0[indcase]=1;pheno0[indctrl]=0;
#modcombat=model.matrix(~1,data=pheno0)
#combat_data=ComBat(dat=expr,batch=confounders[,1]-1,par.prior=TRUE)
exprscaled=t(scale(t(expr)));
k=30
model = PEER()
PEER_setPhenoMean(model,t(as.matrix(expr)))
PEER_setNk(model,k)
PEER_setCovariates(model, as.matrix(confounders))
PEER_update(model)

factors = PEER_getX(model)
weights = PEER_getW(model)
residuals = PEER_getResiduals(model)
#peerres=peer_analysis(expr,k);

selected=which(pheno0==0 )
ctrlexpr=t(residuals)[,selected];rownames(ctrlexpr)<- fData(gse)[match1,1]
write.table(ctrlexpr,"healthy_expr.txt")
