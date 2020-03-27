
dir="~/Downloads/testing_expression/";name="Alz_miRNA";
source(paste(dir,"Rproject/functions.r",sep=""))
#genemap=read.table("gene_map.txt",header=TRUE)
library(GEOquery)
#options('download.file.method.GEOquery'='auto')
gse=getGEO("GSE120584")$GSE120584_series_matrix.txt.gz
#phenoData(gse)
#varLabels(gse)
expr=exprs(gse)
genes=fData(gse)[,3];genes2=featureNames(gse);
patids=gse$title; names(patids)<- colnames(expr);
pheno=gse$characteristics_ch1  #patient_id=gse$title # diagnosis: AD 1021  diagnosis: DLB 169 diagnosis: NC 288; age,sex,apoe4
getage=function(x){spl=unlist(strsplit(x," ")); return(as.numeric(spl[2])); }
confounders=data.frame(gender=gse$characteristics_ch1.2,age=sapply(as.character(gse$characteristics_ch1.1),getage), apoe=gse$characteristics_ch1.3);
indcase=which(pheno=="diagnosis: AD");indctrl=which(pheno=="diagnosis: NC");
save.image(paste(dir,"data/GSE120584_loaded.RData",sep=""),version=2)
