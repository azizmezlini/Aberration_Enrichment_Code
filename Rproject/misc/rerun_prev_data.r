if (length(commandArgs(TRUE))<4){print("Error: incorrect number of arguments");if(NA)print("Error");}
workdir=commandArgs(TRUE)[1];#"/data/mindinformatics/projects/testing_expression/"
reload=commandArgs(TRUE)[2];
k=as.numeric(commandArgs(TRUE)[3]);
k2=as.numeric(commandArgs(TRUE)[4]);


workdir="/data/mindinformatics/projects/testing_expression/"
reload="testexpr/prevresults/breast_mirna100.RData"
load(reload)
k=100;k2=0;

dir=workdir;
dir.create(paste(dir,"results/",name,"/",sep=""));
fname=paste("results/",name,"/",name,"_",k,"_",k2,sep="");
genes2=rownames(expr);patids=colnames(expr);
preprocess_param=list(minexpr=0,perc_minexpr=1,logtransformalpha=0);

source(paste0(workdir,"Rproject/main_disease_preprocess.r"))


#rerun previous data on Rstudio cluster; I used this Jan 2020 
if (FALSE){

dir=""
source(paste(dir,"Rproject/functions.r",sep=""))
fname="prevresults/breast_mirna100"
rg=1:(length(genenames))
nrep=10000000
#Fixing gene names for meth data that dosnt have them (sch and sch_blood)
temp=read.table("prevresults/meth_RA100.txt",header=TRUE);matchedgenes=temp[,2];names(matchedgenes)=temp[,1];names=matchedgenes[genenames]

source(paste(dir,"Rproject/main_disease_analyse.r",sep=""))
}

