if (length(commandArgs(TRUE))<4){print("Error: incorrect number of arguments");if(NA)print("Error");}
workdir=commandArgs(TRUE)[1];
name=commandArgs(TRUE)[2];#name="Alz_miRNA"; | name="Alz_miRNA_Rosmap";
k=as.numeric(commandArgs(TRUE)[3]);
k2=as.numeric(commandArgs(TRUE)[4]);
if (length(commandArgs(TRUE))>4){revperm=as.numeric(commandArgs(TRUE)[5]);}else{revperm=0;}#Controls whether we reverse or randomize cases and controls. 1 means reverse, >1 permute


if (name=="Alz_miRNA")load(paste0(workdir,"data/GSE120584_loaded.RData"))
  
preprocess_param=list(minexpr=-500,perc_minexpr=1,logtransformalpha=0);

if (revperm==1){t=indcase; indcase=indctrl; indctrl=t;name=paste(name,"rev",sep="");}#Switching cases and controls
dir=workdir;
dir.create(paste(dir,"results/",name,"/",sep=""));
fname=paste("results/",name,"/",name,"_",k,"_",k2,sep="")

source(paste0(workdir,"Rproject/main_disease_preprocess.r"))



