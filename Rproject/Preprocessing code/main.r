library(GEOquery)
dir="/home/aziz/Desktop/aziz/testing_expression/"
dir="/hpf/largeprojects/agoldenb/aziz/testing_expression/"
source(paste(dir,"functions.r",sep=""))
#options('download.file.method.GEOquery'='auto')
genemap=read.table("gene_map.txt",header=TRUE)
gse=getGEO(filename=paste(dir,"expression-data/Lupus-GSE11907-GPL96_series_matrix.txt",sep=""))#Asthma-GSE69683_series_matrix.txt Flu-GSE52428_series_matrix.txt MultipleSclerosis-GSE41850_series_matrix.txt Islets-T2D-GSE81608_series_matrix.txt Alzheimer-GSE63063-GPL10558_series_matrix.txt.gz Alzheimer-GSE63063-GPL6947_series_matrix.txt.gz Anxiety-GSE61672_series_matrix.txt Kawasaki-GSE63881_series_matrix.txt Lupus-GSE11907-GPL96_series_matrix.txt Lupus-GSE11907-GPL97_series_matrix.txt RA-moderateToSevere-GSE74143_series_matrix.txt SLE-Tabalumab-GSE88887_series_matrix.txt T1D-Periphandwhole-GSE77350_series_matrix.txt Vaccines-GSE30101_series_matrix.txt
#GSE95674_series_matrix.txt.gz GSE95675_series_matrix.txt.gz
#GSE77164_series_matrix.txt.gz NEPALI CHILD dried blood
#GSE47862_series_matrix.txt.gz breast cancer heritable, peripheral blood
#GSE99039_series_matrix.txt.gz Parkinson disease   GSE95640_series_matrix.txt after-diet 
#breast_mirna_GSE73002_series_matrix.txt.gz  ovarian_mirna_GSE106817_series_matrix.txt.gz
#"GSE29676" proteomic alzheimer

expr=exprs(gse)
phenoData(gse)
names=featureNames(gse)
genes=fData(gse)[,13]
varLabels(gse)



#For data with problems
extractchr=function(x,str){if(grepl(str,x)){return(unlist(strsplit(unlist(strsplit(x,str))[2],"\t"))[1])}else return("-1")}
#extractchrold=function(x,str,b,e){if(grepl(str,x)){return(substring(unlist(strsplit(x,str))[2],b,e))}else return("-1")}

# Anxiety has issues, low sample size ~100/90 if we select nonmissing phenotype and gender.
#need more processing, gender is still missing for a third of the cases. ch1:age, chr1.1:bmi and gender, chr1.2:bmi and ethnicity, chr1.3 ethnicity,gad7,hyb batch;chr1.4 hyb_batch,case/control, gad7; chr1.5 case/control,hyb_batch,rin; chr1.6 hyb_batch,rin; chr1.7 rin
metapheno=paste(gse$characteristics_ch1,gse$characteristics_ch1.1,gse$characteristics_ch1.2,gse$characteristics_ch1.3,gse$characteristics_ch1.4,gse$characteristics_ch1.5,gse$characteristics_ch1.6,gse$characteristics_ch1.7,sep="\t")
data=data.frame(sapply(metapheno,extractchr,str="case\\/control: "), sapply(metapheno,extractchr,str="gad7 score: "), sapply(metapheno,extractchr,str="rin: "))
confounders=data.frame(gender=sapply(metapheno,extractchr,str="Sex: "), ethnicity=sapply(metapheno,extractchr,str="ethnicity: "), age=sapply(metapheno,extractchr,str="age: "),batch= sapply(metapheno,extractchr,str="hybridization batch: "));
table(data[which(as.character(data[,1])=="control"),2])#verify the missing are not at random and gad7 indicates the phenotype (cases are 5>= , contr are <=1, the missing are mostly 2,3,4)
cases=which(as.character(data[,1])=="case");contr=which(as.character(data[,1])=="control");
finalpheno=c(rep(1,length(cases)), rep(0,length(contr)));finalexpr=expr[,c(cases,contr)];

#MultipleSclerosis has issues, need more processing (features mixed together) Control 125  multiple sclerosis 690, 3 visits, gender and treatment are confounders (unless we only take baseline-untreated)
metapheno=paste(gse$characteristics_ch1,gse$characteristics_ch1.1,gse$characteristics_ch1.2,gse$characteristics_ch1.3,gse$characteristics_ch1.4,gse$characteristics_ch1.5,gse$characteristics_ch1.6,sep="\t")
data=data.frame(sapply(metapheno,extractchr,str="disease: "), sapply(metapheno,extractchr,str="visit: "))
confounders=data.frame(sapply(metapheno,extractchr,str="gender: "), sapply(metapheno,extractchr,str="treatment: "), sapply(metapheno,extractchr,str="tissue: ") );
#tissue and treatment have to be fused
confounders2=confounders[,2];ind=which(confounders[,2]=="whole blood");confounders2[ind]= confounders[ind,3];#treatment
cases=which(as.character(data[,1])=="multiple sclerosis" & data[,2]=="Baseline" & confounders2=="untreated");contr=which(as.character(data[,1])=="Control" & data[,2]=="Baseline" & confounders2=="untreated");
finalpheno=c(rep(1,length(cases)), rep(0,length(contr)));finalexpr=expr[,c(cases,contr)];


#SLE has issues, need more processing (features mixed together)(3 visits ch1.1: baseline, week 16,52  ; pat_id is in ch1). 
pheno=gse$characteristics_ch1.2 #Normal 60 SLE 1760+1326, confounders: batch, treatment, race,sex, age,region, sledai,adna,antidsdna (at baseline),c3base, c3 at baseline,c4base, c4 at baseline
#patient_id,time,,batch+treatment, batch+race, race+sex,sex+age, sledai_at_baseline,adna_at_baseline, antidsdna_at_baseline (iu)
metapheno=paste(gse$characteristics_ch1,gse$characteristics_ch1.1,gse$characteristics_ch1.2,gse$characteristics_ch1.3,gse$characteristics_ch1.4,gse$characteristics_ch1.5,gse$characteristics_ch1.6,gse$characteristics_ch1.7,gse$characteristics_ch1.8,gse$characteristics_ch1.9,gse$characteristics_ch1.10,gse$characteristics_ch1.11,sep="\t")
pheno=sapply(metapheno,extractchr,str="group: ")# Normal 60 SLE 1760+1326
time=sapply(metapheno,extractchr,str="time: ")
confoundersraw=data.frame(gender= sapply(metapheno,extractchr,str="Sex: "), batch = sapply(metapheno,extractchr,str="batch: "), race=sapply(metapheno,extractchr,str="race: "),treatment=sapply(metapheno,extractchr,str="treatment (q2w=once in 2 weeks; q4w=once in 4 weeks):"),age = sapply(metapheno,extractchr,str="age_at_baseline: "), region=sapply(metapheno,extractchr,str="region: "))
confounders=data.frame(gender=sapply(confoundersraw$gender,function(x){if (x=="-1")return(NA); return(as.numeric(x=="Female"))}), batch =as.numeric(sapply(as.character(confoundersraw$batch),function(x)(substr(x,nchar(x),nchar(x))=="L"))))
indcase=which(pheno=="IPD");indctrl=which(pheno=="CONTROL");
genes=fData(gse)[,11]
name="parkinson2";


#Vaccines has issues: 39 individuals , many time points, 10 different vaccines, 693 total, confounders: blood source (finger or vein), cell population, gender ch1.4, race: ch1.5, age, ethnicity,sample set


#Asthma   ch1 is phenotype ch1.1 and ch1.2 are gender and ethnicity. Values= cohort: Healthy, non-smoking 87 cohort: Moderate asthma, non-smoking  77 cohort: Severe asthma, non-smoking   246    cohort: Severe asthma, smoking 88    overall 87 controls, 246+77
pheno=gse$characteristics_ch1

#Flu (16 time points -40 individuals)
pheno=gse$characteristics_ch1.2 #clinical phenotype: Asymptomatic  clinical phenotype: Symptomatic. Confounders: Virus ch1.1 challenge virus: H1N1 challenge virus: H3N2, ch1.5 time point 

#Islets-T2D : donor_id (_ch1.1) indicates 12 controls and 6 patient ids ?
pheno=gse$characteristics_ch1.2#condition: non-diabetic  651        condition: T2D 941  , confounders: age ch1.3, ethnicity ch1.4, gender ch1.5, cell subtype ch1.6 (alpha,beta,delta,PP)

#Alzheimer1
pheno=gse$characteristics_ch1 #status: AD 139   status: CTL 134 (control) status: MCI 109 , confounders: age ch1.2, ethnicity ch1.1, gender ch1.3, extra:  included in gwas ch1.4

#Alzheimer2
pheno=gse$characteristics_ch1 # status: AD 145 status: CTL 104 (control) status: MCI  80, confounders: age ch1.2, ethnicity ch1.1,gender ch1.3, extra:  included in gwas ch1.4

#Kawasaki
pheno=gse$characteristics_ch1.2#  phenotype: Mild 146 phenotype: Moderate 9 phenotype: OFI 164  phenotype: Severe 22 , doi ch1.3#Confounders, phase: Acute 171 phase: Convalescent 170 ch1.1,Aneurysm,Dilated,Normal ch1.4 , ivg resistant/responsive ch1.5, age ch1.6, gender ch1.7

#Lupus1   #ch1.4: SLE 120 JIA 46 total 304, Confounders: age ch1, gender  ch1.1, (ethnicity ch1.2, race ch1.3 often not reported)
#Lupus2    #ch1.4: SLE 54 JIA 46 total 242, Confounders: age ch1, gender  ch1.1, (ethnicity ch1.2, race ch1.3 often not reported)
#PROBLEM: VERY FEW PROBES IN COMMON BETWEEN COHORTS (~200). BY GENE NAME 4000 matches 
gse=getGEO(filename="/home/aziz/Desktop/aziz/testing_expression/expression-data/Lupus-GSE11907-GPL96_series_matrix.txt")
gse2=getGEO(filename="/home/aziz/Desktop/aziz/testing_expression/expression-data/Lupus-GSE11907-GPL97_series_matrix.txt")
pheno=c(as.character(gse$characteristics_ch1.4),as.character(gse2$characteristics_ch1.4))
x=as.character(gse2$characteristics_ch1[1])
getage=function(x){spl=unlist(strsplit(x," ")); if (spl[3]%in%c("mths,","mth,"))return(as.numeric(spl[2])/12); if (spl[3]%in%c("yrs,","yr,"))return(as.numeric(spl[2])); return(NA)}
getgender=function(x){ if (x=="Gender: F")return(1);if (x=="Gender: M")return(0);return(NA);}
confounders=data.frame(age=sapply(c(as.character(gse$characteristics_ch1),as.character(gse2$characteristics_ch1)),getage), gender=sapply(c(as.character(gse$characteristics_ch1.1),as.character(gse2$characteristics_ch1.1)),getgender), cohort=c(rep(0,length(gse$characteristics_ch1.4)),rep(1,length(gse2$characteristics_ch1.4))))
indcase=which(pheno=="Illness: SLE");indctrl=which(pheno=="Illness: Healthy" | pheno=="Illness: JIA Systemic Onset" | pheno=="Illness: Melanoma" | pheno=="Illness: Melanoma M1a" | pheno=="Illness: Melanoma M1b" | pheno=="Illness: Melanoma M1c" | pheno=="Illness: Type 1 Diabetes" | pheno=="Illness: Hepatitis C" | pheno=="Illness: Hepatitis B"  |pheno=="Illness: Urinary Tract Infection");

#RA
pheno=gse$characteristics_ch1.2 # Negative 29 Positive 348, Confounders: age ch1.1, gender  ch1

#T1D
pheno=gse$characteristics_ch1.2 #condition: Basal 202  condition: CD4+ 102  condition: CD8+ 84 condition: PMA 6h 202, or hla_dr_type ch1.4, or type1_diabetes: 1 (107) type1_diabetes: 2(483) ch1.5
#Confounders:cell-type ch1.1, gender ch1.3, batchid ch1.7, extra id ch1.6

#Parkinson GSE99039
metapheno=paste(gse$characteristics_ch1,gse$characteristics_ch1.1,gse$characteristics_ch1.2,gse$characteristics_ch1.3,gse$characteristics_ch1.4,gse$characteristics_ch1.5,gse$characteristics_ch1.6,gse$characteristics_ch1.7,gse$characteristics_ch1.8,sep="\t")
pheno=sapply(metapheno,extractchr,str="disease\\ label: ")# IPD 205 CONTROL 233 GPD 41   GENETIC_UNAFFECTED 22
confoundersraw=data.frame(gender= sapply(metapheno,extractchr,str="Sex: "), batch = sapply(metapheno,extractchr,str="batch: "))
confounders=data.frame(gender=sapply(confoundersraw$gender,function(x){if (x=="-1")return(NA); return(as.numeric(x=="Female"))}), batch =as.numeric(sapply(as.character(confoundersraw$batch),function(x)(substr(x,nchar(x),nchar(x))=="L"))))
indcase=which(pheno=="IPD");indctrl=which(pheno=="CONTROL");
genes=fData(gse)[,11]
name="IPD";
#Other Confounders "age at 1st symptoms: " "age at exam: " mostly missing


#Child soldier 
pheno=gse$characteristics_ch1 #132 cases 122  controls , or resilience ch1.7 continuous
#confounders: gender gse$characteristics_ch1.1 ; age gse$characteristics_ch1.2 ; ethnic minority gse$characteristics_ch1.3 ; caste gse$characteristics_ch1.4  ; education gse$characteristics_ch1.5; ptsd: ch1.6; ch1.8 to ch1.15 blood counts

#heritable breast cancer  GSE47862
pheno=gse$characteristics_ch1.2 #developed_breast_cancer: yes 158 no 163
confounders=data.frame(cohort=as.numeric(gse$characteristics_ch1.3=="cohort: Ontario"))
indcase=which(pheno=="developed_breast_cancer: yes");indctrl=which(pheno=="developed_breast_cancer: no")
matchg=match(as.character(fData(gse)[,2]),as.character(genemap[,2]))
genes=as.character(genemap[matchg,1]);genes[which(is.na(matchg))]=paste("X",as.character(fData(gse)[which(is.na(matchg)),2]),sep="")
name="breast_expr";
#other confounders breast_cancer_family_history ch1.1,  carries_brca_mutation  ch1.4
#optional pheno=gse$characteristics_ch1.1; indcase=which(pheno=="breast_cancer_family_history: yes");indctrl=which(pheno=="breast_cancer_family_history: no")
#the people without family history did not have brca mutation , 93 carry mutation (43 cancer 50 no), 226 had family history (106 cancer 120 no), 95 have neither (52 cancer, 43 no)

#expr metastasis breast GSE48091. WARNING: matching controls with gse$"setnr:ch1" variable . Do I use it? setnr is a risk set (case cotrol pair or more)
pheno=gse$characteristics_ch1.2 #case-control status: 0  340   case-control status: 1  166
confounders=data.frame(batch=gse$characteristics_ch1); #maybe training/validation status if batches  gse$"setnr:ch1"
indcase=which(pheno=="case-control status: 1");indctrl=which(pheno=="case-control status: 0")
genes=as.character(fData(gse)[,3]);
name="metastatic_breast_expr";

#expr metastasis breast GSE102484   # not sure how to correct; smaller sample size; maybe use as validation
pheno=gse$characteristics_ch1.7#event_metastasis: 0  582 ; 1 101
#confounders age =gse$characteristics_ch1.3; 3 surgery: (BCT 328, MRM 451), 4 Stage:, 5 t stage:, 6 n stage:



#ovarian-mirna 
pheno=gse$characteristics_ch1;
confounders=NULL;# age missing for the majority of the cohort
indcase=which(pheno!="tumor stage (figo): NA");indctrl=which(pheno=="tumor stage (figo): NA")
genes=fData(gse)[,3]
name="ovarian_mirna";

#breast-mirna 
pheno=gse$characteristics_ch1;
confounders=NULL
indcase=which(pheno=="diagnosis: breast cancer");indctrl=which(pheno=="diagnosis: non-cancer")
genes=fData(gse)[,3]
name="breast_mirna";

#proteomic GSE29676; 2011; weird because the description says 50 AD cases not 350; might be measuring all differentials between all conditions
pheno=gse$characteristics_ch1.4 #Alzheimer's Disease  350; Breast Cancer 30;Older Control 120;Parkinson 29; Younger Control 80
getage=function(x){spl=unlist(strsplit(x," ")); return(as.numeric(spl[2])); }
getgene=function(x){spl=unlist(strsplit(x,"\\(|\\)")); return(spl[2]); return(NA)}
getgender=function(x){ if (x=="gender: Female")return(1);if (x=="gender: Male")return(0);return(NA);}
confounders=data.frame(age=sapply(as.character(gse$characteristics_ch1),getage),gender=sapply(as.character(gse$characteristics_ch1.1),getgender))
indcase=which(pheno=="disease state: Alzheimer's Disease");indctrl=which(pheno!="disease state: Alzheimer's Disease")
genes=sapply(as.character(fData(gse)[,4]),getgene);
name="Alz_proteomic";

#meth ID GSE89353 whole blood
gse=getGEO(filename="data/largedata/GSE89353_series_matrix.txt.gz")
pheno=gse$characteristics_ch1.2
getgender=function(x){ if (x=="gender: Female")return(1);if (x=="gender: Male")return(0);return(NA);}
confounders=data.frame(gender=sapply(as.character(gse$characteristics_ch1.1),getgender))
indctrl=which(sapply(as.character(pheno),function(x)(grepl("*ther$",x)| grepl("*sibling$",x))));indcase=(1:length(pheno))[-indctrl]
genes=as.character(fData(gse)[,22]);
name="meth_ID";
expr=exprs(gse)
save.image(paste(dir,name,".RData",sep=""))

#methylation controls GSE36064 GSE40279 GSE42861 GSE53045 


#meth schizo GSE74193 dorsolateral prefrontal cortex; no data or variable names
pheno=gse$characteristics_ch1.5;
getage=function(x){spl=unlist(strsplit(x,": ")); return(as.numeric(spl[2])); }
getgender=function(x){spl=unlist(strsplit(x,":")); return(spl[2]);}
genderpred=sapply(as.character(gse$characteristics_ch1.11),getgender);
drop=sapply(as.character(gse$characteristics_ch1.2),function(x){if(x=="dropsample (whether to remove the sample for failing quality control): FALSE")return(0);return(NA)});
bestqc=sapply(as.character(gse$characteristics_ch1.12),function(x){if(x=="bestqc (best sample to use when more than 1 array were run on the same subject/brnum): TRUE")return(0);return(NA)});
confounders=data.frame(gender=sapply(as.character(gse$characteristics_ch1.3),getgender), race=gse$characteristics_ch1.4, batch=gse$characteristics_ch1.1, age=sapply(as.character(gse$characteristics_ch1.6),getage),comp1=sapply(as.character(gse$characteristics_ch1.13),getage),comp2=sapply(as.character(gse$characteristics_ch1.14),getage),comp3=sapply(as.character(gse$characteristics_ch1.15),getage),comp4=sapply(as.character(gse$characteristics_ch1.16),getage),comp5=sapply(as.character(gse$characteristics_ch1.17),getage));
confounders[which(is.na(drop) | is.na(bestqc) | genderpred!=confounders[,1] ),1]=NA;#To remove bad samples
indcase=which(pheno=="group: Schizo");indctrl=which(pheno=="group: Control");
name="meth_sch";
exprraw=read.csv(gzfile("/hpf/largeprojects/agoldenb/aziz/testing_expression/data/largedata/meth_sch_cortex_GSE74193_GEO_procData.csv.gz"))
genes=exprraw[,1];
expr=exprraw[,-which( (1:ncol(exprraw))%%2 ==1)];exprraw=NULL;
save.image(paste(dir,name,".RData",sep=""))

#meth schizo GSE80417 whole blood
pheno=gse$characteristics_ch1.2;
getage=function(x){spl=unlist(strsplit(x,": ")); return(as.numeric(spl[2])); }
confounders=data.frame(gender=gse$characteristics_ch1,age=sapply(as.character(gse$characteristics_ch1.1),getage)) 
indcase=which(pheno=="disease status: 1");indctrl=which(pheno=="disease status: 2");
name="meth_sch_blood";
exprraw=read.csv(gzfile("/hpf/largeprojects/agoldenb/aziz/testing_expression/data/largedata/meth_sch_blood_GSE80417_normalizedBetas.csv.gz"))
genes=exprraw[,1];
expr=exprraw[,-1];exprraw=NULL;
save.image(paste(dir,name,".RData",sep=""))

#meth RA  GSE42861
pheno=gse$"disease state:ch1";
confounders=data.frame(gender=gse$"gender:ch1", smoking=gse$"smoking status:ch1",age=as.numeric(gse$"age:ch1"))
indcase=which(pheno=="rheumatoid arthritis");indctrl=which(pheno=="Normal")
genes=as.character(fData(gse)[,1]);
name="meth_RA";
names=fData(gse)[,22]
save.image(paste(dir,name,".RData",sep=""))


#expr IBD GSE73094; WARNING: where is the tissue info? How to correct  table(gse$description)
pheno=gse$characteristics_ch1.1  #patient_id=gse$characteristics_ch1 # disease: CD 608   disease: UC 331; or 
gettissue=function(x){spl=unlist(strsplit(x,"_")); return(spl[length(spl)]);}
tissue=as.factor(sapply(as.character(gse$description),gettissue))
confounders=data.frame(batch=gse$characteristics_ch1.3, tissue=tissue);
noninflamed=which(gse$characteristics_ch1.2=="inflamed (0=no, 1=yes): 0");
uniqnoninflamed=noninflamed[which(!duplicated(gse$characteristics_ch1[noninflamed]))];
indcase=which(pheno=="disease: CD");indctrl=which(pheno=="disease: UC");
indcase=indcase[which(indcase %in% uniqnoninflamed)];indctrl=indctrl[which(indctrl %in% uniqnoninflamed)]
#pheno=gse$characteristics_ch1.2 # inflamed (0=no, 1=yes): 0  609  ; 1 374
#confounders=data.frame(batch=gse$characteristics_ch1.3,disease=gse$characteristics_ch1.1,tissue=tissue)
#indcase=which(pheno=="inflamed (0=no, 1=yes): 1");indctrl=which(pheno=="inflamed (0=no, 1=yes): 0")
genes=as.character(fData(gse)[,3]);
name="IBD_inflammation";

#Asthma exacerbation GSE19301

save.image(paste(dir,name,".RData",sep=""))
genes[which(genes=="")]=paste("X",as.character(fData(gse)[which(genes==""),1]),sep="")
source("main_disease.r")
