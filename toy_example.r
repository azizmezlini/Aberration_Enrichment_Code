
#Load the code : change to the right directory path
source("aziz_test/Rproject/functions.r")

#Create case control/labels: 100 cases followed by 100 controls
y=c(rep(1,100),rep(0,100));

#create variable to be tested
x=runif(200);

#Using our test
res=aziz.test(y,x,rep=10000)#10000 permutations

#Print pvalue
print(res$pval)