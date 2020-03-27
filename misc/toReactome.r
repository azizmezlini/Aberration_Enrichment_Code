#Code selecting associated genes and copying to clipboard for easy copy/paste to reactome webservice
x=read.table("Desktop/testing_expression/prevres/Alzheimer_ctrl_30_0.txt",header=TRUE)
ind=which(x$fdr_esa<0.1 & x$r<0.3);#ind=which(x$approx<0.05/nrow(x) );
#ind=which(x$fdr_w<0.1 );#ind=which(x$wilc<0.05/nrow(x) );
sort(table(x$genenames[ind]),decreasing = TRUE)[1:10]
names=unique(unlist(sapply(as.character(x$genenames[ind]),function(y)strsplit(y,';')[[1]])))
library(clipr)
write_clip(names)
