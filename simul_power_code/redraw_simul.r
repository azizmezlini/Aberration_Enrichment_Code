
cBlind8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")#black,orange,sky blue, bluish green,yellow, blue,vermillion,reddish purple
chosencol=c("white","red","blue","orange","cyan","pink","darkgreen","purple","brown","black")#none, our test, our+T-test, our+ Levene, our+t-test+Levene,our+t-test+wilcoxon(not levene), t-test only, t-test+wilcoxon,all
chosencol=c("white",cBlind8[c(7,6,3)],"darkgray",cBlind8[c(8,5,4)],"brown","black")

m_main=list();length(m_main)<- 2;
m_main[[1]]=rbind(c(9,9,9,9,9,9,9,9,9),#hard
             c(4,9,9,9,9,9,9,9,9),
             c(1,4,9,9,9,9,9,9,9),
             c(1,3,4,4,4,9,9,9,9),
             c(0,1,1,3,4,4,4,9,9),
             c(0,0,1,1,1,3,3,4,4),
             c(0,0,0,0,0,1,1,1,3))
m_main[[2]]=rbind(c(0,3,4,4,9,9,9,9,9),#soft
             c(0,1,1,4,4,4,9,9,9),
             c(0,0,0,1,1,3,4,4,9),
             c(0,0,0,0,0,1,1,4,4),
             c(0,0,0,0,0,0,0,1,3),
             c(0,0,0,0,0,0,0,0,1),
             c(0,0,0,0,0,0,0,0,0))
for (i in 1:length(m_main)){rownames(m_main[[i]])<- rev(c(0.02,0.035,0.05,0.075,0.1,0.15,0.2)); colnames(m_main[[i]])<- c(100,200,300,400,500,600,800,1200,2000);}


m_broad=rbind(c(5,5,5,5,5,5,5,5),
              c(0,2,5,5,5,5,5,9),
              c(0,0,1,5,9,9,9,9),
              c(0,0,0,0,0,2,9,9),
              c(0,0,0,0,0,0,1,4))
rownames(m_broad)<- rev(c(0.2,0.3,0.5,0.75,1)); colnames(m_broad)<- c(20,30,50,75,100,150,200,300);


m_robust=list();length(m_robust)<- 4;
m_robust[[1]]=rbind(c(0,4,4,9,9,9,9,9,9),
                c(0,1,4,4,4,4,9,9,9),
                c(0,0,1,3,3,4,4,4,9),
                c(0,0,0,1,1,1,3,4,4),
                c(0,0,0,0,0,1,1,3,4),
                c(0,0,0,0,0,0,0,1,3),
                c(0,0,0,0,0,0,0,0,1))
m_robust[[2]]=rbind(c(0,3,4,4,9,9,9,9,9),
                c(0,1,1,4,4,4,9,9,9),
                c(0,0,0,1,1,3,4,4,9),
                c(0,0,0,0,0,1,1,4,4),
                c(0,0,0,0,0,0,0,1,3),
                c(0,0,0,0,0,0,0,0,1),
                c(0,0,0,0,0,0,0,0,0))
m_robust[[3]]=rbind(c(0,0,1,2,4,9,9,9,9),
                c(0,0,0,1,1,1,4,9,9),
                c(0,0,0,0,0,0,1,1,4),
                c(0,0,0,0,0,0,0,1,1),
                c(0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0))
m_robust[[4]]=rbind(c(0,0,0,0,1,1,5,9,9),
                c(0,0,0,0,0,0,1,2,9),
                c(0,0,0,0,0,0,0,0,1),
                c(0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0),
                c(0,0,0,0,0,0,0,0,0))
for (i in 1:length(m_robust)){rownames(m_robust[[i]])<- rev(c(0.02,0.035,0.05,0.075,0.1,0.15,0.2)); colnames(m_robust[[i]])<- c(100,200,300,400,500,600,800,1200,2000);}


m_low_d=list();length(m_low_d)<- 4;
m_low_d[[1]]=rbind(c(2,5,5,5,9,9,9,9),
                c(0,0,0,1,2,5,9,9),
                c(0,0,0,0,0,2,5,9),
                c(0,0,0,0,0,0,0,5),
                c(0,0,0,0,0,0,0,0))
m_low_d[[2]]=rbind(c(0,0,5,5,5,5,5,9),
               c(0,0,0,0,0,6,5,5),
               c(0,0,0,0,0,0,0,5),
               c(0,0,0,0,0,0,0,0),
               c(0,0,0,0,0,0,0,0))
m_low_d[[3]]=rbind(c(0,0,0,7,7,5,5,5),
               c(0,0,0,0,0,0,0,5),
               c(0,0,0,0,0,0,0,0),
               c(0,0,0,0,0,0,0,0),
               c(0,0,0,0,0,0,0,0))
m_low_d[[4]]=rbind(c(0,0,0,0,0,0,7,5),
               c(0,0,0,0,0,0,0,0),
               c(0,0,0,0,0,0,0,0),
               c(0,0,0,0,0,0,0,0),
               c(0,0,0,0,0,0,0,0))
for (i in 1:length(m_low_d)){rownames(m_low_d[[i]])<- rev(c(0.02,0.05,0.075,0.1,0.2)); colnames(m_low_d[[i]])<- c(200,300,500,750,1000,1500,2500,5000);}

m_precise=rbind(c(0,1,1,1,3,4,4,4),#hard
                c(0,0,0,1,1,1,1,3),
                c(0,0,0,0,0,0,1,1),
                c(0,0,0,0,0,0,0,0))
rownames(m_precise)<- rev(c(0.005,0.01,0.02,0.03)); colnames(m_precise)<- c(800,1200,1800,2700,4000,6000,9000,14000);


noMargins <- function(..., topkey = FALSE, rightkey = FALSE) {
  nmlist <- list(
    layout.heights = list(
      top.padding = 0,
      main.key.padding = 0,
      key.axis.padding = 0,
      axis.xlab.padding = 0,
      xlab.key.padding = 0,
      key.sub.padding = 0,
      bottom.padding = 0.5
    ),
    layout.widths = list(
      left.padding = 0,
      key.ylab.padding = 0,
      ylab.axis.padding = 0,
      axis.key.padding = 0,
      right.padding = 0.6
    ),
    axis.components = list(
      top = list(pad1 = ifelse(topkey, 2, 1), pad2 = ifelse(topkey, 2, 0)), # padding above top axis
      right = list(pad1 = 0, pad2 = 0)
    )
  )
  c(nmlist, ...)
}

resdirname="Downloads/main_d3-2/redrawn/"
rev_t=function(m)t(m[rev(1:nrow(m)),])
pdf(paste(resdirname,"broad_",3,".pdf",sep=""))
print(levelplot(rev_t(m_broad),xlab="n",ylab="r",col.regions=chosencol, at=c(seq(0,10,0.999),10.01),colorkey=F,interpolate=T))
dev.off()

pdf(paste(resdirname,"precise_",3,".pdf",sep=""))
print(levelplot(rev_t(m_precise),xlab="n",ylab="r",col.regions=chosencol, at=c(seq(0,10,0.999),10.01),colorkey=F,interpolate=T))
dev.off()
#layout(matrix(c(1,2), 1, 2, byrow = TRUE),widths=0.5,heights=c(0.5,0.5))
#par(mfrow=c(2,1)) ; par(mar = c(0, 5, 2, 2))


m=m_main;pdf(paste(resdirname,"main_",3,".pdf",sep=""))
plots=list();length(plots)<- length(m);for (i in 1:length(m)){
  plots[[i]]=levelplot(rev_t(m[[i]]),xlab="n",ylab="r",col.regions=chosencol, at=c(seq(0,10,0.999),10.01),colorkey=F,interpolate=T,par.settings=noMargins())
}
grid.arrange(plots[[1]],plots[[2]], ncol=1)
dev.off()

m=m_robust;pdf(paste(resdirname,"robust_",3,".pdf",sep=""))
ds=c(4,3,2,1.5)
plots=list();length(plots)<- length(m);for (i in 1:length(m)){
  plots[[i]]=levelplot(rev_t(m[[i]]),xlab="n",ylab="r",col.regions=chosencol, at=c(seq(0,10,0.999),10.01),colorkey=F,interpolate=T,par.settings=noMargins(),main=paste("d=",ds[i]))
}
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]], ncol=2)
dev.off()

m=m_low_d;pdf(paste(resdirname,"low_d_",3,".pdf",sep=""));
ds=c(1,0.7,0.5,0.3)
plots=list();length(plots)<- length(m);for (i in 1:length(m)){
  plots[[i]]=levelplot(rev_t(m[[i]]),xlab="n",ylab="r",font.lab=2,cex.lab=2,col.regions=chosencol, at=c(seq(0,10,0.999),10.01),colorkey=F,interpolate=T,par.settings=noMargins(),main=paste("d=",ds[i]))
}
grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]], ncol=2)
dev.off()

