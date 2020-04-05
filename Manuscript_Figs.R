#Setting input and output folders:
basefolder <- "/Users/leekennedy-shaffer/Documents/Harvard/Research/Adjusted ICCs/Code/Lin2018/Manuscript Figures"
figfolder <- basefolder

#Read in functions that implement variance estimation and helpful packages/functions:
source(paste0(basefolder,"/VarCompFns.R"))
library(gee)
expit <- function(x) exp(x)/(1+exp(x))
library(LaplacesDemon)

#Colors for Plots
Res.cols <- c(1,"#1f78b4","#a6cee3")
SS.cols <- c(1,"#33a02c","#b2df8a","#756bb1","#bcbddc")

## Reading in the data:
Lin  <- read.csv(file=paste0(basefolder,"/Lin2018Protozoa.csv"))

## Keeping only those with observed results:
Lin2 <- Lin[Lin$delta_prot=="Observed",]

#Creating cluster-level electricity access variable:
Lin2$ClustElec <- sapply(Lin2$clusterid, FUN=function(x) mean(Lin2$elec[Lin2$clusterid==x]=="Has electricity"))
Lin2$ElecBin <- ifelse(Lin2$ClustElec>.5,1,0)

#Summaries of cluster sizes:
NumObs <- length(Lin2$clusterid)
Clustids <- unique(Lin2$clusterid)
NumClusts <- length(Clustids)
ClustSizes <- sapply(X=Clustids, FUN=function(x) length(Lin2$personid[Lin2$clusterid==x]))
muAdj <- mean(ClustSizes) - 1 ## One less than mean, used for Generalized Poisson Approximation
var <- (sd(ClustSizes))^2
cv <- sqrt(var)/(muAdj+1)
nvec <- 1:50 ## Potential Cluster Sizes used for Generalized Poisson Approximation

#Generalized Poisson Approximation of Cluster Size Distribution:
hist(ClustSizes, breaks=seq(from=.5,to=20.5,by=.5))
points(x=seq(.75,19.75,by=1), y=719*dgpois(seq(0,19,by=1), lambda=muAdj, omega=1-sqrt(muAdj/var)))

#Cluster sizes by electricity access:
ClustSizes.E1 <- sapply(X=unique(Lin2$clusterid[Lin2$ElecBin==1]), FUN=function(x) length(Lin2$personid[Lin2$clusterid==x]))
ClustSizes.E0 <- sapply(X=unique(Lin2$clusterid[Lin2$ElecBin==0]), FUN=function(x) length(Lin2$personid[Lin2$clusterid==x]))

#Data sets for the control and treatment conditions used in our analysis:
CtrlOnly <- Lin2[Lin2$tr=="Control",]
Trt <- Lin2[Lin2$tr %in% c("WSH","Nutrition + WSH"),]
K.total <- length(unique(c(CtrlOnly$clusterid,Trt$clusterid)))

#Fitting unadjusted GEE models within each stratum:
Trt.gee <- gee(posgi~1, id=clusterid, data=Trt, family=binomial, corstr="exchangeable")
CtrlOnly.gee <- gee(posgi~1, id=clusterid, data=CtrlOnly, family=binomial, corstr="exchangeable")
Trt.E1.gee <- gee(posgi~1, id=clusterid, data=Trt, subset=ElecBin==1, family=binomial, corstr="exchangeable")
Trt.E0.gee <- gee(posgi~1, id=clusterid, data=Trt, subset=ElecBin==0, family=binomial, corstr="exchangeable")
Ctrl.E1.gee <- gee(posgi~1, id=clusterid, data=CtrlOnly, subset=ElecBin==1, family=binomial, corstr="exchangeable")
Ctrl.E0.gee <- gee(posgi~1, id=clusterid, data=CtrlOnly, subset=ElecBin==0, family=binomial, corstr="exchangeable")

# ICCs for each stratum:
icc.ctrlonly <- CtrlOnly.gee$working.correlation[1,2]
icc.trt <- Trt.gee$working.correlation[1,2]
icc.trt.E1 <- Trt.E1.gee$working.correlation[1,2]
icc.trt.E0 <- Trt.E0.gee$working.correlation[1,2]
icc.ctrl.E1 <- Ctrl.E1.gee$working.correlation[1,2]
icc.ctrl.E0 <- Ctrl.E0.gee$working.correlation[1,2]
icc.Eratio <- icc.trt.E1/icc.trt.E0

# Prevalence rates for each stratum:
pi.ctrlonly <- unname(expit(CtrlOnly.gee$coefficients))
pi.trt <- unname(expit(Trt.gee$coefficients))
pi.trt.E1 <- unname(expit(Trt.E1.gee$coefficients))
pi.trt.E0 <- unname(expit(Trt.E0.gee$coefficients))
pi.ctrl.E1 <- unname(expit(Ctrl.E1.gee$coefficients))
pi.ctrl.E0 <- unname(expit(Ctrl.E0.gee$coefficients))
pi.Eratio <- pi.trt.E1/pi.trt.E0

#Number of clusters by stratum:
length(unique(CtrlOnly$clusterid))
length(unique(Trt$clusterid))
length(unique(CtrlOnly$clusterid[CtrlOnly$ElecBin==0]))
length(unique(Trt$clusterid[Trt$ElecBin==0]))
length(unique(CtrlOnly$clusterid[CtrlOnly$ElecBin==1]))
length(unique(Trt$clusterid[Trt$ElecBin==1]))

#Length of vectors for plots (higher number means more resolution on plots)
len <- 199

###### One Binary Cluster-Level Covariate ######

### Varying Rho_1 ###
df.VaryRho1 <- data.frame(rho.ratio=seq(from=log(1/4),to=log(4),length.out=len))
df.VaryRho1$rho1s <- exp(df.VaryRho1$rho.ratio)*icc.ctrlonly
ResMat <- sapply(X=df.VaryRho1$rho1s,
                 FUN=function(x) Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                      rho0=icc.ctrlonly,rho1=x,
                                      pi0=pi.ctrlonly,pi1=pi.trt,
                                      K.tot=K.total))
df.VaryRho1$corvar <- unlist(ResMat["Cor",])
df.VaryRho1$indvar <- unlist(ResMat["Ind",])
df.VaryRho1$exchvar <- unlist(ResMat["Exch",])
ResMat.eq <- sapply(X=ResMat["rho.star",],
                    FUN=function(x) Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                         rho0=x,rho1=x,
                                         pi0=pi.ctrlonly,pi1=pi.trt,
                                         K.tot=K.total))
df.VaryRho1$rho.star <- unlist(ResMat["rho.star",])
df.VaryRho1$corvar.eq <- unlist(ResMat.eq["Cor",])
df.VaryRho1$indvar.eq <- unlist(ResMat.eq["Ind",])
df.VaryRho1$exchvar.eq <- unlist(ResMat.eq["Exch",])
df.VaryRho1$exchvar.eq.rho0 <- rep(Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                        rho0=icc.ctrlonly,rho1=icc.ctrlonly,
                                        pi0=pi.ctrlonly,pi1=pi.trt,
                                        K.tot=K.total)$Exch,
                                   len)
ResMat.eq.rho1 <- sapply(X=df.VaryRho1$rho1s,
                         FUN=function(x) Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                              rho0=x,rho1=x,
                                              pi0=pi.ctrlonly,pi1=pi.trt,
                                              K.tot=K.total))
df.VaryRho1$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])


### Varying Pi_1 ###
df.VaryPi1 <- data.frame(pi.ratio=seq(from=log(1/2),to=log(2),length.out=len))
df.VaryPi1$pi1s <- exp(df.VaryPi1$pi.ratio)*pi.ctrlonly
ResMat <- sapply(X=df.VaryPi1$pi1s,
                 FUN=function(x) Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                      rho0=icc.ctrlonly,rho1=icc.trt,
                                      pi0=pi.ctrlonly,pi1=x,
                                      K.tot=K.total))
df.VaryPi1$corvar <- unlist(ResMat["Cor",])
df.VaryPi1$indvar <- unlist(ResMat["Ind",])
df.VaryPi1$exchvar <- unlist(ResMat["Exch",])
ResMat.eq <- sapply(X=df.VaryPi1$pi1s,
                    FUN=function(x) Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                         rho0=unlist(ResMat["rho.star",1]),rho1=unlist(ResMat["rho.star",1]),
                                         pi0=pi.ctrlonly,pi1=x,
                                         K.tot=K.total))
df.VaryPi1$rho.star <- unlist(ResMat["rho.star",])
df.VaryPi1$corvar.eq <- unlist(ResMat.eq["Cor",])
df.VaryPi1$indvar.eq <- unlist(ResMat.eq["Ind",])
df.VaryPi1$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq.rho0 <- sapply(X=df.VaryPi1$pi1s,
                         FUN=function(x) Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                              rho0=icc.ctrlonly,rho1=icc.ctrlonly,
                                              pi0=pi.ctrlonly,pi1=x,
                                              K.tot=K.total))
df.VaryPi1$exchvar.eq.rho0 <- unlist(ResMat.eq.rho0["Exch",])
ResMat.eq.rho1 <- sapply(X=df.VaryPi1$pi1s,
                         FUN=function(x) Vars(n0vec=ClustSizes,n1vec=ClustSizes,
                                              rho0=icc.trt,rho1=icc.trt,
                                              pi0=pi.ctrlonly,pi1=x,
                                              K.tot=K.total))
df.VaryPi1$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])

Plot.Res4(resdf1=df.VaryRho1, resdf2=df.VaryPi1,
          xlab1=bquote(rho["1"]/rho["0"]), xlab2=bquote(pi["1"]/pi["0"]),
          xvar1="rho.ratio", xvar2="pi.ratio",
          xtitle1=bquote("Ratio of"~rho["1"]~"to"~rho["0"]), xtitle2=bquote("Ratio of"~pi["1"]~"to"~pi["0"]),
          fileout="RhoPiRatios",
          truth1=log(icc.trt/icc.ctrlonly), truth2=log(pi.trt/pi.ctrlonly),
          ylimit1=c(0.075,0.125), ylimit2=c(90,101),
          colors=Res.cols,
          logX1=TRUE, logX2=TRUE)


### Varying Cluster Size CVs ###
df.VaryCVs <- data.frame(CVs=seq(from=.1,to=.7,length.out=len))
df.VaryCVs$Vars <- (df.VaryCVs$CVs^2)*(muAdj+1)^2
nwts <- sapply(X=df.VaryCVs$Vars,
               FUN=function(x) sapply(X=nvec,
                                      FUN=function(y) dgpois(x=y-1,lambda=muAdj, omega=1-sqrt(muAdj/x))))
nwts <- ifelse(is.nan(nwts),0,nwts)
nwts <- prop.table(nwts,2)
df.VaryCVs$truemeans <- apply(nwts, 2, FUN=function(x) sum(x*nvec))
df.VaryCVs$truevars <- apply(nwts, 2, FUN=function(x) sum(x*(nvec-sum(x*nvec))^2))
df.VaryCVs$trueCVs <- round(sqrt(df.VaryCVs$truevars)/df.VaryCVs$truemeans,digits=3)
ResMat <- apply(X=nwts, MARGIN=2,
                FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                            rho0=icc.ctrlonly,rho1=icc.trt,
                                            pi0=pi.ctrlonly,pi1=pi.trt,
                                            n0wts=x,n1wts=x,
                                            K.tot=K.total)))
df.VaryCVs$corvar <- unlist(ResMat["Cor",])
df.VaryCVs$indvar <- unlist(ResMat["Ind",])
df.VaryCVs$exchvar <- unlist(ResMat["Exch",])
ResMat.eq <- apply(X=nwts, MARGIN=2,
                   FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                               rho0=unlist(ResMat["rho.star",1]),rho1=unlist(ResMat["rho.star",1]),
                                               pi0=pi.ctrlonly,pi1=pi.trt,
                                               n0wts=x,n1wts=x,
                                               K.tot=K.total)))
df.VaryCVs$rho.star <- unlist(ResMat["rho.star",])
df.VaryCVs$corvar.eq <- unlist(ResMat.eq["Cor",])
df.VaryCVs$indvar.eq <- unlist(ResMat.eq["Ind",])
df.VaryCVs$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq.rho0 <- apply(X=nwts, MARGIN=2,
                        FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                                    rho0=icc.ctrlonly,rho1=icc.ctrlonly,
                                                    pi0=pi.ctrlonly,pi1=pi.trt,
                                                    n0wts=x,n1wts=x,
                                                    K.tot=K.total)))
df.VaryCVs$exchvar.eq.rho0 <- unlist(ResMat.eq.rho0["Exch",])
ResMat.eq.rho1 <- apply(X=nwts, MARGIN=2,
                        FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                                    rho0=icc.trt,rho1=icc.trt,
                                                    pi0=pi.ctrlonly,pi1=pi.trt,
                                                    n0wts=x,n1wts=x,
                                                    K.tot=K.total)))
df.VaryCVs$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])


### Varying Cluster Size Means ###
df.VaryMeans <- data.frame(Means=seq(from=5,to=20,length.out=len))
df.VaryMeans$Vars <- (df.VaryMeans$Means*cv)^2
nwts <- apply(X=as.matrix(df.VaryMeans), MARGIN=1,
              FUN=function(x) sapply(X=nvec,
                                     FUN=function(y) dgpois(x=y-1,lambda=x["Means"]-1,
                                                            omega=1-sqrt((x["Means"]-1)/x["Vars"]))))
nwts <- ifelse(is.nan(nwts),0,nwts)
nwts <- prop.table(nwts,2)
df.VaryMeans$truemeans <- apply(nwts, 2, FUN=function(x) sum(x*nvec))
df.VaryMeans$truevars <- apply(nwts, 2, FUN=function(x) sum(x*(nvec-sum(x*nvec))^2))
df.VaryMeans$trueCVs <- round(sqrt(df.VaryMeans$truevars)/df.VaryMeans$truemeans,digits=3)
df.VaryMeans$truemeans <- round(df.VaryMeans$truemeans, digits=3)
ResMat <- apply(X=nwts, MARGIN=2,
                FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                            rho0=icc.ctrlonly,rho1=icc.trt,
                                            pi0=pi.ctrlonly,pi1=pi.trt,
                                            n0wts=x,n1wts=x,
                                            K.tot=K.total)))
df.VaryMeans$corvar <- unlist(ResMat["Cor",])
df.VaryMeans$indvar <- unlist(ResMat["Ind",])
df.VaryMeans$exchvar <- unlist(ResMat["Exch",])
ResMat.eq <- apply(X=nwts, MARGIN=2,
                   FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                               rho0=unlist(ResMat["rho.star",1]),rho1=unlist(ResMat["rho.star",1]),
                                               pi0=pi.ctrlonly,pi1=pi.trt,
                                               n0wts=x,n1wts=x,
                                               K.tot=K.total)))
df.VaryMeans$rho.star <- unlist(ResMat["rho.star",])
df.VaryMeans$corvar.eq <- unlist(ResMat.eq["Cor",])
df.VaryMeans$indvar.eq <- unlist(ResMat.eq["Ind",])
df.VaryMeans$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq.rho0 <- apply(X=nwts, MARGIN=2,
                        FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                                    rho0=icc.ctrlonly,rho1=icc.ctrlonly,
                                                    pi0=pi.ctrlonly,pi1=pi.trt,
                                                    n0wts=x,n1wts=x,
                                                    K.tot=K.total)))
df.VaryMeans$exchvar.eq.rho0 <- unlist(ResMat.eq.rho0["Exch",])
ResMat.eq.rho1 <- apply(X=nwts, MARGIN=2,
                        FUN=function(x) unlist(Vars(n0vec=nvec,n1vec=nvec,
                                                    rho0=icc.trt,rho1=icc.trt,
                                                    pi0=pi.ctrlonly,pi1=pi.trt,
                                                    n0wts=x,n1wts=x,
                                                    K.tot=K.total)))
df.VaryMeans$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])


Plot.Res4(resdf1=df.VaryMeans, resdf2=df.VaryCVs,
          xlab1="Mean Cluster Size", xlab2="CV of Cluster Size Distribution",
          xvar1="Means", xvar2="CVs",
          xtitle1="Mean Cluster Size", xtitle2="CV of Cluster Size Distribution",
          fileout="SizeDistn",
          truth1=mean(ClustSizes), truth2=sqrt(var)/(muAdj+1),
          ylimit1=c(0.075,0.125), ylimit2=c(90,101),
          colors=Res.cols)


## Sample Size Plotting:
setEPS()
postscript(file=paste0(figfolder,"/SS_Single4.eps"),
           width=12, height=12, paper="special")
par(mfrow=c(2,2))
Plot.SS.Full(resdf=df.VaryRho1,bquote(rho["1"]/rho["0"]), xvar="rho.ratio",
             title=bquote("a."~"Ratio of "~rho["1"]~"to"~rho["0"]), truth=log(icc.trt/icc.ctrlonly),
             ylimit=c(90,110), colors=SS.cols, logX=TRUE)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                        expression("Common Exch."~rho[0]),
                        expression("Independence"),
                        expression("Common Exch."~rho[1]),
                        expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df.VaryPi1,bquote(pi["1"]/pi["0"]), xvar="pi.ratio",
             title=bquote("b."~"Ratio of "~pi["1"]~"to"~pi["0"]), truth=log(pi.trt/pi.ctrlonly),
             ylimit=c(90,110), colors=SS.cols, logX=TRUE)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                        expression("Common Exch."~rho[0]),
                        expression("Independence"),
                        expression("Common Exch."~rho[1]),
                        expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df.VaryMeans,xlab=bquote("Mean"~"Cluster Size"), xvar="Means",
             title=bquote("c."~"Mean Cluster Size"), truth=muAdj+1,
             ylimit=c(90,110), colors=SS.cols)
legend(x="top",legend=c(expression("Correctly Specified"),
                        expression("Common Exch."~rho[0]),
                        expression("Independence"),
                        expression("Common Exch."~rho[1]),
                        expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df.VaryCVs,xlab=bquote("CV"~"of Cluster Size"), xvar="CVs",
             title=bquote("d."~"CV of Cluster Size"), truth=sqrt(var)/(muAdj+1),
             ylimit=c(90,110), colors=SS.cols)
legend(x="top",legend=c(expression("Correctly Specified"),
                        expression("Common Exch."~rho[0]),
                        expression("Independence"),
                        expression("Common Exch."~rho[1]),
                        expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
dev.off()




###### Two Binary Cluster-Level Covariates ######
### Varying Rho Ratio ###
df2.RhoRatio <- data.frame(rho.ratio=seq(from=log(1/4),to=log(4),length.out=len))
df2.RhoRatio$rho00s <- icc.ctrl.E0
df2.RhoRatio$rho10s <- icc.trt.E0
df2.RhoRatio$rho01s <- df2.RhoRatio$rho00s*exp(df2.RhoRatio$rho.ratio)
df2.RhoRatio$rho11s <- df2.RhoRatio$rho10s*exp(df2.RhoRatio$rho.ratio)
ResMat <- apply(X=df2.RhoRatio, MARGIN=1,
                 FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n01vec=ClustSizes, n10vec=ClustSizes, n11vec=ClustSizes,
                                       pi00=pi.ctrl.E0, pi10=pi.trt.E1, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                       rho00=x["rho00s"], rho10=x["rho10s"], rho01=x["rho01s"], rho11=x["rho11s"],
                                       n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                       n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                       K.tot=K.total)))
df2.RhoRatio$corvar <- unlist(ResMat["Cor",])
df2.RhoRatio$indvar <- unlist(ResMat["Ind",])
df2.RhoRatio$exchvar <- unlist(ResMat["Exch",])
df2.RhoRatio$rho.star <- unlist(ResMat["rho.star",])
ResMat.eq <- sapply(X=df2.RhoRatio$rho.star,
                    FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                 pi00=pi.ctrl.E0, pi10=pi.trt.E1, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                rho00=x, rho10=x, rho01=x, rho11=x,
                                                n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                K.tot=K.total)))
df2.RhoRatio$corvar.eq <- unlist(ResMat.eq["Cor",])
df2.RhoRatio$indvar.eq <- unlist(ResMat.eq["Ind",])
df2.RhoRatio$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq.rho0 <- apply(X=df2.RhoRatio, MARGIN=1,
                           FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                        pi00=pi.ctrl.E0, pi10=pi.trt.E1, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                        rho00=0.4*x["rho00s"]+0.6*x["rho01s"], rho10=0.4*x["rho00s"]+0.6*x["rho01s"], 
                                                        rho01=0.4*x["rho00s"]+0.6*x["rho01s"], rho11=0.4*x["rho00s"]+0.6*x["rho01s"],
                                                        n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                        n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                        K.tot=K.total)))
df2.RhoRatio$exchvar.eq.rho0 <- unlist(ResMat.eq.rho0["Exch",])
ResMat.eq.rho1 <- apply(X=df2.RhoRatio, MARGIN=1,
                           FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                        pi00=pi.ctrl.E0, pi10=pi.trt.E1, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                        rho00=0.4*x["rho10s"]+0.6*x["rho11s"], rho10=0.4*x["rho10s"]+0.6*x["rho11s"], 
                                                        rho01=0.4*x["rho10s"]+0.6*x["rho11s"], rho11=0.4*x["rho10s"]+0.6*x["rho11s"],
                                                        n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                        n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                        K.tot=K.total)))
df2.RhoRatio$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])


### Varying Pi Ratio ###
df2.pi.Eratio <- data.frame(pi.ratio=seq(from=log(1/2),to=log(2),length.out=len))
df2.pi.Eratio$pi00s <- pi.ctrl.E0
df2.pi.Eratio$pi10s <- pi.trt.E0
df2.pi.Eratio$pi01s <- df2.pi.Eratio$pi00s*exp(df2.pi.Eratio$pi.ratio)
df2.pi.Eratio$pi11s <- df2.pi.Eratio$pi10s*exp(df2.pi.Eratio$pi.ratio)
ResMat <- apply(X=df2.pi.Eratio, MARGIN=1,
                FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                             pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                             rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                             n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                             n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                             K.tot=K.total)))
df2.pi.Eratio$corvar <- unlist(ResMat["Cor",])
df2.pi.Eratio$indvar <- unlist(ResMat["Ind",])
df2.pi.Eratio$exchvar <- unlist(ResMat["Exch",])
df2.pi.Eratio$rho.star <- unlist(ResMat["rho.star",])
ResMat.eq <- apply(X=df2.pi.Eratio, MARGIN=1,
                    FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                 pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                 rho00=df2.pi.Eratio$rho.star[1], rho10=df2.pi.Eratio$rho.star[1], rho01=df2.pi.Eratio$rho.star[1], rho11=df2.pi.Eratio$rho.star[1],
                                                 n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                 n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                 K.tot=K.total)))
df2.pi.Eratio$rho.star <- rep(unname(unlist(ResMat["rho.star",1])),dim(df2.pi.Eratio)[1])
df2.pi.Eratio$corvar.eq <- unlist(ResMat.eq["Cor",])
df2.pi.Eratio$indvar.eq <- unlist(ResMat.eq["Ind",])
df2.pi.Eratio$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq.rho0 <- apply(X=df2.pi.Eratio, MARGIN=1,
                   FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                rho00=icc.ctrlonly, rho10=icc.ctrlonly, rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                                n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                K.tot=K.total)))
df2.pi.Eratio$exchvar.eq.rho0 <- unlist(ResMat.eq.rho0["Exch",])
ResMat.eq.rho1 <- apply(X=df2.pi.Eratio, MARGIN=1,
                          FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                       pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                       rho00=icc.trt, rho10=icc.trt, rho01=icc.trt, rho11=icc.trt,
                                                       n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                       n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                       K.tot=K.total)))
df2.pi.Eratio$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])

Plot.Res4(resdf1=df2.RhoRatio, resdf2=df2.pi.Eratio,
          xlab1=bquote(rho["01"]/rho["00"]==rho["11"]/rho["10"]), xlab2=bquote(pi["01"]/pi["00"]==pi["11"]/pi["10"]),
          xvar1="rho.ratio", xvar2="pi.ratio",
          xtitle1=bquote("Ratio of"~rho["01"]~"to"~rho["00"]~"and"~rho["11"]~"to"~rho["10"]), xtitle2=bquote("Ratio of"~pi["01"]~"to"~pi["00"]~"and"~pi["11"]~"to"~pi["10"]),
          fileout="2_RhoPiRatios",
          truth1=log(icc.Eratio), truth2=log(pi.Eratio),
          ylimit1=c(0.075,0.125), ylimit2=c(85,100),
          colors=Res.cols,
          logX1=TRUE, logX2=TRUE,
          truth1b=log(icc.ctrl.E1/icc.ctrl.E0), truth2b=log(pi.ctrl.E1/pi.ctrl.E0))


### Varying Cluster Size CV ###
df2.VaryCVs <- data.frame(CVs=seq(from=.1,to=.7,length.out=len))
df2.VaryCVs$Vars <- (df2.VaryCVs$CVs^2)*(muAdj+1)^2
nwts <- sapply(X=df2.VaryCVs$Vars,
               FUN=function(x) sapply(X=nvec,
                                      FUN=function(y) dgpois(x=y-1,lambda=muAdj, omega=1-sqrt(muAdj/x))))
nwts <- ifelse(is.nan(nwts),0,nwts)
nwts <- prop.table(nwts,2)
df2.VaryCVs$truemeans <- apply(nwts, 2, FUN=function(x) sum(x*nvec))
df2.VaryCVs$truevars <- apply(nwts, 2, FUN=function(x) sum(x*(nvec-sum(x*nvec))^2))
df2.VaryCVs$trueCVs <- round(sqrt(df2.VaryCVs$truevars)/df2.VaryCVs$truemeans,digits=3)
ResMat <- apply(X=nwts, MARGIN=2,
                FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                             pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                             rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                             n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                             K.tot=K.total)))
df2.VaryCVs$corvar <- unlist(ResMat["Cor",])
df2.VaryCVs$indvar <- unlist(ResMat["Ind",])
df2.VaryCVs$exchvar <- unlist(ResMat["Exch",])
df2.VaryCVs.rho.star <- unname(unlist(ResMat["rho.star",1]))
ResMat.eq <- apply(X=nwts, MARGIN=2,
                   FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                                pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                rho00=df2.VaryCVs.rho.star, rho10=df2.VaryCVs.rho.star, rho01=df2.VaryCVs.rho.star, rho11=df2.VaryCVs.rho.star,
                                                n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                                K.tot=K.total)))
df2.VaryCVs$rho.star <- rep(df2.VaryCVs.rho.star,dim(df2.VaryCVs)[1])
df2.VaryCVs$corvar.eq <- unlist(ResMat.eq["Cor",])
df2.VaryCVs$indvar.eq <- unlist(ResMat.eq["Ind",])
df2.VaryCVs$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq.rho0 <- apply(X=nwts, MARGIN=2,
                   FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                                pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                rho00=icc.ctrlonly, rho10=icc.ctrlonly, rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                                n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                                K.tot=K.total)))
df2.VaryCVs$exchvar.eq.rho0 <- unlist(ResMat.eq.rho0["Exch",])
ResMat.eq.rho1 <- apply(X=nwts, MARGIN=2,
                          FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                                       pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                       rho00=icc.trt, rho10=icc.trt, rho01=icc.trt, rho11=icc.trt,
                                                       n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                                       K.tot=K.total)))
df2.VaryCVs$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])

### Varying Cluster Size Mean ###
df2.VaryMeans <- data.frame(Means=seq(from=5,to=20,length.out=len))
df2.VaryMeans$Vars <- (df2.VaryMeans$Means*cv)^2
nwts <- apply(X=as.matrix(df2.VaryMeans), MARGIN=1,
              FUN=function(x) sapply(X=nvec,
                                     FUN=function(y) dgpois(x=y-1,lambda=x["Means"]-1,
                                                            omega=1-sqrt((x["Means"]-1)/x["Vars"]))))
nwts <- ifelse(is.nan(nwts),0,nwts)
nwts <- prop.table(nwts,2)
df2.VaryMeans$truemeans <- apply(nwts, 2, FUN=function(x) sum(x*nvec))
df2.VaryMeans$truevars <- apply(nwts, 2, FUN=function(x) sum(x*(nvec-sum(x*nvec))^2))
df2.VaryMeans$trueCVs <- round(sqrt(df2.VaryMeans$truevars)/df2.VaryMeans$truemeans,digits=3)
df2.VaryMeans$truemeans <- round(df2.VaryMeans$truemeans, digits=3)
ResMat <- apply(X=nwts, MARGIN=2,
                FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                             pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                             rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                             n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                             K.tot=K.total)))
df2.VaryMeans$corvar <- unlist(ResMat["Cor",])
df2.VaryMeans$indvar <- unlist(ResMat["Ind",])
df2.VaryMeans$exchvar <- unlist(ResMat["Exch",])
df2.VaryMeans.rho.star <- unname(unlist(ResMat["rho.star",1]))
ResMat.eq <- apply(X=nwts, MARGIN=2,
                   FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                                pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                rho00=df2.VaryMeans.rho.star, rho10=df2.VaryMeans.rho.star, rho01=df2.VaryMeans.rho.star, rho11=df2.VaryMeans.rho.star,
                                                n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                                K.tot=K.total)))
df2.VaryMeans$rho.star <- rep(df2.VaryMeans.rho.star,dim(df2.VaryMeans)[1])
df2.VaryMeans$corvar.eq <- unlist(ResMat.eq["Cor",])
df2.VaryMeans$indvar.eq <- unlist(ResMat.eq["Ind",])
df2.VaryMeans$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq.rho0 <- apply(X=nwts, MARGIN=2,
                   FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                                pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                rho00=icc.ctrlonly, rho10=icc.ctrlonly, rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                                n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                                K.tot=K.total)))
df2.VaryMeans$exchvar.eq.rho0 <- unlist(ResMat.eq.rho0["Exch",])
ResMat.eq.rho1 <- apply(X=nwts, MARGIN=2,
                          FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                                       pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                       rho00=icc.trt, rho10=icc.trt, rho01=icc.trt, rho11=icc.trt,
                                                       n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                                       K.tot=K.total)))
df2.VaryMeans$exchvar.eq.rho1 <- unlist(ResMat.eq.rho1["Exch",])

Plot.Res4(resdf1=df2.VaryMeans, resdf2=df2.VaryCVs,
          xlab1="Mean Cluster Size", xlab2="CV of Cluster Size Distribution",
          xvar1="truemeans", xvar2="trueCVs",
          xtitle1="Mean Cluster Size", xtitle2="CV of Cluster Size Distribution",
          fileout="2_SizeDistn",
          truth1=mean(ClustSizes), truth2=sqrt(var)/(muAdj+1),
          ylimit1=c(0.075,0.125), ylimit2=c(85,100),
          colors=Res.cols)


## Sample Size Plotting:
setEPS()
postscript(file=paste0(figfolder,"/SS_2Covar4.eps"),
           width=12, height=12, paper="special")
par(mfrow=c(2,2))
Plot.SS.Full(resdf=df2.RhoRatio,bquote(rho["01"]/rho["00"]==rho["11"]/rho["10"]), xvar="rho.ratio",
             title=bquote("a."~"Ratio of"~rho["01"]~"to"~rho["00"]~"and"~rho["11"]~"to"~rho["10"]),
             truth=log(icc.Eratio),
             ylimit=c(80,120), colors=SS.cols, logX=TRUE, truthb=log(icc.ctrl.E1/icc.ctrl.E0))
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.pi.Eratio,bquote(pi["01"]/pi["00"]==pi["11"]/pi["10"]), xvar="pi.ratio",
        title=bquote("b."~"Ratio of"~pi["01"]~"to"~pi["00"]~"and"~pi["11"]~"to"~pi["10"]),
        truth=log(pi.Eratio),
        ylimit=c(80,120), colors=SS.cols, logX=TRUE, truthb=log(pi.ctrl.E1/pi.ctrl.E0))
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.VaryMeans,xlab=bquote("Mean"~"Cluster Size"), xvar="truemeans",
        title=bquote("c."~"Mean Cluster Size"), truth=mean(ClustSizes),
        ylimit=c(80,120), colors=SS.cols)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.VaryCVs,xlab=bquote("CV"~"of Cluster Size Distribution"), xvar="trueCVs",
        title=bquote("d."~"CV of Cluster Size Distribution"), truth=sqrt(var)/(muAdj+1),
        ylimit=c(80,120), colors=SS.cols)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
dev.off()




###################

### Varying Ratio of ICCs due to Electrification Differentially 
###  Between Intervention (Arm 1) and Control (Arm 0) Arms
df2.RhoArm1 <- data.frame(rho.ratio=seq(from=log(1/4),to=log(4),length.out=len))
df2.RhoArm1$rho00s <- icc.ctrl.E0
df2.RhoArm1$rho10s <- icc.trt.E0
df2.RhoArm1$rho01s <- icc.ctrl.E1
df2.RhoArm1$rho11s <- df2.RhoArm1$rho10s*exp(df2.RhoArm1$rho.ratio)

df2.RhoArm0 <- data.frame(rho.ratio=seq(from=log(1/4),to=log(4),length.out=len))
df2.RhoArm0$rho00s <- icc.ctrl.E0
df2.RhoArm0$rho10s <- icc.trt.E0
df2.RhoArm0$rho01s <- df2.RhoArm0$rho00s*exp(df2.RhoArm0$rho.ratio)
df2.RhoArm0$rho11s <- icc.trt.E1

ResMat <- apply(X=df2.RhoArm1, MARGIN=1,
                FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                             pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                             rho00=x["rho00s"], rho10=x["rho10s"], rho01=x["rho01s"], rho11=x["rho11s"],
                                             n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                             n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                             K.tot=K.total)))
df2.RhoArm1$corvar <- unlist(ResMat["Cor",])
df2.RhoArm1$indvar <- unlist(ResMat["Ind",])
df2.RhoArm1$exchvar <- unlist(ResMat["Exch",])
df2.RhoArm1$rho.star <- unlist(ResMat["rho.star",])
ResMat0 <- apply(X=df2.RhoArm0, MARGIN=1,
                 FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                              pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                              rho00=x["rho00s"], rho10=x["rho10s"], rho01=x["rho01s"], rho11=x["rho11s"],
                                              n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                              n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                              K.tot=K.total)))
df2.RhoArm0$corvar <- unlist(ResMat0["Cor",])
df2.RhoArm0$indvar <- unlist(ResMat0["Ind",])
df2.RhoArm0$exchvar <- unlist(ResMat0["Exch",])
df2.RhoArm0$rho.star <- unlist(ResMat0["rho.star",])
ResMat.eq <- sapply(X=df2.RhoArm1$rho.star,
                    FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                 pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                 rho00=x, rho10=x, rho01=x, rho11=x,
                                                 n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                 n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                 K.tot=K.total)))
df2.RhoArm1$corvar.eq <- unlist(ResMat.eq["Cor",])
df2.RhoArm1$indvar.eq <- unlist(ResMat.eq["Ind",])
df2.RhoArm1$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq0 <- sapply(X=df2.RhoArm0$rho.star,
                     FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                  pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                  rho00=x, rho10=x, rho01=x, rho11=x,
                                                  n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                  n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                  K.tot=K.total)))
df2.RhoArm0$corvar.eq <- unlist(ResMat.eq0["Cor",])
df2.RhoArm0$indvar.eq <- unlist(ResMat.eq0["Ind",])
df2.RhoArm0$exchvar.eq <- unlist(ResMat.eq0["Exch",])
ResMat.eq1.rho0 <- apply(X=df2.RhoArm1, MARGIN=1,
                         FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                      pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                      rho00=.4*x["rho00s"]+.6*x["rho01s"], rho10=.4*x["rho00s"]+.6*x["rho01s"], 
                                                      rho01=.4*x["rho00s"]+.6*x["rho01s"], rho11=.4*x["rho00s"]+.6*x["rho01s"],
                                                      n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                      n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                      K.tot=K.total)))
df2.RhoArm1$exchvar.eq.rho0 <- unlist(ResMat.eq1.rho0["Exch",])
ResMat.eq1.rho1 <- apply(X=df2.RhoArm1, MARGIN=1,
                         FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                      pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                      rho00=.4*x["rho10s"]+.6*x["rho11s"], rho10=.4*x["rho10s"]+.6*x["rho11s"], 
                                                      rho01=.4*x["rho10s"]+.6*x["rho11s"], rho11=.4*x["rho10s"]+.6*x["rho11s"],
                                                      n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                      n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                      K.tot=K.total)))
df2.RhoArm1$exchvar.eq.rho1 <- unlist(ResMat.eq1.rho1["Exch",])
ResMat.eq0.rho0 <- apply(X=df2.RhoArm0, MARGIN=1,
                     FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                  pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                  rho00=.4*x["rho00s"]+.6*x["rho01s"], rho10=.4*x["rho00s"]+.6*x["rho01s"], 
                                                  rho01=.4*x["rho00s"]+.6*x["rho01s"], rho11=.4*x["rho00s"]+.6*x["rho01s"],
                                                  n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                  n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                  K.tot=K.total)))
df2.RhoArm0$exchvar.eq.rho0 <- unlist(ResMat.eq0.rho0["Exch",])
ResMat.eq0.rho1 <- apply(X=df2.RhoArm0, MARGIN=1,
                         FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                      pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                      rho00=.4*x["rho10s"]+.6*x["rho11s"], rho10=.4*x["rho10s"]+.6*x["rho11s"], 
                                                      rho01=.4*x["rho10s"]+.6*x["rho11s"], rho11=.4*x["rho10s"]+.6*x["rho11s"],
                                                      n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                      n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                      K.tot=K.total)))
df2.RhoArm0$exchvar.eq.rho1 <- unlist(ResMat.eq0.rho1["Exch",])

Plot.Res4(resdf1=df2.RhoArm0, resdf2=df2.RhoArm1,
          xlab1=bquote(rho["01"]/rho["00"]), xlab2=bquote(rho["11"]/rho["10"]),
          xvar1="rho.ratio", xvar2="rho.ratio",
          xtitle1=bquote("Ratio of"~rho["01"]~"to"~rho["00"]), xtitle2=bquote("Ratio of"~rho["11"]~"to"~rho["10"]),
          fileout="2_ArmsRhoRatio",
          truth1=log(icc.ctrl.E1/icc.ctrl.E0), truth2=log(icc.trt.E1/icc.trt.E0),
          ylimit1=c(0.075,0.125), ylimit2=c(85,100),
          colors=Res.cols, logX1=TRUE, logX2=TRUE)


### Varying Ratio of ICCs due to Electrification Differentially 
###  Between Intervention (Arm 1) and Control (Arm 0) Arms
df2.PiArm1 <- data.frame(pi.ratio=seq(from=log(1/2),to=log(2),length.out=len))
df2.PiArm1$pi00s <- pi.ctrl.E0
df2.PiArm1$pi10s <- pi.trt.E0
df2.PiArm1$pi01s <- pi.ctrl.E1
df2.PiArm1$pi11s <- df2.PiArm1$pi10s*exp(df2.PiArm1$pi.ratio)

df2.PiArm0 <- data.frame(pi.ratio=seq(from=log(1/2),to=log(2),length.out=len))
df2.PiArm0$pi00s <- pi.ctrl.E0
df2.PiArm0$pi10s <- pi.trt.E0
df2.PiArm0$pi01s <- df2.PiArm0$pi00s*exp(df2.PiArm0$pi.ratio)
df2.PiArm0$pi11s <- pi.trt.E1

ResMat <- apply(X=df2.PiArm1, MARGIN=1,
                FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                             pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                             rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                             n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                             n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                             K.tot=K.total)))
df2.PiArm1$corvar <- unlist(ResMat["Cor",])
df2.PiArm1$indvar <- unlist(ResMat["Ind",])
df2.PiArm1$exchvar <- unlist(ResMat["Exch",])
df2.PiArm1$rho.star <- unlist(ResMat["rho.star",])
ResMat0 <- apply(X=df2.PiArm0, MARGIN=1,
                 FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                              pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                              n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                              K.tot=K.total)))
df2.PiArm0$corvar <- unlist(ResMat0["Cor",])
df2.PiArm0$indvar <- unlist(ResMat0["Ind",])
df2.PiArm0$exchvar <- unlist(ResMat0["Exch",])
df2.PiArm0$rho.star <- unlist(ResMat0["rho.star",])
ResMat.eq <- apply(X=df2.PiArm1, MARGIN=1,
                    FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                 pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                 rho00=x["rho.star"], rho10=x["rho.star"], rho01=x["rho.star"], rho11=x["rho.star"],
                                                 n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                 n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                 K.tot=K.total)))
df2.PiArm1$corvar.eq <- unlist(ResMat.eq["Cor",])
df2.PiArm1$indvar.eq <- unlist(ResMat.eq["Ind",])
df2.PiArm1$exchvar.eq <- unlist(ResMat.eq["Exch",])
ResMat.eq0 <- apply(X=df2.PiArm0, MARGIN=1,
                     FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                  pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                  rho00=x["rho.star"], rho10=x["rho.star"], rho01=x["rho.star"], rho11=x["rho.star"],
                                                  n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                  n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                  K.tot=K.total)))
df2.PiArm0$corvar.eq <- unlist(ResMat.eq0["Cor",])
df2.PiArm0$indvar.eq <- unlist(ResMat.eq0["Ind",])
df2.PiArm0$exchvar.eq <- unlist(ResMat.eq0["Exch",])
ResMat.eq1.rho0 <- apply(X=df2.PiArm1, MARGIN=1,
                         FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                      pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                      rho00=icc.ctrlonly, rho10=icc.ctrlonly, rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                                      n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                      n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                      K.tot=K.total)))
df2.PiArm1$exchvar.eq.rho0 <- unlist(ResMat.eq1.rho0["Exch",])
ResMat.eq1.rho1 <- apply(X=df2.PiArm1, MARGIN=1,
                         FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                      pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                      rho00=icc.trt, rho10=icc.trt, rho01=icc.trt, rho11=icc.trt,
                                                      n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                      n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                      K.tot=K.total)))
df2.PiArm1$exchvar.eq.rho1 <- unlist(ResMat.eq1.rho1["Exch",])
ResMat.eq0.rho0 <- apply(X=df2.PiArm0, MARGIN=1,
                         FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                      pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                      rho00=icc.ctrlonly, rho10=icc.ctrlonly, rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                                      n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                      n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                      K.tot=K.total)))
df2.PiArm0$exchvar.eq.rho0 <- unlist(ResMat.eq0.rho0["Exch",])
ResMat.eq0.rho1 <- apply(X=df2.PiArm0, MARGIN=1,
                         FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                                      pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                                      rho00=icc.trt, rho10=icc.trt, rho01=icc.trt, rho11=icc.trt,
                                                      n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                                      n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                                      K.tot=K.total)))
df2.PiArm0$exchvar.eq.rho1 <- unlist(ResMat.eq0.rho1["Exch",])

Plot.Res4(resdf1=df2.PiArm0, resdf2=df2.PiArm1,
          xlab1=bquote(pi["01"]/pi["00"]), xlab2=bquote(pi["11"]/pi["10"]),
          xvar1="pi.ratio", xvar2="pi.ratio",
          xtitle1=bquote("Ratio of"~pi["01"]~"to"~pi["00"]), xtitle2=bquote("Ratio of"~pi["11"]~"to"~pi["10"]),
          fileout="2_ArmsPiRatio",
          truth1=log(pi.ctrl.E1/pi.ctrl.E0), truth2=log(pi.trt.E1/pi.trt.E0),
          ylimit1=c(0.075,0.125), ylimit2=c(85,100),
          colors=Res.cols, logX1=TRUE, logX2=TRUE)


### Varying Mean Cluster Size among High- and Low-Electrification Clusters ###
nwts <- sapply(X=nvec, FUN=function(y) dgpois(x=y-1,lambda=muAdj, omega=1-sqrt(muAdj/var)))
df2.VaryE1 <- data.frame(Means=seq(from=5,to=20,length.out=len))
df2.VaryE1$Vars <- (df2.VaryE1$Means*cv)^2
nwtsVary <- apply(X=as.matrix(df2.VaryE1), MARGIN=1,
               FUN=function(x) sapply(X=nvec,FUN=function(y) dgpois(x=y-1,lambda=x["Means"]-1,
                                                                    omega=1-sqrt((x["Means"]-1)/x["Vars"]))))
nwtsVary <- ifelse(is.nan(nwtsVary),0,nwtsVary)
nwtsVary <- prop.table(nwtsVary,2)
df2.VaryE1$truemeans <- apply(nwtsVary, 2, FUN=function(x) sum(x*nvec))
df2.VaryE1$truevars <- apply(nwtsVary, 2, FUN=function(x) sum(x*(nvec-sum(x*nvec))^2))
df2.VaryE1$trueCVs <- round(sqrt(df2.VaryE1$truevars)/df2.VaryE1$truemeans,digits=3)
df2.VaryE1$truemeans <- round(df2.VaryE1$truemeans, digits=3)
df2.VaryE0 <- df2.VaryE1
df2.VaryE1$truetotmean <- .6*df2.VaryE1$truemeans+.4*(muAdj+1)
df2.VaryE0$truetotmean <- .6*(muAdj+1)+.4*df2.VaryE0$truemeans
ResMat <- apply(X=nwtsVary, MARGIN=2,
                FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                             pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                             rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                             n00wts=nwts*.2, n10wts=nwts*.2, n01wts=x*.3, n11wts=x*.3,
                                             K.tot=K.total)))
df2.VaryE1$corvar <- unlist(ResMat["Cor",])
df2.VaryE1$indvar <- unlist(ResMat["Ind",])
df2.VaryE1$exchvar <- unlist(ResMat["Exch",])
df2.VaryE1$rho.star <- unlist(ResMat["rho.star",])
ResMat0 <- apply(X=nwtsVary, MARGIN=2,
                FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                             pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                             rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                             n00wts=x*.2, n10wts=x*.2, n01wts=nwts*.3, n11wts=nwts*.3,
                                             K.tot=K.total)))
df2.VaryE0$corvar <- unlist(ResMat0["Cor",])
df2.VaryE0$indvar <- unlist(ResMat0["Ind",])
df2.VaryE0$exchvar <- unlist(ResMat0["Exch",])
df2.VaryE0$rho.star <- unlist(ResMat0["rho.star",])
ResMat.eq <- matrix(data=NA, nrow=4, ncol=len)
row.names(ResMat.eq) <- c("rho.star","Cor","Ind","Exch")
ResMat.eq.rho0 <- ResMat.eq
ResMat.eq.rho1 <- ResMat.eq
for (i in 1:len) {
  ResMat.eq[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                rho00=df2.VaryE1$rho.star[i], rho10=df2.VaryE1$rho.star[i],
                                rho01=df2.VaryE1$rho.star[i], rho11=df2.VaryE1$rho.star[i],
                                n00wts=nwts*.2, n10wts=nwts*.2,
                                n01wts=nwtsVary[,i]*.3, n11wts=nwtsVary[,i]*.3,
                                K.tot=K.total))
  ResMat.eq.rho0[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                     pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                     rho00=icc.ctrlonly, rho10=icc.ctrlonly, 
                                     rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                     n00wts=nwts*.2, n10wts=nwts*.2,
                                     n01wts=nwtsVary[,i]*.3, n11wts=nwtsVary[,i]*.3,
                                     K.tot=K.total))
  ResMat.eq.rho1[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                     pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                     rho00=icc.trt, rho10=icc.trt, 
                                     rho01=icc.trt, rho11=icc.trt,
                                     n00wts=nwts*.2, n10wts=nwts*.2,
                                     n01wts=nwtsVary[,i]*.3, n11wts=nwtsVary[,i]*.3,
                                     K.tot=K.total))
}
df2.VaryE1$corvar.eq <- ResMat.eq["Cor",]
df2.VaryE1$indvar.eq <- ResMat.eq["Ind",]
df2.VaryE1$exchvar.eq <- ResMat.eq["Exch",]
df2.VaryE1$exchvar.eq.rho0 <- ResMat.eq.rho0["Exch",]
df2.VaryE1$exchvar.eq.rho1 <- ResMat.eq.rho1["Exch",]

ResMat.eq0 <- matrix(data=NA, nrow=4, ncol=len)
row.names(ResMat.eq0) <- c("rho.star","Cor","Ind","Exch")
ResMat.eq0.rho0 <- ResMat.eq0
ResMat.eq0.rho1 <- ResMat.eq0
for (i in 1:len) {
  ResMat.eq0[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                 pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                 rho00=df2.VaryE0$rho.star[i], rho10=df2.VaryE0$rho.star[i],
                                 rho01=df2.VaryE0$rho.star[i], rho11=df2.VaryE0$rho.star[i],
                                 n00wts=nwtsVary[,i]*.2, n10wts=nwtsVary[,i]*.2,
                                 n01wts=nwts*.3, n11wts=nwts*.3,
                                 K.tot=K.total))
  ResMat.eq0.rho0[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                      pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                      rho00=icc.ctrlonly, rho10=icc.ctrlonly,
                                      rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                      n00wts=nwtsVary[,i]*.2, n10wts=nwtsVary[,i]*.2,
                                      n01wts=nwts*.3, n11wts=nwts*.3,
                                      K.tot=K.total))
  ResMat.eq0.rho1[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                      pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                      rho00=icc.trt, rho10=icc.trt,
                                      rho01=icc.trt, rho11=icc.trt,
                                      n00wts=nwtsVary[,i]*.2, n10wts=nwtsVary[,i]*.2,
                                      n01wts=nwts*.3, n11wts=nwts*.3,
                                      K.tot=K.total))
}
df2.VaryE0$corvar.eq <- ResMat.eq0["Cor",]
df2.VaryE0$indvar.eq <- ResMat.eq0["Ind",]
df2.VaryE0$exchvar.eq <- ResMat.eq0["Exch",]
df2.VaryE0$exchvar.eq.rho0 <- ResMat.eq0.rho0["Exch",]
df2.VaryE0$exchvar.eq.rho1 <- ResMat.eq0.rho1["Exch",]

Plot.Res4(resdf1=df2.VaryE0, resdf2=df2.VaryE1,
          xlab1=bquote("Mean Cluster Size when"~Z["2i"]==0), xlab2=bquote("Mean Cluster Size when"~Z["2i"]==1),
          xvar1="truemeans", xvar2="truemeans",
          xtitle1=bquote("Mean Cluster Size when"~Z["2i"]==0), xtitle2=bquote("Mean Cluster Size when"~Z["2i"]==1),
          fileout="2_EsSizeMeans",
          truth1=muAdj+1, truth2=muAdj+1,
          ylimit1=c(0.075,0.125), ylimit2=c(85,100),
          colors=Res.cols)


### Varying Cluster Size CV among High- and Low-Electrification Clusters ###
df2.VaryE1CV <- data.frame(CVs=seq(from=.1,to=.7,length.out=len))
df2.VaryE1CV$Vars <- (df2.VaryE1CV$CVs^2)*(muAdj+1)^2
nwtsVaryCV <- sapply(X=df2.VaryE1CV$Vars,
               FUN=function(x) sapply(X=nvec,
                                      FUN=function(y) dgpois(x=y-1,lambda=muAdj, omega=1-sqrt(muAdj/x))))
nwtsVaryCV <- ifelse(is.nan(nwtsVaryCV),0,nwtsVaryCV)
nwtsVaryCV <- prop.table(nwtsVaryCV,2)
df2.VaryE1CV$truemeans <- apply(nwtsVaryCV, 2, FUN=function(x) sum(x*nvec))
df2.VaryE1CV$truevars <- apply(nwtsVaryCV, 2, FUN=function(x) sum(x*(nvec-sum(x*nvec))^2))
df2.VaryE1CV$trueCVs <- round(sqrt(df2.VaryE1CV$truevars)/df2.VaryE1CV$truemeans,digits=3)
df2.VaryE0CV <- df2.VaryE1CV
ResMat <- apply(X=nwtsVaryCV, MARGIN=2,
                 FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                              pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=nwts*.2, n10wts=nwts*.2, n01wts=x*.3, n11wts=x*.3,
                                              K.tot=K.total)))
df2.VaryE1CV$corvar <- unlist(ResMat["Cor",])
df2.VaryE1CV$indvar <- unlist(ResMat["Ind",])
df2.VaryE1CV$exchvar <- unlist(ResMat["Exch",])
df2.VaryE1CV$rho.star <- unlist(ResMat["rho.star",])
ResMat0 <- apply(X=nwtsVaryCV, MARGIN=2,
                 FUN=function(x) unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                              pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=x*.2, n10wts=x*.2, n01wts=nwts*.3, n11wts=nwts*.3,
                                              K.tot=K.total)))
df2.VaryE0CV$corvar <- unlist(ResMat0["Cor",])
df2.VaryE0CV$indvar <- unlist(ResMat0["Ind",])
df2.VaryE0CV$exchvar <- unlist(ResMat0["Exch",])
df2.VaryE0CV$rho.star <- unlist(ResMat0["rho.star",])
ResMat.eq <- matrix(data=NA, nrow=4, ncol=len)
row.names(ResMat.eq) <- c("rho.star","Cor","Ind","Exch")
ResMat.eq.rho0 <- ResMat.eq
ResMat.eq.rho1 <- ResMat.eq
for (i in 1:len) {
  ResMat.eq[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                rho00=df2.VaryE1CV$rho.star[i], rho10=df2.VaryE1CV$rho.star[i],
                                rho01=df2.VaryE1CV$rho.star[i], rho11=df2.VaryE1CV$rho.star[i],
                                n00wts=nwts*.2, n10wts=nwts*.2, n01wts=nwtsVaryCV[,i]*.3, n11wts=nwtsVaryCV[,i]*.3,
                                K.tot=K.total))
  ResMat.eq.rho0[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                     pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                     rho00=icc.ctrlonly, rho10=icc.ctrlonly,
                                     rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                     n00wts=nwts*.2, n10wts=nwts*.2, n01wts=nwtsVaryCV[,i]*.3, n11wts=nwtsVaryCV[,i]*.3,
                                     K.tot=K.total))
  ResMat.eq.rho1[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                     pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                     rho00=icc.trt, rho10=icc.trt,
                                     rho01=icc.trt, rho11=icc.trt,
                                     n00wts=nwts*.2, n10wts=nwts*.2, n01wts=nwtsVaryCV[,i]*.3, n11wts=nwtsVaryCV[,i]*.3,
                                     K.tot=K.total))
}
df2.VaryE1CV$corvar.eq <- ResMat.eq["Cor",]
df2.VaryE1CV$indvar.eq <- ResMat.eq["Ind",]
df2.VaryE1CV$exchvar.eq <- ResMat.eq["Exch",]
df2.VaryE1CV$exchvar.eq.rho0 <- ResMat.eq.rho0["Exch",]
df2.VaryE1CV$exchvar.eq.rho1 <- ResMat.eq.rho1["Exch",]

ResMat.eq0 <- matrix(data=NA, nrow=4, ncol=len)
row.names(ResMat.eq0) <- c("rho.star","Cor","Ind","Exch")
ResMat.eq0.rho0 <- ResMat.eq0
ResMat.eq0.rho1 <- ResMat.eq0
for (i in 1:len) {
  ResMat.eq0[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                 pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                 rho00=df2.VaryE0CV$rho.star[i], rho10=df2.VaryE0CV$rho.star[i],
                                 rho01=df2.VaryE0CV$rho.star[i], rho11=df2.VaryE0CV$rho.star[i],
                                 n00wts=nwtsVaryCV[,i]*.2, n10wts=nwtsVaryCV[,i]*.2,
                                 n01wts=nwts*.3, n11wts=nwts*.3,
                                 K.tot=K.total))
  ResMat.eq0.rho0[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                      pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                      rho00=icc.ctrlonly, rho10=icc.ctrlonly,
                                      rho01=icc.ctrlonly, rho11=icc.ctrlonly,
                                      n00wts=nwtsVaryCV[,i]*.2, n10wts=nwtsVaryCV[,i]*.2,
                                      n01wts=nwts*.3, n11wts=nwts*.3,
                                      K.tot=K.total))
  ResMat.eq0.rho1[,i] <- unlist(Vars2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                      pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                      rho00=icc.trt, rho10=icc.trt,
                                      rho01=icc.trt, rho11=icc.trt,
                                      n00wts=nwtsVaryCV[,i]*.2, n10wts=nwtsVaryCV[,i]*.2,
                                      n01wts=nwts*.3, n11wts=nwts*.3,
                                      K.tot=K.total))
}
df2.VaryE0CV$corvar.eq <- ResMat.eq0["Cor",]
df2.VaryE0CV$indvar.eq <- ResMat.eq0["Ind",]
df2.VaryE0CV$exchvar.eq <- ResMat.eq0["Exch",]
df2.VaryE0CV$exchvar.eq.rho0 <- ResMat.eq0.rho0["Exch",]
df2.VaryE0CV$exchvar.eq.rho1 <- ResMat.eq0.rho1["Exch",]

Plot.Res4(resdf1=df2.VaryE0CV, resdf2=df2.VaryE1CV,
          xlab1=bquote("CV of Cluster Size Distribution when"~Z["1i"]==0),
          xlab2=bquote("CV of Cluster Size Distribution when"~Z["1i"]==1),
          xvar1="trueCVs", xvar2="trueCVs",
          xtitle1=bquote("CV of Cluster Size Distribution when"~Z["1i"]==0),
          xtitle2=bquote("CV of Cluster Size Distribution when"~Z["1i"]==1),
          fileout="2_EsSizeCVs",
          truth1=cv, truth2=cv,
          ylimit1=c(0.075,0.125), ylimit2=c(85,100),
          colors=Res.cols)

## Sample Size Plotting:
#Varying the Rho & Pi Ratios Due to Electrification Parameters Differentially Between Treatment Arms:
setEPS()
postscript(file=paste0(figfolder,"/SS_2CovarArmsRatios.eps"),
           width=12, height=12, paper="special")
par(mfrow=c(2,2))
Plot.SS.Full(resdf=df2.RhoArm0,bquote(rho["01"]/rho["00"]), xvar="rho.ratio",
             title=bquote("a."~"Ratio of"~rho["01"]~"to"~rho["00"]),
             truth=log(icc.ctrl.E1/icc.ctrl.E0),
             ylimit=c(80,120), colors=SS.cols, logX=TRUE)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.RhoArm1,bquote(rho["11"]/rho["10"]), xvar="rho.ratio",
             title=bquote("b."~"Ratio of"~rho["11"]~"to"~rho["10"]),
             truth=log(icc.trt.E1/icc.trt.E0),
             ylimit=c(80,120), colors=SS.cols, logX=TRUE)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.PiArm0,bquote(pi["01"]/pi["00"]), xvar="pi.ratio",
             title=bquote("c."~"Ratio of"~pi["01"]~"to"~pi["00"]),
             truth=log(pi.ctrl.E1/pi.ctrl.E0),
             ylimit=c(80,120), colors=SS.cols, logX=TRUE)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.PiArm1,bquote(pi["11"]/pi["10"]), xvar="pi.ratio",
             title=bquote("d."~"Ratio of"~pi["11"]~"to"~pi["10"]),
             truth=log(pi.trt.E1/pi.trt.E0),
             ylimit=c(80,120), colors=SS.cols, logX=TRUE)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
dev.off()

#Varying the Size Parameters Differentially Between Low- and High-Electrification Clusters:
setEPS()
postscript(file=paste0(figfolder,"/SS_2CovarEsSize.eps"),
           width=12, height=12, paper="special")
par(mfrow=c(2,2))
Plot.SS.Full(resdf=df2.VaryE0,xlab=bquote("Mean Cluster Size when"~Z["2i"]==0),
             xvar="truemeans", title=bquote("a."~"Mean Cluster Size when"~Z["2i"]==0),
             truth=muAdj+1, ylimit=c(80,120), colors=SS.cols)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.VaryE1,xlab=bquote("Mean Cluster Size when"~Z["2i"]==1),
             xvar="truemeans", title=bquote("b."~"Mean Cluster Size when"~Z["2i"]==1),
             truth=muAdj+1, ylimit=c(80,120), colors=SS.cols)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.VaryE0CV,xlab=bquote("CV of Cluster Size Distribution when"~Z["2i"]==0),
             xvar="trueCVs", title=bquote("c."~"CV of Cluster Size Distribution when"~Z["2i"]==0),
             truth=cv, ylimit=c(80,120), colors=SS.cols)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
Plot.SS.Full(resdf=df2.VaryE1CV,xlab=bquote("CV of Cluster Size Distribution when"~Z["2i"]==1),
             xvar="trueCVs", title=bquote("d."~"CV of Cluster Size Distribution when"~Z["2i"]==1),
             truth=cv, ylimit=c(80,120), colors=SS.cols)
legend(x="bottom",legend=c(expression("Correctly Specified"),
                           expression("Common Exch."~rho[0]),
                           expression("Independence"),
                           expression("Common Exch."~rho[1]),
                           expression("Common Exch."~rho*"*")),
       col=SS.cols[c(1,4,3,5,2)], lty=c(1,5,2,6,3), lwd=rep(4,3), cex=.95, bg="white",
       ncol=3)
dev.off()


