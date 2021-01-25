
###### Read in Variance Calculation Functions ######
source("Variance_Functions.R")

###### Load Required Packages ######
require(gee)
require(LaplacesDemon)
require(tidyverse)
require(patchwork)

###### Read Data ######
Lin  <- read.csv(file="Lin2018.csv")
Lin2 <- Lin[Lin$delta_prot=="Observed",] # Only use those with observed results

## Create cluster-level electricification variable: ##
Lin2$ClustElec <- sapply(Lin2$clusterid, FUN=function(x) mean(Lin2$elec[Lin2$clusterid==x]=="Has electricity"))
Lin2$ElecBin <- ifelse(Lin2$ClustElec>.5,1,0)

## Sumarries of cluster sizes ##
NumObs <- length(Lin2$clusterid)
Clustids <- unique(Lin2$clusterid)
NumClusts <- length(Clustids)
ClustSizes <- sapply(X=Clustids, FUN=function(x) length(Lin2$personid[Lin2$clusterid==x]))
muAdj <- mean(ClustSizes) - 1 # One less than mean, used for Generalized Poisson Approximation
var <- (sd(ClustSizes))^2
cv <- sqrt(var)/(muAdj+1)
nvec <- 1:50 # Potential Cluster Sizes used for Generalized Poisson Approximation

## Assess generalized Poisson approximation of cluster size distribution: ##
hist(ClustSizes, breaks=seq(from=.5,to=20.5,by=.5))
points(x=seq(.75,19.75,by=1), y=719*dgpois(seq(0,19,by=1), lambda=muAdj, omega=1-sqrt(muAdj/var)))

## Cluster sizes by electrification: ##
ClustSizes.E1 <- sapply(X=unique(Lin2$clusterid[Lin2$ElecBin==1]), FUN=function(x) length(Lin2$personid[Lin2$clusterid==x]))
ClustSizes.E0 <- sapply(X=unique(Lin2$clusterid[Lin2$ElecBin==0]), FUN=function(x) length(Lin2$personid[Lin2$clusterid==x]))

## Data sets for the control and treatment conditions used in our analysis: ##
CtrlOnly <- Lin2[Lin2$tr=="Control",]
Trt <- Lin2[Lin2$tr %in% c("WSH","Nutrition + WSH"),]
K.total <- length(unique(c(CtrlOnly$clusterid,Trt$clusterid)))

## Overall adjusted GEE model: ##
Full <- rbind(CtrlOnly,Trt)
Full$arm <- ifelse(Full$tr=="Control",0,1)
fullgee <- gee(posgi~arm+ElecBin, id=clusterid, data=Full, family=binomial, corstr="exchangeable")

## Fitting unadjusted GEE models within each stratum to get parameters: ##
Trt.gee <- gee(posgi~1, id=clusterid, data=Trt, family=binomial, corstr="exchangeable")
CtrlOnly.gee <- gee(posgi~1, id=clusterid, data=CtrlOnly, family=binomial, corstr="exchangeable")
Trt.E1.gee <- gee(posgi~1, id=clusterid, data=Trt, subset=ElecBin==1, family=binomial, corstr="exchangeable")
Trt.E0.gee <- gee(posgi~1, id=clusterid, data=Trt, subset=ElecBin==0, family=binomial, corstr="exchangeable")
Ctrl.E1.gee <- gee(posgi~1, id=clusterid, data=CtrlOnly, subset=ElecBin==1, family=binomial, corstr="exchangeable")
Ctrl.E0.gee <- gee(posgi~1, id=clusterid, data=CtrlOnly, subset=ElecBin==0, family=binomial, corstr="exchangeable")

## Observed ICCs for each stratum: ##
icc.ctrlonly <- CtrlOnly.gee$working.correlation[1,2]
icc.trt <- Trt.gee$working.correlation[1,2]
icc.trt.E1 <- Trt.E1.gee$working.correlation[1,2]
icc.trt.E0 <- Trt.E0.gee$working.correlation[1,2]
icc.ctrl.E1 <- Ctrl.E1.gee$working.correlation[1,2]
icc.ctrl.E0 <- Ctrl.E0.gee$working.correlation[1,2]
icc.Eratio <- icc.trt.E1/icc.trt.E0

## Observed prevalence rates for each stratum: ##
pi.ctrlonly <- unname(expit(CtrlOnly.gee$coefficients))
pi.trt <- unname(expit(Trt.gee$coefficients))
pi.trt.E1 <- unname(expit(Trt.E1.gee$coefficients))
pi.trt.E0 <- unname(expit(Trt.E0.gee$coefficients))
pi.ctrl.E1 <- unname(expit(Ctrl.E1.gee$coefficients))
pi.ctrl.E0 <- unname(expit(Ctrl.E0.gee$coefficients))
pi.Eratio <- pi.trt.E1/pi.trt.E0

## Number of clusters by stratum: ##
length(unique(CtrlOnly$clusterid))
length(unique(Trt$clusterid))
length(unique(CtrlOnly$clusterid[CtrlOnly$ElecBin==0]))
length(unique(Trt$clusterid[Trt$ElecBin==0]))
length(unique(CtrlOnly$clusterid[CtrlOnly$ElecBin==1]))
length(unique(Trt$clusterid[Trt$ElecBin==1]))


###### Plotting Functions ######
len <- 401 #Length of vectors to generate plot lines (higher means more resolution)

Res.cols <- c(1,"#1f78b4","#a6cee3")
SS.cols <- c(1,"#bae4b3","#756bb1","#54278f")

Plot.Res4.gg <- function(resdf1,resdf2,xlab1,xlab2,
                         xtitle1,xtitle2,
                         xvar1,xvar2,
                         fileout,
                         xbreaks1=NULL, xbreaks2=NULL,
                         xlabels1=NULL, xlabels2=NULL,
                         ylimSE=NULL, ylimARE=NULL,
                         truth1=NULL, truth2=NULL,
                         truecol="red", truelty=2, linesize=1,
                         truth1b=NULL, truth2b=NULL, trueltyb=2) {
  
  tib1 <- tibble(X=rep(resdf1[,xvar1],3),
                     Type=rep(c("(i) Correct","(ii) Independent","(iii) Exchangeable"), 
                              each=length(resdf1[,xvar1])),
                     SE=sqrt(c(resdf1$corvar,resdf1$indvar,resdf1$exchvar)),
                     ARE=c(resdf1$corvar/resdf1$corvar*100,
                           resdf1$corvar/resdf1$indvar*100,
                           resdf1$corvar/resdf1$exchvar*100))
  
  tib2 <- tibble(X=rep(resdf2[,xvar2],3),
                 Type=rep(c("(i) Correct","(ii) Independent","(iii) Exchangeable"), 
                          each=length(resdf2[,xvar2])),
                 SE=sqrt(c(resdf2$corvar,resdf2$indvar,resdf2$exchvar)),
                 ARE=c(resdf2$corvar/resdf2$corvar*100,
                       resdf2$corvar/resdf2$indvar*100,
                       resdf2$corvar/resdf2$exchvar*100))
  
  g.a <- ggplot(tib1, aes(x=X, y=SE, color=Type, linetype=Type)) +
    geom_line(size=linesize) +
    scale_color_manual(name="Working Correlation Structure",
                       values=Res.cols) +
    scale_linetype_manual(name="Working Correlation Structure",
                          values=c(1,4,3)) +
    guides(color=guide_legend(override.aes=list(size=1))) +
    theme_light()
  if (!is.null(truth1)) {
    g.a <- g.a + geom_vline(aes(xintercept=truth1), color=truecol, linetype=truelty)
  }
  if (!is.null(truth1b)) {
    g.a <- g.a + geom_vline(aes(xintercept=truth1b), color=truecol, linetype=trueltyb)
  }
  if (is.null(xbreaks1)) {
    g.a <- g.a + scale_x_continuous(name=xlab1)
  } else if (is.null(xlabels1)) {
    g.a <- g.a + scale_x_continuous(name=xlab1,
                                    breaks=xbreaks1,
                                    minor_breaks=NULL,
                                    expand=expansion(mult=0))
  } else {
    g.a <- g.a + scale_x_continuous(name=xlab1,
                                    breaks=xbreaks1,
                                    minor_breaks=NULL,
                                    labels=xlabels1,
                                    expand=expansion(mult=0))
  }
  if (is.null(ylimSE)) {
    g.a <- g.a + scale_y_continuous(name="Standard Error of Estimate")
  } else {
    g.a <- g.a + scale_y_continuous(name="Standard Error of Estimate",
                                    limits=ylimSE)
  }
  if (!is.null(xtitle1)) {
    g.a <- g.a + labs(subtitle=bquote("SE"~"vs."~.(xtitle1)),
                      tag="A") +
      theme(plot.subtitle=element_text(hjust=0.5))
  }
  
  g.c <- ggplot(tib1, aes(x=X, y=ARE, color=Type, linetype=Type)) +
    geom_line(size=linesize) +
    scale_color_manual(name="Working Correlation Structure",
                       values=Res.cols) +
    scale_linetype_manual(name="Working Correlation Structure",
                          values=c(1,4,3)) +
    guides(color=guide_legend(override.aes=list(size=1))) +
    theme_light()
  if (!is.null(truth1)) {
    g.c <- g.c + geom_vline(aes(xintercept=truth1), color=truecol, linetype=truelty)
  }
  if (!is.null(truth1b)) {
    g.c <- g.c + geom_vline(aes(xintercept=truth1b), color=truecol, linetype=trueltyb)
  }
  if (is.null(xbreaks1)) {
    g.c <- g.c + scale_x_continuous(name=xlab1)
  } else if (is.null(xlabels1)) {
    g.c <- g.c + scale_x_continuous(name=xlab1,
                                    breaks=xbreaks1,
                                    minor_breaks=NULL,
                                    expand=expansion(mult=0))
  } else {
    g.c <- g.c + scale_x_continuous(name=xlab1,
                                    breaks=xbreaks1,
                                    minor_breaks=NULL,
                                    labels=xlabels1,
                                    expand=expansion(mult=0))
  }
  if (is.null(ylimSE)) {
    g.c <- g.c + scale_y_continuous(name="ARE (%)")
  } else {
    g.c <- g.c + scale_y_continuous(name="ARE (%)",
                                    limits=ylimARE)
  }
  if (!is.null(xtitle1)) {
    g.c <- g.c + labs(subtitle=bquote("ARE"~"vs."~.(xtitle1)),
                      tag="C") +
      theme(plot.subtitle=element_text(hjust=0.5))
  }
  
  g.b <- ggplot(tib2, aes(x=X, y=SE, color=Type, linetype=Type)) +
    geom_line(size=linesize) +
    scale_color_manual(name="Working Correlation Structure",
                       values=Res.cols) +
    scale_linetype_manual(name="Working Correlation Structure",
                          values=c(1,4,3)) +
    guides(color=guide_legend(override.aes=list(size=1))) +
    theme_light()
  if (!is.null(truth2)) {
    g.b <- g.b + geom_vline(aes(xintercept=truth2), color=truecol, linetype=truelty)
  }
  if (!is.null(truth2b)) {
    g.b <- g.b + geom_vline(aes(xintercept=truth2b), color=truecol, linetype=trueltyb)
  }
  if (is.null(xbreaks2)) {
    g.b <- g.b + scale_x_continuous(name=xlab2)
  } else if (is.null(xlabels2)) {
    g.b <- g.b + scale_x_continuous(name=xlab2,
                                    breaks=xbreaks2,
                                    minor_breaks=NULL,
                                    expand=expansion(mult=0))
  } else {
    g.b <- g.b + scale_x_continuous(name=xlab2,
                                    breaks=xbreaks2,
                                    minor_breaks=NULL,
                                    labels=xlabels2,
                                    expand=expansion(mult=0))
  }
  if (is.null(ylimSE)) {
    g.b <- g.b + scale_y_continuous(name="Standard Error of Estimate")
  } else {
    g.b <- g.b + scale_y_continuous(name="Standard Error of Estimate",
                                    limits=ylimSE)
  }
  if (!is.null(xtitle2)) {
    g.b <- g.b + labs(subtitle=bquote("SE"~"vs."~.(xtitle2)),
                      tag="B") +
      theme(plot.subtitle=element_text(hjust=0.5))
  }
  
  g.d <- ggplot(tib2, aes(x=X, y=ARE, color=Type, linetype=Type)) +
    geom_line(size=linesize) +
    scale_color_manual(name="Working Correlation Structure",
                       values=Res.cols) +
    scale_linetype_manual(name="Working Correlation Structure",
                          values=c(1,4,3)) +
    guides(color=guide_legend(override.aes=list(size=1))) +
    theme_light()
  if (!is.null(truth2)) {
    g.d <- g.d + geom_vline(aes(xintercept=truth2), color=truecol, linetype=truelty)
  }
  if (!is.null(truth2b)) {
    g.d <- g.d + geom_vline(aes(xintercept=truth2b), color=truecol, linetype=trueltyb)
  }
  if (is.null(xbreaks2)) {
    g.d <- g.d + scale_x_continuous(name=xlab2)
  } else if (is.null(xlabels2)) {
    g.d <- g.d + scale_x_continuous(name=xlab2,
                                    breaks=xbreaks2,
                                    minor_breaks=NULL,
                                    expand=expansion(mult=0))
  } else {
    g.d <- g.d + scale_x_continuous(name=xlab2,
                                    breaks=xbreaks2,
                                    minor_breaks=NULL,
                                    labels=xlabels2,
                                    expand=expansion(mult=0))
  }
  if (is.null(ylimSE)) {
    g.d <- g.d + scale_y_continuous(name="ARE (%)")
  } else {
    g.d <- g.d + scale_y_continuous(name="ARE (%)",
                                    limits=ylimARE)
  }
  if (!is.null(xtitle2)) {
    g.d <- g.d + labs(subtitle=bquote("ARE"~"vs."~.(xtitle2)),
                      tag="D") +
      theme(plot.subtitle=element_text(hjust=0.5))
  }
  
  ggsave(filename=fileout, 
         plot=g.a+g.b+g.c+g.d+
           plot_layout(ncol=2,nrow=2,byrow=TRUE, guides="collect") & theme(legend.position = "bottom"),
         device="eps", width=7, height=7, units="in")
}

Plot.SS.single <- function(resdf,xlab,xvar,
                           title=NULL,truth=NULL,
                           ylimit=NULL,
                           xbreaks=NULL, xlabels=NULL,
                           truecol="red", truelty=2,
                           linesize=1,
                           truthb=NULL, trueltyb=2) {
  tib <- tibble(X=rep(resdf[,xvar],4),
                Type=factor(rep(c("Corr","rhostar","rho0","rho1"),each=length(resdf[,xvar])),
                            levels=c("Corr","rhostar","rho0","rho1")),
                SS.ratio=c(resdf$SSratio.cor,resdf$SSratio.exch.rhost,
                           resdf$SSratio.exch.rho0,resdf$SSratio.exch.rho1)*100)
  g <- ggplot(tib,
              aes(x=X, y=SS.ratio, color=Type, linetype=Type)) +
    geom_line(size=linesize) +
    scale_color_manual(name="Correlation Specifications",
                       values=SS.cols,
                       labels=c(expression("Correctly Specified"),
                                expression("Common Exchangeable,"~rho*"*"),
                                expression("Common Exchangeable,"~rho[0]),
                                expression("Common Exchangeable,"~rho[1]))) +
    scale_linetype_manual(name="Correlation Specifications",
                          values=c(1,4,3,3),
                          labels=c(expression("Correctly Specified"),
                                   expression("Common Exchangeable,"~rho*"*"),
                                   expression("Common Exchangeable,"~rho[0]),
                                   expression("Common Exchangeable,"~rho[1]))) +
    theme_light()
  if (!is.null(truth)) {
    g <- g + geom_vline(aes(xintercept=truth), color=truecol, linetype=truelty)
  }
  if (!is.null(truthb)) {
    g <- g + geom_vline(aes(xintercept=truthb), color=truecol, linetype=trueltyb)
  }
  if (is.null(xbreaks)) {
    g <- g + scale_x_continuous(name=xlab)
  } else if (is.null(xlabels)) {
    g <- g + scale_x_continuous(name=xlab,
                                breaks=xbreaks,
                                minor_breaks=NULL,
                                expand=expansion(mult=0))
  } else {
    g <- g + scale_x_continuous(name=xlab,
                                breaks=xbreaks,
                                minor_breaks=NULL,
                                labels=xlabels,
                                expand=expansion(mult=0))
  }
  if (is.null(ylimit)) {
    g <- g + scale_y_continuous(name="Estimated SS/Required SS (%)")
  } else {
    g <- g + scale_y_continuous(name="Estimated SS/Required SS (%)",
                                    limits=ylimit)
  }
  if (!is.null(title)) {
    g <- g + labs(subtitle=title) +
      theme(plot.subtitle=element_text(hjust=0.5))
  }
  return(g)
}

#### Sample Size Calculation Parameters ####
alpha <- 0.05
power <- 0.80

#### Names of sample size columns ####
SSmat.names <- c("SSratio.cor","SSratio.ind","SSratio.exch.rhost",
                 "SSratio.exch.rho0","SSratio.exch.rho1")

###### Main Text Figures ######

### Figure 1 ###

## Varying ICC ##
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
df.VaryRho1$rho.star <- unlist(ResMat["rho.star",])
ResMat.SS <- sapply(X=df.VaryRho1$rho1s,
                    FUN=function(x) unlist(SSest(n0vec=ClustSizes,n1vec=ClustSizes,
                                                 rho0=icc.ctrlonly,rho1=x,
                                                 pi0=pi.ctrlonly,pi1=pi.trt,
                                                 n0wts=NULL,n1wts=NULL,
                                                 sig=alpha,pwr=power)))
rownames(ResMat.SS) <- SSmat.names 
df.VaryRho1 <- cbind(df.VaryRho1,t(ResMat.SS))

## Varying Outcome Probability ##
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
df.VaryPi1$rho.star <- unlist(ResMat["rho.star",])
ResMat.SS <- sapply(X=df.VaryPi1$pi1s,
                    FUN=function(x) unlist(SSest(n0vec=ClustSizes,n1vec=ClustSizes,
                                                 rho0=icc.ctrlonly,rho1=icc.trt,
                                                 pi0=pi.ctrlonly,pi1=x,
                                                 n0wts=NULL,n1wts=NULL,
                                                 sig=alpha,pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df.VaryPi1 <- cbind(df.VaryPi1,t(ResMat.SS))

## Plotting Figure 1 ##
Plot.Res4.gg(resdf1=df.VaryRho1, resdf2=df.VaryPi1,
             xlab1=bquote(rho["1"]/rho["0"]), xlab2=bquote(pi["1"]/pi["0"]),
             xtitle1=bquote("Ratio of"~rho["1"]~"to"~rho["0"]), xtitle2=bquote("Ratio of"~pi["1"]~"to"~pi["0"]),
             xvar1="rho.ratio", xvar2="pi.ratio",
             fileout="Figure1.eps",
             xbreaks1=log(c(1/4,1/2,3/4,1,3/2,2,3,4)), 
             xbreaks2=log(c(1/4,1/2,3/4,1,3/2,2,3,4)),
             xlabels1=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"), 
             xlabels2=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=log(icc.trt/icc.ctrlonly), truth2=log(pi.trt/pi.ctrlonly),
             truecol="red", truelty=2, linesize=1)


### Figure 2 ###

## Varying Cluster Size CVs ##
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
df.VaryCVs$rho.star <- unlist(ResMat["rho.star",])
ResMat.SS <- apply(X=nwts, MARGIN=2,
                    FUN=function(x) unlist(SSest(n0vec=nvec,n1vec=nvec,
                                                 rho0=icc.ctrlonly,rho1=icc.trt,
                                                 pi0=pi.ctrlonly,pi1=pi.trt,
                                                 n0wts=x,n1wts=x,
                                                 sig=alpha,pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df.VaryCVs <- cbind(df.VaryCVs,t(ResMat.SS))

## Varying Cluster Size Means ##
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
df.VaryMeans$rho.star <- unlist(ResMat["rho.star",])
ResMat.SS <- apply(X=nwts, MARGIN=2,
                   FUN=function(x) unlist(SSest(n0vec=nvec,n1vec=nvec,
                                                rho0=icc.ctrlonly,rho1=icc.trt,
                                                pi0=pi.ctrlonly,pi1=pi.trt,
                                                n0wts=x,n1wts=x,
                                                sig=alpha,pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df.VaryMeans <- cbind(df.VaryMeans,t(ResMat.SS))

## Plotting Figure 2 ##
Plot.Res4.gg(resdf1=df.VaryMeans, resdf2=df.VaryCVs,
             xlab1="Mean Cluster Size", xlab2="CV of Cluster Sizes",
             xtitle1="Mean Cluster Size", xtitle2="CV of Cluster Sizes",
             xvar1="Means", xvar2="CVs",
             fileout="Figure2.eps",
             xbreaks1=seq(5,20,by=5), 
             xbreaks2=seq(0.1,0.7,by=0.1),
             xlabels1=NULL, 
             xlabels2=NULL,
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=mean(ClustSizes), truth2=sqrt(var)/(muAdj+1),
             truecol="red", truelty=2, linesize=1)


### Figure 3 ###

g.3.a <- Plot.SS.single(resdf=df.VaryRho1,
                        xlab=expression(rho[1]/rho[0]),
                        xvar="rho.ratio",
                        title=expression("Ratio of"~rho[1]~"to"~rho[0]),
                        truth=log(icc.trt/icc.ctrlonly),
                        ylimit=c(80,120),
                        xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                        xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                        truecol="red", truelty=2, linesize=1.5) +
  labs(tag="A")
g.3.b <- Plot.SS.single(resdf=df.VaryPi1,
                        xlab=expression(pi[1]/pi[0]),
                        xvar="pi.ratio",
                        title=expression("Ratio of"~pi[1]~"to"~pi[0]),
                        truth=log(pi.trt/pi.ctrlonly),
                        ylimit=c(80,120),
                        xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                        xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                        truecol="red", truelty=2, linesize=1.5) +
  labs(tag="B")
g.3.c <- Plot.SS.single(resdf=df.VaryMeans,
                        xlab=expression("Mean"~"Cluster Size"),
                        xvar="Means",
                        title=expression("Mean Cluster Size"),
                        truth=muAdj+1,
                        ylimit=c(80,120),
                        xbreaks=seq(5,20,by=5),
                        xlabels=NULL,
                        truecol="red", truelty=2, linesize=1.5) +
  labs(tag="C")
g.3.d <- Plot.SS.single(resdf=df.VaryCVs,
                        xlab=expression("CV"~"of Cluster Size"),
                        xvar="CVs",
                        title=expression("CV of Cluster Sizes"),
                        truth=sqrt(var)/(muAdj+1),
                        ylimit=c(80,120),
                        xbreaks=seq(0.1,0.7,by=0.1),
                        xlabels=NULL,
                        truecol="red", truelty=2, linesize=1.5) +
  labs(tag="D")
ggsave(filename="Figure3.eps",
       plot=g.3.a+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.3.b+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.3.c+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.3.d+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         plot_layout(ncol=2,nrow=2,byrow=TRUE, guides="collect") & theme(legend.position = "bottom"),
       device="eps", width=7, height=7, units="in")


##### Supplemental Figures #####

###### Two Binary Cluster-Level Covariates ######

### Varying Rho Ratio ###
df2.RhoRatio <- data.frame(rho.ratio=seq(from=log(1/4),to=log(4),length.out=len))
df2.RhoRatio$rho00s <- icc.ctrl.E0
df2.RhoRatio$rho10s <- icc.trt.E0
df2.RhoRatio$rho01s <- df2.RhoRatio$rho00s*exp(df2.RhoRatio$rho.ratio)
df2.RhoRatio$rho11s <- df2.RhoRatio$rho10s*exp(df2.RhoRatio$rho.ratio)
ResMat <- apply(X=df2.RhoRatio, MARGIN=1,
                 FUN=function(x) unlist(Vars2(n00vec=ClustSizes,n01vec=ClustSizes, n10vec=ClustSizes, n11vec=ClustSizes,
                                       pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                       rho00=x["rho00s"], rho10=x["rho10s"], rho01=x["rho01s"], rho11=x["rho11s"],
                                       n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                       n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                       K.tot=K.total)))
df2.RhoRatio$corvar <- unlist(ResMat["Cor",])
df2.RhoRatio$indvar <- unlist(ResMat["Ind",])
df2.RhoRatio$exchvar <- unlist(ResMat["Exch",])
df2.RhoRatio$rho.star <- unlist(ResMat["rho.star",])
ResMat.SS <- apply(X=df2.RhoRatio, MARGIN=1,
                FUN=function(x) unlist(SSest2(n00vec=ClustSizes,n01vec=ClustSizes, n10vec=ClustSizes, n11vec=ClustSizes,
                                              pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                              rho00=x["rho00s"], rho10=x["rho10s"], rho01=x["rho01s"], rho11=x["rho11s"],
                                              n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                              n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                              sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.RhoRatio <- cbind(df2.RhoRatio, t(ResMat.SS))



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
ResMat.SS <- apply(X=df2.pi.Eratio, MARGIN=1,
                FUN=function(x) unlist(SSest2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                              pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                              n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                              sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.pi.Eratio <- cbind(df2.pi.Eratio, t(ResMat.SS))

Plot.Res4.gg(resdf1=df2.RhoRatio, resdf2=df2.pi.Eratio,
             xlab1=bquote(rho["01"]/rho["00"]==rho["11"]/rho["10"]),
             xlab2=bquote(pi["01"]/pi["00"]==pi["11"]/pi["10"]),
             xtitle1=bquote("Ratio of"~rho["01"]~"to"~rho["00"]~"and"~rho["11"]~"to"~rho["10"]), 
             xtitle2=bquote("Ratio of"~pi["01"]~"to"~pi["00"]~"and"~pi["11"]~"to"~pi["10"]),
             xvar1="rho.ratio", xvar2="pi.ratio",
             fileout="FigureS1.eps",
             xbreaks1=log(c(1/4,1/2,3/4,1,3/2,2,3,4)), 
             xbreaks2=log(c(1/4,1/2,3/4,1,3/2,2,3,4)),
             xlabels1=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"), 
             xlabels2=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=log(icc.Eratio), truth2=log(pi.Eratio),
             truecol="red", truelty=2, linesize=1, 
             truth1b=log(icc.ctrl.E1/icc.ctrl.E0), truth2b=log(pi.ctrl.E1/pi.ctrl.E0), trueltyb=5)


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
df2.VaryCVs$rho.star <- rep(df2.VaryCVs.rho.star,dim(df2.VaryCVs)[1])
ResMat.SS <- apply(X=nwts, MARGIN=2,
                FUN=function(x) unlist(SSest2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                              pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                              sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.VaryCVs <- cbind(df2.VaryCVs, t(ResMat.SS))

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
df2.VaryMeans$rho.star <- rep(df2.VaryMeans.rho.star,dim(df2.VaryMeans)[1])
ResMat.SS <- apply(X=nwts, MARGIN=2,
                FUN=function(x) unlist(SSest2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                              pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=x*.2, n10wts=x*.2, n01wts=x*.3, n11wts=x*.3,
                                              sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.VaryMeans <- cbind(df2.VaryMeans,t(ResMat.SS))

Plot.Res4.gg(resdf1=df2.VaryMeans, resdf2=df2.VaryCVs,
             xlab1="Mean Cluster Size", xlab2="CV of Cluster Sizes",
             xtitle1="Mean Cluster Size", xtitle2="CV of Cluster Sizes",
             xvar1="truemeans", xvar2="trueCVs",
             fileout="FigureS2.eps",
             xbreaks1=seq(5,20,by=5), 
             xbreaks2=seq(0.1,0.7,by=0.1),
             xlabels1=NULL, 
             xlabels2=NULL,
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=mean(ClustSizes), sqrt(var)/(muAdj+1),
             truecol="red", truelty=2, linesize=1)


### Sample Size Plotting for Two Covariate Cases ###
g.s7.a <- Plot.SS.single(resdf=df2.RhoRatio,
                         xlab=bquote(rho["01"]/rho["00"]==rho["11"]/rho["10"]),
                         xvar="rho.ratio",
                         title=bquote("Ratio of"~rho["01"]~"to"~rho["00"]~"and"~rho["11"]~"to"~rho["10"]),
                         truth=log(icc.Eratio),
                         ylimit=c(80,120),
                         xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                         xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                         truecol="red", truelty=2, linesize=1.5,
                         truthb=log(icc.ctrl.E1/icc.ctrl.E0), trueltyb=5) +
  labs(tag="A")
g.s7.b <- Plot.SS.single(resdf=df2.pi.Eratio,
                         xlab=bquote(pi["01"]/pi["00"]==pi["11"]/pi["10"]),
                         xvar="pi.ratio",
                         title=bquote("Ratio of"~pi["01"]~"to"~pi["00"]~"and"~pi["11"]~"to"~pi["10"]),
                         truth=log(pi.Eratio),
                         ylimit=c(80,120),
                         xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                         xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                         truecol="red", truelty=2, linesize=1.5,
                         truthb=log(pi.ctrl.E1/pi.ctrl.E0), trueltyb=5) +
  labs(tag="B")
g.s7.c <- Plot.SS.single(resdf=df2.VaryMeans,
                         xlab=bquote("Mean"~"Cluster Size"),
                         xvar="truemeans",
                         title=bquote("Mean Cluster Size"),
                         truth=mean(ClustSizes),
                         ylimit=c(80,120),
                         xbreaks=seq(5,20,by=5),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="C")
g.s7.d <- Plot.SS.single(resdf=df2.VaryCVs,
                         xlab=bquote("CV"~"of Cluster Sizes"),
                         xvar="trueCVs",
                         title=bquote("CV of Cluster Sizes"),
                         truth=sqrt(var)/(muAdj+1),
                         ylimit=c(80,120),
                         xbreaks=seq(0.1,0.7,by=0.1),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="D")
ggsave(filename="FigureS7.eps",
       plot=g.s7.a+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s7.b+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s7.c+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s7.d+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         plot_layout(ncol=2,nrow=2,byrow=TRUE, guides="collect") & theme(legend.position = "bottom"),
       device="eps", width=7, height=7, units="in")


###### Varying Parameters Differentially Between Arms ######

### Varying Ratio of ICCs due to Electrification Differentially ###
###  Between Treated (Arm 1) and Control (Arm 0) Arms ###
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
ResMat.SS <- apply(X=df2.RhoArm1, MARGIN=1,
                FUN=function(x) unlist(SSest2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                             pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                             rho00=x["rho00s"], rho10=x["rho10s"], rho01=x["rho01s"], rho11=x["rho11s"],
                                             n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                             n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                             sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.RhoArm1 <- cbind(df2.RhoArm1, t(ResMat.SS))
ResMat.SS0 <- apply(X=df2.RhoArm0, MARGIN=1,
                 FUN=function(x) unlist(SSest2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                               pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                               rho00=x["rho00s"], rho10=x["rho10s"], rho01=x["rho01s"], rho11=x["rho11s"],
                                               n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                               n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                               sig=alpha, pwr=power)))
rownames(ResMat.SS0) <- SSmat.names
df2.RhoArm0 <- cbind(df2.RhoArm0, t(ResMat.SS0))

Plot.Res4.gg(resdf1=df2.RhoArm0, resdf2=df2.RhoArm1,
             xlab1=bquote(rho["01"]/rho["00"]), xlab2=bquote(rho["11"]/rho["10"]),
             xtitle1=bquote("Ratio of"~rho["01"]~"to"~rho["00"]), 
             xtitle2=bquote("Ratio of"~rho["11"]~"to"~rho["10"]),
             xvar1="rho.ratio", xvar2="rho.ratio",
             fileout="FigureS3.eps",
             xbreaks1=log(c(1/4,1/2,3/4,1,3/2,2,3,4)), 
             xbreaks2=log(c(1/4,1/2,3/4,1,3/2,2,3,4)),
             xlabels1=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"), 
             xlabels2=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=log(icc.ctrl.E1/icc.ctrl.E0), truth2=log(icc.trt.E1/icc.trt.E0),
             truecol="red", truelty=2, linesize=1)


### Varying Ratio of Pis due to Electrification Differentially ###
###  Between Treated (Arm 1) and Control (Arm 0) Arms ###
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
ResMat.SS <- apply(X=df2.PiArm1, MARGIN=1,
                FUN=function(x) unlist(SSest2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                              pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                              n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                              sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.PiArm1 <- cbind(df2.PiArm1, t(ResMat.SS))
ResMat.SS0 <- apply(X=df2.PiArm0, MARGIN=1,
                 FUN=function(x) unlist(SSest2(n00vec=ClustSizes,n10vec=ClustSizes, n01vec=ClustSizes, n11vec=ClustSizes,
                                               pi00=x["pi00s"], pi10=x["pi10s"], pi01=x["pi01s"], pi11=x["pi11s"],
                                               rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                               n00wts=rep(.2,length(ClustSizes)), n10wts=rep(.2,length(ClustSizes)),
                                               n01wts=rep(.3,length(ClustSizes)), n11wts=rep(.3,length(ClustSizes)),
                                               sig=alpha, pwr=power)))
rownames(ResMat.SS0) <- SSmat.names
df2.PiArm0 <- cbind(df2.PiArm0, t(ResMat.SS0))

Plot.Res4.gg(resdf1=df2.PiArm0, resdf2=df2.PiArm1,
             xlab1=bquote(pi["01"]/pi["00"]), xlab2=bquote(pi["11"]/pi["10"]),
             xtitle1=bquote("Ratio of"~pi["01"]~"to"~pi["00"]), 
             xtitle2=bquote("Ratio of"~pi["11"]~"to"~pi["10"]),
             xvar1="pi.ratio", xvar2="pi.ratio",
             fileout="FigureS4.eps",
             xbreaks1=log(c(1/4,1/2,3/4,1,3/2,2,3,4)), 
             xbreaks2=log(c(1/4,1/2,3/4,1,3/2,2,3,4)),
             xlabels1=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"), 
             xlabels2=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=log(pi.ctrl.E1/pi.ctrl.E0), truth2=log(pi.trt.E1/pi.trt.E0),
             truecol="red", truelty=2, linesize=1)


### Varying Mean Cluster Size Differentially ###
### Between High-Electrification (E1) and Low-Electrification (E0) Clusters ###
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

ResMat.SS <- apply(X=nwtsVary, MARGIN=2,
                   FUN=function(x) unlist(SSest2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                                 pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                                 rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                                 n00wts=nwts*.2, n10wts=nwts*.2, n01wts=x*.3, n11wts=x*.3,
                                                 sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.VaryE1 <- cbind(df2.VaryE1, t(ResMat.SS))
ResMat.SS0 <- apply(X=nwtsVary, MARGIN=2,
                 FUN=function(x) unlist(SSest2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                               pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                               rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                               n00wts=x*.2, n10wts=x*.2, n01wts=nwts*.3, n11wts=nwts*.3,
                                               sig=alpha, pwr=power)))
rownames(ResMat.SS0) <- SSmat.names
df2.VaryE0 <- cbind(df2.VaryE0, t(ResMat.SS0))

Plot.Res4.gg(resdf1=df2.VaryE0, resdf2=df2.VaryE1,
             xlab1=bquote("Mean Cluster Size when"~Z["2i"]==0), 
             xlab2=bquote("Mean Cluster Size when"~Z["2i"]==1),
             xtitle1=bquote("Mean Cluster Size when"~Z["2i"]==0), 
             xtitle2=bquote("Mean Cluster Size when"~Z["2i"]==1),
             xvar1="truemeans", xvar2="truemeans",
             fileout="FigureS5.eps",
             xbreaks1=seq(5,20,by=5), 
             xbreaks2=seq(5,20,by=5), 
             xlabels1=NULL, 
             xlabels2=NULL,
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=muAdj+1, truth2=muAdj+1,
             truecol="red", truelty=2, linesize=1)


### Varying Cluster Size CV Differentially ###
### Between High-Electrification (E1) and Low-Electrification (E0) Clusters ###
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

ResMat.SS <- apply(X=nwtsVaryCV, MARGIN=2,
                FUN=function(x) unlist(SSest2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                              pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                              rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                              n00wts=nwts*.2, n10wts=nwts*.2, n01wts=x*.3, n11wts=x*.3,
                                              sig=alpha, pwr=power)))
rownames(ResMat.SS) <- SSmat.names
df2.VaryE1CV <- cbind(df2.VaryE1CV, t(ResMat.SS))
ResMat.SS0 <- apply(X=nwtsVaryCV, MARGIN=2,
                 FUN=function(x) unlist(SSest2(n00vec=nvec, n10vec=nvec, n01vec=nvec, n11vec=nvec,
                                               pi00=pi.ctrl.E0, pi10=pi.trt.E0, pi01=pi.ctrl.E1, pi11=pi.trt.E1,
                                               rho00=icc.ctrl.E0, rho10=icc.trt.E0, rho01=icc.ctrl.E1, rho11=icc.trt.E1,
                                               n00wts=x*.2, n10wts=x*.2, n01wts=nwts*.3, n11wts=nwts*.3,
                                               sig=alpha, pwr=power)))
rownames(ResMat.SS0) <- SSmat.names
df2.VaryE0CV <- cbind(df2.VaryE0CV, t(ResMat.SS0))

Plot.Res4.gg(resdf1=df2.VaryE0CV, resdf2=df2.VaryE1CV,
             xlab1=bquote("CV of Cluster Sizes when"~Z["2i"]==0),
             xlab2=bquote("CV of Cluster Sizes when"~Z["2i"]==1),
             xtitle1=bquote("CV of Cluster Sizes when"~Z["2i"]==0),
             xtitle2=bquote("CV of Cluster Sizes when"~Z["2i"]==1),
             xvar1="trueCVs", xvar2="trueCVs",
             fileout="FigureS6.eps",
             xbreaks1=seq(0.1,0.7,by=0.1), 
             xbreaks2=seq(0.1,0.7,by=0.1),
             xlabels1=NULL, 
             xlabels2=NULL,
             ylimSE=c(0.075,0.125), ylimARE=c(90,100),
             truth1=cv, truth2=cv,
             truecol="red", truelty=2, linesize=1)

### Sample Size Plotting for Differentially Varying ###
### Rho and Pi Ratios Between Treatment Arms ###
g.s8.a <- Plot.SS.single(resdf=df2.RhoArm0,
                         xlab=bquote(rho["01"]/rho["00"]),
                         xvar="rho.ratio",
                         title=bquote("Ratio of"~rho["01"]~"to"~rho["00"]),
                         truth=log(icc.ctrl.E1/icc.ctrl.E0),
                         ylimit=c(80,120),
                         xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                         xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="A")
g.s8.b <- Plot.SS.single(resdf=df2.RhoArm1,
                         xlab=bquote(rho["11"]/rho["10"]),
                         xvar="rho.ratio",
                         title=bquote("Ratio of"~rho["11"]~"to"~rho["10"]),
                         truth=log(icc.trt.E1/icc.trt.E0),
                         ylimit=c(80,120),
                         xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                         xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="B")
g.s8.c <- Plot.SS.single(resdf=df2.PiArm0,
                         xlab=bquote(pi["01"]/pi["00"]),
                         xvar="pi.ratio",
                         title=bquote("Ratio of"~pi["01"]~"to"~pi["00"]),
                         truth=log(pi.ctrl.E1/pi.ctrl.E0),
                         ylimit=c(80,120),
                         xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                         xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="C")
g.s8.d <- Plot.SS.single(resdf=df2.PiArm1,
                         xlab=bquote(pi["11"]/pi["10"]),
                         xvar="pi.ratio",
                         title=bquote("Ratio of"~pi["11"]~"to"~pi["10"]),
                         truth=log(pi.trt.E1/pi.trt.E0),
                         ylimit=c(80,120),
                         xbreaks=log(c(0.25,0.5,0.75,1.0,1.5,2,3,4)),
                         xlabels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="D")
ggsave(filename="FigureS8.eps",
       plot=g.s8.a+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s8.b+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s8.c+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s8.d+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         plot_layout(ncol=2,nrow=2,byrow=TRUE, guides="collect") & theme(legend.position = "bottom"),
       device="eps", width=7, height=7, units="in")


### Sample Size Plotting for Differentially Varying ###
### Cluster Size Distribution Between High- and Low-Electrification Clusters ###
g.s9.a <- Plot.SS.single(resdf=df2.VaryE0,
                         xlab=bquote("Mean Cluster Size when"~Z["2i"]==0),
                         xvar="truemeans",
                         title=bquote("Mean Cluster Size when"~Z["2i"]==0),
                         truth=mean(ClustSizes),
                         ylimit=c(80,120),
                         xbreaks=seq(5,20,by=5),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="A")
g.s9.b <- Plot.SS.single(resdf=df2.VaryE1,
                         xlab=bquote("Mean Cluster Size when"~Z["2i"]==1),
                         xvar="truemeans",
                         title=bquote("Mean Cluster Size when"~Z["2i"]==1),
                         truth=mean(ClustSizes),
                         ylimit=c(80,120),
                         xbreaks=seq(5,20,by=5),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="B")
g.s9.c <- Plot.SS.single(resdf=df2.VaryE0CV,
                         xlab=bquote("CV of Cluster Sizes when"~Z["2i"]==0),
                         xvar="trueCVs",
                         title=bquote("CV of Cluster Sizes when"~Z["2i"]==0),
                         truth=sqrt(var)/(muAdj+1),
                         ylimit=c(80,120),
                         xbreaks=seq(0.1,0.7,by=0.1),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="C")
g.s9.d <- Plot.SS.single(resdf=df2.VaryE1CV,
                         xlab=bquote("CV of Cluster Sizes when"~Z["2i"]==1),
                         xvar="trueCVs",
                         title=bquote("CV of Cluster Sizes when"~Z["2i"]==1),
                         truth=sqrt(var)/(muAdj+1),
                         ylimit=c(80,120),
                         xbreaks=seq(0.1,0.7,by=0.1),
                         truecol="red", truelty=2, linesize=1.5) +
  labs(tag="D")
ggsave(filename="FigureS9.eps",
       plot=g.s9.a+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s9.b+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s9.c+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         g.s9.d+guides(color=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=1)))+theme(legend.position="bottom")+
         plot_layout(ncol=2,nrow=2,byrow=TRUE, guides="collect") & theme(legend.position = "bottom"),
       device="eps", width=7, height=7, units="in")


