#### One Binary Cluster-Level Covariate #####

VarFactors <- function(n0vec, n1vec, rho0, rho1, n0wts=NULL, n1wts=NULL) {
  if (is.null(n0wts) + is.null(n1wts) == 1) {
    print("Error: Specify Weights for both n0 and n1 or Neither")
  } else if (is.null(n0wts) & is.null(n1wts)) {
    n0wts <- rep(1,length(n0vec))
    n1wts <- rep(1,length(n1vec))
  }
  K0 <- sum(n0wts)
  K1 <- sum(n1wts)
  rhost <- (sum(n0wts*n0vec*(n0vec-1)/2)*rho0 + sum(n1wts*n1vec*(n1vec-1)/2)*rho1)/(sum(n0wts*n0vec*(n0vec-1)/2) + sum(n1wts*n1vec*(n1vec-1)/2))
  CorDenom <- c(sum(n0wts*n0vec/(1+(n0vec-1)*rho0))/K0,
                sum(n1wts*n1vec/(1+(n1vec-1)*rho1))/K1)
  n0mean <- sum(n0wts*n0vec)/K0
  n1mean <- sum(n1wts*n1vec)/K1
  n0var <- sum(n0wts*(n0vec-n0mean)^2)/K0
  n1var <- sum(n1wts*(n1vec-n1mean)^2)/K1
  IndVal <- c((1+((n0var/(n0mean^2)+1)*n0mean-1)*rho0)/n0mean,
              (1+((n1var/(n1mean^2)+1)*n1mean-1)*rho1)/n1mean)
  ExchNum <- c(sum(n0wts*n0vec*(1+(n0vec-1)*rho0)/((1+(n0vec-1)*rhost)^2))/K0,
               sum(n1wts*n1vec*(1+(n1vec-1)*rho1)/((1+(n1vec-1)*rhost)^2))/K1)
  ExchDenom <- c((sum(n0wts*n0vec/(1+(n0vec-1)*rhost))/K0)^2,
                 (sum(n1wts*n1vec/(1+(n1vec-1)*rhost))/K1)^2)
  Psivec <- c(K0,K1)/(K0+K1)
  return(list(rho.star = rhost,
              Cor = 1/(Psivec*CorDenom),
              Ind = IndVal/Psivec,
              Exch = ExchNum/(ExchDenom*Psivec)))
}

Vars <- function(n0vec,n1vec,rho0,rho1,pi0,pi1,
                 n0wts=NULL,n1wts=NULL,K.tot=1) {
  Factors <- VarFactors(n0vec=n0vec,n1vec=n1vec,
                        rho0=rho0,rho1=rho1,
                        n0wts=n0wts,n1wts=n1wts)
  vvec <- c(pi0*(1-pi0),pi1*(1-pi1))
  return(list(rho.star = Factors$rho.star,
              Cor = sum(Factors$Cor/vvec)/K.tot,
              Ind = sum(Factors$Ind/vvec)/K.tot,
              Exch = sum(Factors$Exch/vvec)/K.tot))
}

Plot.Res <- function(resdf,xlab,xvar,xtitle,fileout, truth=NULL,
                     ylimit1=NULL, ylimit2=NULL, horiz=FALSE, 
                     colors=c(1,"#1f78b4","#a6cee3"), 
                     truth2=NULL, 
                     logX=FALSE) {
  setEPS()
  if (horiz) {
    postscript(file=paste0(figfolder,"/",fileout,".eps"),
               width=12, height=6, paper="special")
    par(mfrow=c(1,2))
  } else {
    postscript(file=paste0(figfolder,"/",fileout,".eps"),
               width=6, height=12, paper="special")
    par(mfrow=c(2,1))
  }
  if (is.null(ylimit1)) {
    ylimit1 <- c(floor(min(resdf[,c("indvar","corvar","exchvar")])),ceiling(max(resdf[,c("indvar","corvar","exchvar")])))
  }
  if (!logX) {
    xvals <- resdf[,xvar]
    logval <- ""
    xaxtval=NULL
  } else {
    xvals <- exp(resdf[,xvar])
    if (!is.null(truth)) {
      truth <- exp(truth)
    }
    if (!is.null(truth2)) {
      truth2 <- exp(truth2)
    }
    logval <- "x"
    xaxtval="n"
  }
  plot(x=xvals,y=sqrt(resdf$indvar),type="l",col=colors[3],lty=2,lwd=6,
       xlab=xlab, ylab="Standard Error of Estimate",
       ylim=ylimit1,
       main=bquote("a."~"SE vs."~.(xtitle)~"by WCS"),
       log=logval, xaxt=xaxtval
  )
  if (logX) {
    axis(1, at=c(1/4,1/2,3/4,1,3/2,2,3,4), 
         labels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"))
  }
  lines(x=xvals,y=sqrt(resdf$corvar),type="l",col=colors[1],lty=1,lwd=6)
  lines(x=xvals,y=sqrt(resdf$exchvar),type="l",col=colors[2],lty=3,lwd=6)
  if (!is.null(truth)) {
    abline(v=truth, col=2, lty=4, lwd=3)
  }
  if (!is.null(truth2)) {
    abline(v=truth2, col=2, lty=1, lwd=3)
  }
  legend(x="topleft",legend=c("Correct","Exchangeable","Independent"),
         lty=c(1,3,2),lwd=rep(4,3),col=colors, bg="white")
  
  if (is.null(ylimit2)) {
    ylimit2 <- c(floor(min(c(resdf$corvar/resdf$exchvar*100,resdf$corvar/resdf$indvar*100))),101)
  }
  plot(x=xvals,y=resdf$corvar/resdf$indvar*100,type="l",col=colors[3],lty=2,lwd=6,
       xlab=xlab, ylab="ARE (%)", 
       ylim=ylimit2,
       main=bquote("b."~"ARE vs."~.(xtitle)~"by WCS"),
       log=logval, xaxt=xaxtval
  )
  if (logX) {
    axis(1, at=c(1/4,1/2,3/4,1,3/2,2,3,4), 
         labels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"))
  }
  lines(x=xvals,y=resdf$corvar/resdf$corvar*100,type="l",col=colors[1],lty=1,lwd=6)
  lines(x=xvals,y=resdf$corvar/resdf$exchvar*100,type="l",col=colors[2],lty=3,lwd=6)
  if (!is.null(truth)) {
    abline(v=truth, col=2, lty=4, lwd=3)
  }
  if (!is.null(truth2)) {
    abline(v=truth2, col=2, lty=1, lwd=3)
  }
  legend(x="bottomleft",legend=c("Correct","Exchangeable","Independent"),
         lty=c(1,3,2),lwd=rep(4,3),col=colors, bg="white")
  dev.off()
}

Plot.Res4 <- function(resdf1,resdf2,xlab1,xlab2,xvar1,xvar2,xtitle1,xtitle2,fileout, 
                      truth1=NULL, truth2=NULL, 
                      ylimit1=NULL, ylimit2=NULL, 
                      colors=c(1,"#1f78b4","#a6cee3"),
                      logX1=FALSE, logX2=FALSE,
                      truth1b=NULL, truth2b=NULL) {
  setEPS()
  postscript(file=paste0(figfolder,"/",fileout,".eps"),
               width=12, height=12, paper="special")
  par(mfrow=c(2,2))
  if (is.null(ylimit1)) {
    ylimit1 <- c(floor(min(min(resdf1[,c("indvar","corvar","exchvar")]),min(resdf2[,c("indvar","corvar","exchvar")]))),
                 ceiling(max(max(resdf1[,c("indvar","corvar","exchvar")]),max(resdf2[,c("indvar","corvar","exchvar")]))))
  }
  if (!logX1) {
    xvals1 <- resdf1[,xvar1]
    logval1 <- ""
    xaxtval1 <- NULL
  } else {
    xvals1 <- exp(resdf1[,xvar1])
    if (!is.null(truth1)) {
      truth1 <- exp(truth1)
    }
    if (!is.null(truth1b)) {
      truth1b <- exp(truth1b)
    }
    logval1 <- "x"
    xaxtval1 <- "n"
  }
  plot(x=xvals1,y=sqrt(resdf1$indvar),type="l",col=colors[3],lty=2,lwd=6,
       xlab=xlab1, ylab="Standard Error of Estimate",
       ylim=ylimit1,
       main=bquote("a."~"SE vs."~.(xtitle1)~"by WCS"),
       log=logval1, xaxt=xaxtval1
  )
  if (logX1) {
    axis(1, at=c(1/4,1/2,3/4,1,3/2,2,3,4), 
         labels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"))
  }
  lines(x=xvals1,y=sqrt(resdf1$corvar),type="l",col=colors[1],lty=1,lwd=6)
  lines(x=xvals1,y=sqrt(resdf1$exchvar),type="l",col=colors[2],lty=3,lwd=6)
  if (!is.null(truth1)) {
    abline(v=truth1, col=2, lty=4, lwd=3)
  }
  if (!is.null(truth1b)) {
    abline(v=truth1b, col=2, lty=1, lwd=3)
  }
  legend(x="topright",legend=c("Correct","Exchangeable","Independent"),
         lty=c(1,3,2),lwd=rep(4,3),col=colors,  
         cex=1.5, bg="white")
  
  if (!logX2) {
    xvals2 <- resdf2[,xvar2]
    logval2 <- ""
    xaxtval2 <- NULL
  } else {
    xvals2 <- exp(resdf2[,xvar2])
    if (!is.null(truth2)) {
      truth2 <- exp(truth2)
    }
    if (!is.null(truth2b)) {
      truth2b <- exp(truth2b)
    }
    logval2 <- "x"
    xaxtval2 <- "n"
  }
  plot(x=xvals2,y=sqrt(resdf2$indvar),type="l",col=colors[3],lty=2,lwd=6,
       xlab=xlab2, ylab="Standard Error of Estimate",
       ylim=ylimit1,
       main=bquote("b."~"SE vs."~.(xtitle2)~"by WCS"),
       log=logval2, xaxt=xaxtval2
  )
  if (logX2) {
    axis(1, at=c(1/4,1/2,3/4,1,3/2,2,3,4), 
         labels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"))
  }
  lines(x=xvals2,y=sqrt(resdf2$corvar),type="l",col=colors[1],lty=1,lwd=6)
  lines(x=xvals2,y=sqrt(resdf2$exchvar),type="l",col=colors[2],lty=3,lwd=6)
  if (!is.null(truth2)) {
    abline(v=truth2, col=2, lty=4, lwd=3)
  }
  if (!is.null(truth2b)) {
    abline(v=truth2b, col=2, lty=1, lwd=3)
  }
  legend(x="topright",legend=c("Correct","Exchangeable","Independent"),
         lty=c(1,3,2),lwd=rep(4,3),col=colors,  
         cex=1.5, bg="white")
  
  if (is.null(ylimit2)) {
    ylimit2 <- c(floor(min(c(resdf1$corvar/resdf1$exchvar*100,resdf1$corvar/resdf1$indvar*100,
                             resdf2$corvar/resdf2$exchvar*100,resdf2$corvar/resdf2$indvar*100))),101)
  }
  plot(x=xvals1,y=resdf1$corvar/resdf1$indvar*100,type="l",col=colors[3],lty=2,lwd=6,
       xlab=xlab1, ylab="ARE (%)", 
       ylim=ylimit2,
       main=bquote("c."~"ARE vs."~.(xtitle1)~"by WCS"),
       log=logval1, xaxt=xaxtval1
  )
  if (logX1) {
    axis(1, at=c(1/4,1/2,3/4,1,3/2,2,3,4), 
         labels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"))
  }
  lines(x=xvals1,y=resdf1$corvar/resdf1$corvar*100,type="l",col=colors[1],lty=1,lwd=6)
  lines(x=xvals1,y=resdf1$corvar/resdf1$exchvar*100,type="l",col=colors[2],lty=3,lwd=6)
  if (!is.null(truth1)) {
    abline(v=truth1, col=2, lty=4, lwd=3)
  }
  if (!is.null(truth1b)) {
    abline(v=truth1b, col=2, lty=1, lwd=3)
  }
  legend(x="bottomleft",legend=c("Correct","Exchangeable","Independent"),
         lty=c(1,3,2),lwd=rep(4,3),col=colors,  
         cex=1.5, bg="white")
  
  plot(x=xvals2,y=resdf2$corvar/resdf2$indvar*100,type="l",col=colors[3],lty=2,lwd=6,
       xlab=xlab2, ylab="ARE (%)", 
       ylim=ylimit2,
       main=bquote("d."~"ARE vs."~.(xtitle2)~"by WCS"),
       log=logval2, xaxt=xaxtval2
  )
  if (logX2) {
    axis(1, at=c(1/4,1/2,3/4,1,3/2,2,3,4), 
         labels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"))
  }
  lines(x=xvals2,y=resdf2$corvar/resdf2$corvar*100,type="l",col=colors[1],lty=1,lwd=6)
  lines(x=xvals2,y=resdf2$corvar/resdf2$exchvar*100,type="l",col=colors[2],lty=3,lwd=6)
  if (!is.null(truth2)) {
    abline(v=truth2, col=2, lty=4, lwd=3)
  }
  if (!is.null(truth2b)) {
    abline(v=truth2b, col=2, lty=1, lwd=3)
  }
  legend(x="bottomleft",legend=c("Correct","Exchangeable","Independent"),
         lty=c(1,3,2),lwd=rep(4,3),col=colors,  
         cex=1.5, bg="white")
  dev.off()
}

Plot.SS <- function(resdf,xlab,xvar,title=NULL,truth=NULL,
                    ylimit=NULL, colors=c(1,"#33a02c","#b2df8a"),
                    logX=FALSE, truthb=NULL) {
  if (is.null(ylimit)) {
    ylimit <- c(floor(min(c(resdf$exchvar.eq/resdf$exchvar*100,resdf$indvar.eq/resdf$indvar*100,99))),
                ceiling(max(c(resdf$exchvar.eq/resdf$exchvar*100,resdf$indvar.eq/resdf$indvar*100,101))))
  }
  if (!logX) {
    xvals <- resdf[,xvar]
    logval <- ""
    xaxtval <- NULL
  } else {
    xvals <- exp(resdf[,xvar])
    logval <- "x"
    if (!is.null(truth)) {
      truth <- exp(truth)
    }
    if (!is.null(truthb)) {
      truthb <- exp(truthb)
    }
    xaxtval <- "n"
  }
  plot(x=xvals,y=resdf$exchvar.eq/resdf$exchvar*100,type="l",col=colors[2],lty=3,lwd=6,
       xlab=xlab, ylab="Estimated SS/Required SS (%)",
       ylim=ylimit,
       main=title, log=logval, xaxt=xaxtval)
  if (logX) {
    axis(1, at=c(1/4,1/2,3/4,1,3/2,2,3,4), 
         labels=c("0.25","0.5","0.75","1.0","1.5","2.0","3.0","4.0"))
  }
  lines(x=xvals,y=resdf$indvar.eq/resdf$indvar*100,type="l",col=colors[3],lty=2,lwd=6)
  lines(x=xvals,y=resdf$corvar/resdf$corvar*100,type="l",col=colors[1],lty=1,lwd=6)
  if (!is.null(truth)) {
    abline(v=truth, col=2, lty=4, lwd=3)
  }
  if (!is.null(truthb)) {
    abline(v=truthb, col=2, lty=1, lwd=3)
  }
}

Plot.SS.Full <- function(resdf,xlab,xvar,title=NULL,truth=NULL,
                         ylimit=NULL, colors=c(1,"#33a02c","#b2df8a","#756bb1","#bcbddc"),
                         logX=FALSE, truthb=NULL) {
  Plot.SS(resdf,xlab,xvar,title,truth,ylimit,colors[1:3],logX,truthb)
  if (is.null(ylimit)) {
    ylimit <- c(floor(min(c(resdf$exchvar.eq/resdf$exchvar*100,resdf$indvar.eq/resdf$indvar*100,99))),
                ceiling(max(c(resdf$exchvar.eq/resdf$exchvar*100,resdf$indvar.eq/resdf$indvar*100,101))))
  }
  if (!logX) {
    xvals <- resdf[,xvar]
    logval <- ""
    xaxtval <- NULL
  } else {
    xvals <- exp(resdf[,xvar])
    logval <- "x"
    if (!is.null(truth)) {
      truth <- exp(truth)
    }
    if (!is.null(truthb)) {
      truthb <- exp(truthb)
    }
    xaxtval <- "n"
  }
  lines(x=xvals,y=resdf$exchvar.eq.rho0/resdf$exchvar*100,type="l",col=colors[4],lty=5,lwd=6)
  lines(x=xvals,y=resdf$exchvar.eq.rho1/resdf$exchvar*100,type="l",col=colors[5],lty=6,lwd=6)
  if (!is.null(truth)) {
    abline(v=truth, col=2, lty=4, lwd=3)
  }
  if (!is.null(truthb)) {
    abline(v=truthb, col=2, lty=1, lwd=3)
  }
}

##### Two Binary Cluster-Level Covariates #####

Vars2 <- function(n00vec, n10vec, n01vec, n11vec, 
                  pi00, pi10, pi01, pi11,
                  rho00, rho10, rho01, rho11,
                  n00wts=NULL, n10wts=NULL, n01wts=NULL, n11wts=NULL,
                  K.tot=NULL) {
  ## Note: parameter of interest corresponds to first 0/1 index
  if (is.null(n00wts) + is.null(n10wts) + is.null(n01wts) + is.null(n11wts) < 4 & is.null(n00wts) + is.null(n01wts) + is.null(n10wts) + is.null(n11wts) > 0) {
    print("Error: Specify Weights for All Cluster Size Vectors or None")
  } else if (is.null(n00wts)) {
    n00wts <- rep(1,length(n00vec))
    n10wts <- rep(1,length(n10vec))
    n01wts <- rep(1,length(n01vec))
    n11wts <- rep(1,length(n11vec))
  }
  
  K00 <- sum(n00wts)
  K10 <- sum(n10wts)
  K01 <- sum(n01wts)
  K11 <- sum(n11wts)
  Kvec <- c(K00,K10,K01,K11)
  Psivec <- Kvec/sum(Kvec)
  
  meanvec <- c(sum(n00vec*n00wts)/K00,sum(n10vec*n10wts)/K10,
               sum(n01vec*n01wts)/K01,sum(n11vec*n11wts)/K11)
  varvec <- c(sum(n00wts*(n00vec-meanvec[1])^2)/K00,sum(n10wts*(n10vec-meanvec[2])^2)/K10,
              sum(n01wts*(n01vec-meanvec[3])^2)/K01,sum(n11wts*(n11vec-meanvec[4])^2)/K11)
  
  vvec <- unname(unlist(c(pi00*(1-pi00),pi10*(1-pi10),pi01*(1-pi01),pi11*(1-pi11))))
  
  qvec <- c(sum(n00wts*n00vec/(1+(n00vec-1)*rho00)),sum(n10wts*n10vec/(1+(n10vec-1)*rho10)),
            sum(n01wts*n01vec/(1+(n10vec-1)*rho01)),
            sum(n11wts*n11vec/(1+(n11vec-1)*rho11)))/Kvec
  Cor <- ((sum(1/(vvec[1:2]*Psivec[1:2]*qvec[1:2])))^(-1) + (sum(1/(vvec[3:4]*Psivec[3:4]*qvec[3:4])))^(-1))^(-1)
  
  bvecInd <- vvec*Psivec*meanvec
  cvecInd <- bvecInd*(1+(((varvec/(meanvec^2))+1)*meanvec-1)*unname(unlist(c(rho00,rho10,rho01,rho11))))

  pairs00 <- n00vec*(n00vec-1)/2
  pairs10 <- n10vec*(n10vec-1)/2
  pairs01 <- n01vec*(n01vec-1)/2
  pairs11 <- n11vec*(n11vec-1)/2
  rhost <- unname(unlist((sum(n00wts*pairs00)*rho00+sum(n10wts*pairs10)*rho10+sum(n01wts*pairs01)*rho01+sum(n11wts*pairs11)*rho11)/(sum(n00wts*pairs00)+sum(n10wts*pairs10)+sum(n01wts*pairs01)+sum(n11wts*pairs11))))
  
  bvecExch <- vvec*Psivec*c(sum(n00wts*n00vec/(1+(n00vec-1)*rhost)),
                            sum(n10wts*n10vec/(1+(n10vec-1)*rhost)),
                            sum(n01wts*n01vec/(1+(n01vec-1)*rhost)),
                            sum(n11wts*n11vec/(1+(n11vec-1)*rhost)))/Kvec
  cvecExch <- vvec*Psivec*c(sum(n00wts*n00vec*(1+(n00vec-1)*rho00)/(1+(n00vec-1)*rhost)^2),
                            sum(n10wts*n10vec*(1+(n10vec-1)*rho10)/(1+(n10vec-1)*rhost)^2),
                            sum(n01wts*n01vec*(1+(n01vec-1)*rho01)/(1+(n01vec-1)*rhost)^2),
                            sum(n11wts*n11vec*(1+(n11vec-1)*rho11)/(1+(n11vec-1)*rhost)^2))/Kvec
  
  bcFunc <- function(b,c) {
    Num1 <- (c[1]*b[2]^2+c[2]*b[1]^2)*(b[3]+b[4])^2
    Num2 <- (c[3]*b[4]^2+c[4]*b[3]^2)*(b[1]+b[2])^2
    Denom <- (b[1]+b[2])*b[3]*b[4]+(b[3]+b[4])*b[1]*b[2]
    return((Num1+Num2)/(Denom^2))
  }
  
  Ind <- bcFunc(bvecInd,cvecInd)
  Exch <- bcFunc(bvecExch,cvecExch)
  
  return(list(rho.star = rhost,
              Cor = Cor/K.tot,
              Ind = Ind/K.tot,
              Exch = Exch/K.tot))
}