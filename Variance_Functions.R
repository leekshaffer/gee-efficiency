#### Useful Functions ####
expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

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

SSest <- function(n0vec,n1vec,rho0,rho1,pi0,pi1,
                  n0wts=NULL,n1wts=NULL,
                  sig=0.05,pwr=.8) {
  if (is.null(n0wts) + is.null(n1wts) == 1) {
    print("Error: Specify Weights for both n0 and n1 or Neither")
  } else if (is.null(n0wts) & is.null(n1wts)) {
    n0wts <- rep(1,length(n0vec))
    n1wts <- rep(1,length(n1vec))
  }
  AltVarRes <- Vars(n0vec,n1vec,rho0,rho1,pi0,pi1,
                    n0wts,n1wts,K.tot=1)
  AltVars <- unlist(AltVarRes[c("Cor","Ind","Exch")])
  rhost <- AltVarRes$rho.star
  NullVars <- unlist(Vars(n0vec,n1vec,rhost,rhost,
                          pi0,pi0,n0wts,n1wts,K.tot=1)[c("Cor","Ind","Exch")])
  AltHyp <- log((pi1/(1-pi1))/(pi0/(1-pi0)))
  K.req <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars)+qnorm(pwr)*sqrt(AltVars))^2
  
  AltVars.rho0 <- unlist(Vars(n0vec,n1vec,rho0,rho0,pi0,pi1,
                              n0wts,n1wts,K.tot=1)[c("Cor","Ind","Exch")])
  NullVars.rho0 <- unlist(Vars(n0vec,n1vec,rho0,rho0,
                               pi0,pi0,n0wts,n1wts,K.tot=1)[c("Cor","Ind","Exch")])
  K.est.rho0 <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars.rho0)+qnorm(pwr)*sqrt(AltVars.rho0))^2
  
  AltVars.rho1 <- unlist(Vars(n0vec,n1vec,rho1,rho1,pi0,pi1,
                              n0wts,n1wts,K.tot=1)[c("Cor","Ind","Exch")])
  NullVars.rho1 <- unlist(Vars(n0vec,n1vec,rho1,rho1,
                               pi0,pi0,n0wts,n1wts,K.tot=1)[c("Cor","Ind","Exch")])
  K.est.rho1 <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars.rho1)+qnorm(pwr)*sqrt(AltVars.rho1))^2
  
  AltVars.rhost <- unlist(Vars(n0vec,n1vec,rhost,rhost,pi0,pi1,
                              n0wts,n1wts,K.tot=1)[c("Cor","Ind","Exch")])
  NullVars.rhost <- unlist(Vars(n0vec,n1vec,rhost,rhost,
                               pi0,pi0,n0wts,n1wts,K.tot=1)[c("Cor","Ind","Exch")])
  K.est.rhost <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars.rhost)+qnorm(pwr)*sqrt(AltVars.rhost))^2
  
  return(list(SSratio.cor=K.req["Cor"]/K.req["Cor"],
              SSratio.ind=K.est.rhost["Ind"]/K.req["Ind"],
              SSratio.exch.rhost=K.est.rhost["Exch"]/K.req["Exch"],
              SSratio.exch.rho0=K.est.rho0["Exch"]/K.req["Exch"],
              SSratio.exch.rho1=K.est.rho1["Exch"]/K.req["Exch"]))
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
  
  zetavec <- c(sum(n00wts*n00vec/(1+(n00vec-1)*rho00)),
            sum(n10wts*n10vec/(1+(n10vec-1)*rho10)),
            sum(n01wts*n01vec/(1+(n01vec-1)*rho01)),
            sum(n11wts*n11vec/(1+(n11vec-1)*rho11)))/Kvec
  Cor <- ((sum(1/(vvec[1:2]*Psivec[1:2]*zetavec[1:2])))^(-1) + (sum(1/(vvec[3:4]*Psivec[3:4]*zetavec[3:4])))^(-1))^(-1)
  
  bvecInd <- vvec*Psivec*meanvec
  cvecInd <- bvecInd*(1+(((varvec/(meanvec^2))+1)*meanvec-1)*unname(unlist(c(rho00,rho10,rho01,rho11))))

  pairs00 <- n00vec*(n00vec-1)/2
  pairs10 <- n10vec*(n10vec-1)/2
  pairs01 <- n01vec*(n01vec-1)/2
  pairs11 <- n11vec*(n11vec-1)/2
  rhost <- unname(unlist((sum(n00wts*pairs00)*rho00+sum(n10wts*pairs10)*rho10+sum(n01wts*pairs01)*rho01+sum(n11wts*pairs11)*rho11)/(sum(n00wts*pairs00)+sum(n10wts*pairs10)+sum(n01wts*pairs01)+sum(n11wts*pairs11))))
  
  tauvecExch <- vvec*Psivec*c(sum(n00wts*n00vec/(1+(n00vec-1)*rhost)),
                            sum(n10wts*n10vec/(1+(n10vec-1)*rhost)),
                            sum(n01wts*n01vec/(1+(n01vec-1)*rhost)),
                            sum(n11wts*n11vec/(1+(n11vec-1)*rhost)))/Kvec
  omegavecExch <- vvec*Psivec*c(sum(n00wts*n00vec*(1+(n00vec-1)*rho00)/(1+(n00vec-1)*rhost)^2),
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
  Exch <- bcFunc(tauvecExch,omegavecExch)
  
  return(list(rho.star = rhost,
              Cor = Cor/K.tot,
              Ind = Ind/K.tot,
              Exch = Exch/K.tot))
}

SSest2 <- function(n00vec, n10vec, n01vec, n11vec, 
                   pi00, pi10, pi01, pi11,
                   rho00, rho10, rho01, rho11,
                   n00wts=NULL, n10wts=NULL, n01wts=NULL, n11wts=NULL,
                   sig=0.05,pwr=0.80) {
  if (is.null(n00wts) + is.null(n10wts) + is.null(n01wts) + is.null(n11wts) < 4 & is.null(n00wts) + is.null(n01wts) + is.null(n10wts) + is.null(n11wts) > 0) {
    print("Error: Specify Weights for All Cluster Size Vectors or None")
  } else if (is.null(n00wts)) {
    n00wts <- rep(1,length(n00vec))
    n10wts <- rep(1,length(n10vec))
    n01wts <- rep(1,length(n01vec))
    n11wts <- rep(1,length(n11vec))
  }
  AltVarRes <- Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi10, pi01, pi11,
                     rho00, rho10, rho01, rho11, n00wts, n10wts, n01wts, n11wts,K.tot=1)
  AltVars <- unlist(AltVarRes[c("Cor","Ind","Exch")])
  rhost <- AltVarRes$rho.star
  rhost.0 <- (rho00*sum(n00vec*(n00vec-1)*n00wts)+rho10*sum(n10vec*(n10vec-1)*n10wts))/(sum(n00vec*(n00vec-1)*n00wts)+sum(n10vec*(n10vec-1)*n10wts))
  rhost.1 <- (rho01*sum(n01vec*(n01vec-1)*n01wts)+rho11*sum(n11vec*(n11vec-1)*n11wts))/(sum(n01vec*(n01vec-1)*n01wts)+sum(n11vec*(n11vec-1)*n11wts))
  rho0.st <- (rho00*sum(n00vec*(n00vec-1)*n00wts)+rho01*sum(n01vec*(n01vec-1)*n01wts))/(sum(n00vec*(n00vec-1)*n00wts)+sum(n01vec*(n01vec-1)*n01wts))
  rho1.st <- (rho10*sum(n10vec*(n10vec-1)*n10wts)+rho11*sum(n11vec*(n11vec-1)*n11wts))/(sum(n10vec*(n10vec-1)*n10wts)+sum(n11vec*(n11vec-1)*n11wts))
  NullVars <- unlist(Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi00, pi01, pi01,
                           rhost.0, rhost.0, rhost.1, rhost.1, n00wts, n10wts, n01wts, n11wts,
                           K.tot=1)[c("Cor","Ind","Exch")])
  beta.0 <- log((pi10/(1-pi10))/(pi00/(1-pi00)))
  beta.1 <- log((pi11/(1-pi11))/(pi01/(1-pi01)))
  AltHyp <- (beta.0*sum(n00wts,n10wts)+beta.1*sum(n01wts,n11wts))/sum(n00wts,n10wts,n01wts,n11wts)
  K.req <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars)+qnorm(pwr)*sqrt(AltVars))^2
  
  AltVars.rho0 <- unlist(Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi10, pi01, pi11,
                               rho0.st, rho0.st, rho0.st, rho0.st, n00wts, n10wts, n01wts, n11wts,
                               K.tot=1)[c("Cor","Ind","Exch")])
  NullVars.rho0 <- unlist(Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi00, pi01, pi01,
                                rho0.st, rho0.st, rho0.st, rho0.st, n00wts, n10wts, n01wts, n11wts,
                                K.tot=1)[c("Cor","Ind","Exch")])
  K.est.rho0 <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars.rho0)+qnorm(pwr)*sqrt(AltVars.rho0))^2
  
  AltVars.rho1 <- unlist(Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi10, pi01, pi11,
                               rho1.st, rho1.st, rho1.st, rho1.st, n00wts, n10wts, n01wts, n11wts,
                               K.tot=1)[c("Cor","Ind","Exch")])
  NullVars.rho1 <- unlist(Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi00, pi01, pi01,
                                rho1.st, rho1.st, rho1.st, rho1.st, n00wts, n10wts, n01wts, n11wts,
                                K.tot=1)[c("Cor","Ind","Exch")])
  K.est.rho1 <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars.rho1)+qnorm(pwr)*sqrt(AltVars.rho1))^2
  
  AltVars.rhost <- unlist(Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi10, pi01, pi11,
                                rhost, rhost, rhost, rhost, n00wts, n10wts, n01wts, n11wts,
                                K.tot=1)[c("Cor","Ind","Exch")])
  NullVars.rhost <- unlist(Vars2(n00vec, n10vec, n01vec, n11vec, pi00, pi00, pi01, pi01,
                                 rhost, rhost, rhost, rhost, n00wts, n10wts, n01wts, n11wts,
                                 K.tot=1)[c("Cor","Ind","Exch")])
  K.est.rhost <- 1/(AltHyp^2)*(qnorm(1-sig/2)*sqrt(NullVars.rhost)+qnorm(pwr)*sqrt(AltVars.rhost))^2
  
  return(list(SSratio.cor=K.req["Cor"]/K.req["Cor"],
              SSratio.ind=K.est.rhost["Ind"]/K.req["Ind"],
              SSratio.exch.rhost=K.est.rhost["Exch"]/K.req["Exch"],
              SSratio.exch.rho0=K.est.rho0["Exch"]/K.req["Exch"],
              SSratio.exch.rho1=K.est.rho1["Exch"]/K.req["Exch"]))
}
