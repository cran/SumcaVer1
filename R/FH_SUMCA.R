#' MSPE estimation in FH model using SUMCA method. Calculate the mspe of Fay-Herriot model in SAE using Sumca method.
#'
#' @param m number of small areas
#' @param p number of fixed model parameters
#' @param X covariates
#' @param beta regression coefficients
#' @param A variance of area-specific random effects
#' @param D sampling variance
#' @param K number of Monte Carlo for the SUMCA method
#' @param R number of simulation runs
#'
#' @return Par: return estimation of model parameters
#' @return MSPE.TRUE.Final: return empirical MSPE of small area predictor
#' @return mspe.Sumca.Final: return mspe of small area predictor using the SUMCA method
#' @return RB.SUMCA: return relative bias (RB) of mspe of small area predictor using the SUMCA method
#' @export
#'
#' @importFrom stats rnorm
#' @importFrom psych tr
#' @importFrom stats runif
#'
#' @examples mspe_FH_Sumca(20,3,matrix(runif(60,0,1),nrow=20,byrow=TRUE),c(1,1,1),10,2.5,10,10)
#'
#'
mspe_FH_Sumca=function(m,p,X,beta,A,D,K,R)
{
#=============================================
#Define the constants of the data and model

n=m/2



Di.a=runif(n,3.5,4.5)
Di.b=runif(n,0.5,1.5)
D<-c(Di.a, Di.b)
SD<-sqrt(D)


BETA.HAT<-matrix(0,R,p)
A.HAT<-rep(0,R)
MSPE<-matrix(0,R,m)
mspe.new<-matrix(0,R,m)

Est.A<-0
Est.BETA<-rep(0,p)
#-----------------------------------
for(r in 1:R){

  random1<-rnorm(m,0,1)
  random2<-rnorm(m,0,1)
  yib<-rep(0,m)
  mu<-rep(0,m)
  for(i in 1:n){
    yib[i]<-as.numeric(X[i,]%*%beta)+sqrt(A)*random1[i]+SD[i]*random2[i]
    mu[i]<-as.numeric(X[i,]%*%beta)+sqrt(A)*random1[i]
  }
  for(i in (n+1):m){
    yib[i]<-as.numeric(X[i,]%*%beta)+sqrt(A)*random1[i]+SD[i]*random2[i]
    mu[i]<-as.numeric(X[i,]%*%beta)+sqrt(A)*random1[i]
  }
#----------------------------------------------------model
  PX<-X%*% solve(t(X)%*%X)%*%t(X)
  PXI<-diag(1,m)-PX
  D1<-diag(SD^2)
  #####We use the following estimators (based on the PR and EBLUE) approach:
  A.hat<-(t(yib)%*%PXI%*%yib-tr(PXI%*%D1))/(m-p)
  A.HAT[r]<-A.hat
  A.hat<-ifelse(A.hat<=0,0.0001,A.hat)

  sumx.tem<-matrix(0,p,p)
  sumy.tem<-matrix(0,p,1)
  tem1=0
  for(i in 1:m){

    sumx.tem<-sumx.tem+( X[i,]%*%t(X[i,] ) ) /drop(A.hat+SD[i]^2)
    sumy.tem<-sumy.tem+( X[i,]*yib[i] )/drop(A.hat+SD[i]^2)
    tem1<-tem1+(A.hat+D[i])^2
  }
  BETA.hat<-solve(sumx.tem)%*%sumy.tem
  BETA.HAT[r,]<-BETA.hat


  Est.A<-Est.A+A.hat
  Est.BETA<-Est.BETA+BETA.hat

  ############################################MSPE:
  for(i in 1:m){
    MSPE[r,i]<-((A.hat*yib[i])/(A.hat+D[i]) +(D[i]*X[i,]%*%BETA.hat)/(A.hat+D[i])-mu[i])^2
  }
  ###############################################
  #MSPE estimation using the SUMCA method:

  bias.MC<-matrix(0,K,m)

  for(k in 1:K){
    random1<-rnorm(m,0,1)
    random2<-rnorm(m,0,1)
    yib.k<-rep(0,m)
    for(i in 1:n){
      yib.k[i]<-as.numeric(X[i,]%*%BETA.hat)+sqrt(A.hat)*random1[i]+SD[i]*random2[i]
    }
    for(i in (n+1):m){
      yib.k[i]<-as.numeric(X[i,]%*%BETA.hat)+sqrt(A.hat)*random1[i]+SD[i]*random2[i]
    }

    A.hat.k<-(t(yib.k)%*%PXI%*%yib.k-tr(PXI%*%D1))/(m-p)
    A.hat.k<-ifelse(A.hat.k<=0,0.0001,A.hat.k)


    sumx.tem.k<-matrix(0,p,p)
    sumy.tem.k<-matrix(0,p,1)
    for(i in 1:m){
      sumx.tem.k<-sumx.tem.k+( X[i,]%*%t(X[i,] ) ) /drop(A.hat.k+SD[i]^2)
      sumy.tem.k<-sumy.tem.k+( X[i,]*yib.k[i] )/drop(A.hat.k+SD[i]^2)
    }
    BETA.hat.k<-solve(sumx.tem.k)%*%sumy.tem.k
    ############################################
    for(i in 1:m){
      bias.MC[k,i]<-bias.MC[k,i]+D[i]*(A.hat/(A.hat+D[i])-A.hat.k/(A.hat.k+D[i]))+( (A.hat.k/(A.hat.k+D[i])-A.hat/(A.hat+D[i]))*yib.k[i]+D[i]*(t(X[i,])%*%BETA.hat.k/(A.hat.k+D[i])-t(X[i,])%*%BETA.hat/(A.hat+D[i])))^2
    }#loop m
  }#loop K

  bias.co<-apply(bias.MC,2,mean)

  for(i in 1:m){
    mspe.new[r,i]<-(A.hat*D[i])/(A.hat+D[i])+bias.co[i]
  }

}#loop R (simulation runs)
######################################################
BETA.HAT<-apply(BETA.HAT,2,mean)
A.HAT<-mean(A.HAT)
MSPE.TRUE<-apply(MSPE,2,mean)
mspe.Sumca<-apply(mspe.new,2,mean)

Par<-c(BETA.HAT, A.HAT)
MSPE.TRUE.Final<-c(MSPE.TRUE)
mspe.Sumca.Final<-c(mspe.Sumca)

RB.SUMCA<-(mspe.Sumca.Final/MSPE.TRUE.Final)-1

list (Par=Par, MSPE.TRUE.Final=MSPE.TRUE.Final, mspe.Sumca.Final=mspe.Sumca.Final, RB.SUMCA=RB.SUMCA) #Result=Result)

}#loop function
