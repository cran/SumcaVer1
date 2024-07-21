#' MSPE estimation in FH model using Prasad-Rao method. Calculate the mspe of Fay-Herriot model in SAE using Prasad-Rao method.
#'
#' @param m number of small areas
#' @param p number of fixed model parameters
#' @param X Covariates
#' @param beta regression coefficients
#' @param A variance of area-specific random effects
#' @param D sampling variance
#' @param R number of simulation runs
#'
#' @return Par: return estimation of model parameters
#' @return MSPE.TRUE.Final: return empirical MSPE of small area predictor
#' @return mspe.PR.Final: return mspe of small area predictor using the Prasad-Rao method
#' @return RB.PR: return relative bias (RB) of mspe of small area predictor using the Prasad-Rao method
#' @export
#'
#' @importFrom stats  rnorm
#' @importFrom psych tr
#' @importFrom stats runif
#'
#' @examples mspe_FH_PR(20,3,matrix(runif(60,0,1),nrow=20,byrow=TRUE),c(1,1,1),10,2.5,10)
#'
#'
mspe_FH_PR=function(m,p,X,beta,A,D,R)
{
#=============================================
#Define the constants of the data and model

n=m/2

Di.a=runif(n,3.5,4.5)
Di.b=runif(n,0.5,1.5)
D<-c(Di.a, Di.b)
SD<-sqrt(D)


BETA.HAT.PR<-matrix(0,R,p)
A.HAT.PR<-rep(0,R)
MSPE<-matrix(0,R,m)
mspe.pr<-matrix(0,R,m)

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
  A.HAT.PR[r]<-A.hat
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
  BETA.HAT.PR[r,]<-BETA.hat


  Est.A<-Est.A+A.hat
  Est.BETA<-Est.BETA+BETA.hat

  ############################################MSPE:
  for(i in 1:m){
    MSPE[r,i]<-((A.hat*yib[i])/(A.hat+D[i]) +(D[i]*X[i,]%*%BETA.hat)/(A.hat+D[i])-mu[i])^2
  }
  ###############################################
  #MSPE estimation using the PR method:
  for(i in 1:m){
    mspe.pr[r,i]<-(A.hat*D[i])/(A.hat+D[i])+((D[i]/(A.hat+D[i]))^2)*(t(X[i,])%*%solve(sumx.tem)%*%(X[i,]))+((4*D[i]^2)/((A.hat+D[i])^3*m^2))*tem1
   }
}#loop R (simulation runs)
######################################################
BETA.HAT<-apply(BETA.HAT.PR,2,mean)
A.HAT<-mean(A.HAT.PR)
MSPE.TRUE<-apply(MSPE,2,mean)
mspe.PR<-apply(mspe.pr,2,mean)

Par<-c(BETA.HAT, A.HAT)
MSPE.TRUE.Final<-c(MSPE.TRUE)
mspe.PR.Final<-c(mspe.PR)

RB.PR<-(mspe.PR.Final/MSPE.TRUE.Final)-1

list (Par=Par, MSPE.TRUE.Final=MSPE.TRUE.Final, mspe.PR.Final=mspe.PR.Final, RB.PR=RB.PR) #Result=Result)

}#loop function
