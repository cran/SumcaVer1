#' Post model selection MSPE estimation in FH model with mean mis-specification using Datta-Hall-Mandal method. Calculate the post-model selection mspe of Fay-Herriot model with mean mis-specification using Datta-Hall-Mandal method.
#'
#' @param m number of small areas
#' @param p number of fixed model parameters
#' @param X covariates
#' @param beta1 regression coefficients
#' @param beta2 regression coefficients
#' @param A variance of area-specific random effects
#' @param D sampling variance
#' @param R number of simulation runs
#'
#' @return Par: return estimation of model parameters
#' @return MSPE.TRUE.Final: return empirical MSPE of small area predictor
#' @return mspe.DHM.Final: return mspe of small area predictor using the Datta-Hall-Mandal method
#' @return RB.DHM: return relative bias (RB) of mspe of small area predictor using the Datta-Hall-Mandal method
#' @return Rate: return the probability of rejection (nominal level= 0.2)
#' @export
#'
#' @importFrom stats rnorm
#' @importFrom stats qchisq
#' @importFrom psych tr
#' @importFrom stats runif
#'
#' @examples mspe_PMS_Mis_FH_DHM(20,3,matrix(runif(60,0,1),nrow=20,byrow=TRUE),
#' c(1,1,1),c(1,1,1),10,2.5,10)
#'
mspe_PMS_Mis_FH_DHM=function(m,p,X,beta1,beta2,A,D,R)
{

#=============================================
#Define the constants of the data and model

n=m/2

Di.a=runif(n,0.5,1.5)
Di.b=runif(n,0.5,1.5)
D<-c(Di.a, Di.b)
SD<-sqrt(D)

BETA.HAT.PR<-matrix(0,R,p)
A.HAT.PR<-rep(0,R)
MSPE<-matrix(0,R,m)
mspe.dhm<-matrix(0,R,m)

T.DHM<-rep(0,R)

Est.A<-0
Est.BETA<-rep(0,p)
#-----------------------------------
for(r in 1:R){

  random1<-rnorm(m,0,1)
  random2<-rnorm(m,0,1)
  yib<-rep(0,m)
  mu<-rep(0,m)
  for(i in 1:n){
    yib[i]<-as.numeric(X[i,]%*%beta1)+SD[i]*random2[i] +sqrt(A)*random1[i]
    mu[i]<-as.numeric(X[i,]%*%beta1) +sqrt(A)*random1[i]
  }
  for(i in (n+1):m){
    yib[i]<-as.numeric(X[i,]%*%beta2)+SD[i]*random2[i]+sqrt(A)*random1[i]
    mu[i]<-as.numeric(X[i,]%*%beta2) +sqrt(A)*random1[i]
  }
  ###################################
  PX<-X%*% solve(t(X)%*%X)%*%t(X)
  PXI<-diag(1,m)-PX
  D1<-diag(SD^2)
  #################Datta approach######################
  T.DHM[r]<-t(yib)%*%(solve(D1)-solve(D1)%*%X%*%(solve(t(X)%*%solve(D1)%*%X)) %*%t(X)%*%solve(D1))%*%yib
  #####We use the following estimators (based on the PR and EBLUE approach):

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
  ###############################EBLUP
  mu.hat<-rep(0,m)

  if(T.DHM[r]>qchisq(.80, df=m-p)){
    for(i in 1:m){
      mu.hat[i]<-(A.hat*yib[i])/(A.hat+D[i]) +(D[i]*X[i,]%*%BETA.hat)/(A.hat+D[i])}}

  if(T.DHM[r]<= qchisq(.80, df=m-p)){
    for(i in 1:m){
      mu.hat[i]<-t(X[i,])%*%(solve(t(X)%*%solve(D1)%*%X))%*%t(X)%*%solve(D1)%*%yib}}
  ############################################MSPE:
  for(i in 1:m){
    MSPE[r,i]<-(mu.hat[i]-mu[i])^2
  }
  ###############################################
  #MSPE estimation using the DHM method:

  if(T.DHM[r]>qchisq(.80, df=m-p)){
    for(i in 1:m){
      mspe.dhm[r,i]<-(A.hat*D[i])/(A.hat+D[i])+((D[i]/(A.hat+D[i]))^2)*(t(X[i,])%*%solve(sumx.tem)%*%(X[i,]))+((4*D[i]^2)/((A.hat+D[i])^3*m^2))*tem1
    }
  }

  if(T.DHM[r]<= qchisq(.80, df=m-p)){
    for(i in 1:m){
      mspe.dhm[r,i]<-t(X[i,])%*%(solve(t(X)%*%solve(D1)%*%X))%*%(X[i,])
    }
  }
  #---------------------------------------------
    }#loop R (simulation runs)
######################################################
BETA.HAT<-apply(BETA.HAT.PR,2,mean)
A.HAT<-mean(A.HAT.PR)
MSPE.TRUE<-apply(MSPE,2,mean)
mspe.DHM<-apply(mspe.dhm,2,mean)

Par<-c(BETA.HAT, A.HAT)
MSPE.TRUE.Final<-c(MSPE.TRUE)
mspe.DHM.Final<-c(mspe.DHM)

RB.DHM<-(mspe.DHM.Final/MSPE.TRUE.Final)-1

#rate of alpha:
numm<-0
for(r in 1:R){if(T.DHM[r]>qchisq(.80, df=m-p)){numm<-numm+1}}
Rate<-numm/R

list (Par=Par, MSPE.TRUE.Final=MSPE.TRUE.Final, mspe.DHM.Final=mspe.DHM.Final, RB.DHM=RB.DHM, Rate=Rate) #Result=Result)

}#loop function

