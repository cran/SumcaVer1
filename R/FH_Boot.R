#' MSPE estimation in FH model using double-phase bootstrap method.Calculate the mspe of Fay-Herriot model in SAE using double-phase bootstrap method.
#'
#' @param m number of small areas
#' @param p number of fixed model parameters
#' @param X covariates
#' @param beta regression coefficients
#' @param A variance of area-specific random effects
#' @param D sampling variance
#' @param B1 number of first-phase bootstrap method
#' @param B2 number of second-phase bootstrap method
#' @param R number of simulation runs
#'
#' @return Par: return estimation of model parameters
#' @return MSPE.TRUE.Final: return empirical MSPE of small area predictor
#' @return mspe.Boot1.Final: return mspe of small area predictor using the bootstrap method 1
#' @return mspe.Boot2.Final: return mspe of small area predictor using the bootstrap method 2
#' @return RB.Boot1: return relative bias (RB) of mspe of small area predictor using the bootstrap method 1
#' @return RB.Boot2: return relative bias (RB) of mspe of small area predictor using the bootstrap method 2
#' @export
#'
#' @importFrom stats rnorm
#' @importFrom psych tr
#' @importFrom stats runif
#'
#' @examples mspe_FH_Boot(20,3,matrix(runif(60,0,1),nrow=20,byrow=TRUE),c(1,1,1),10,2.5,20,20,10)
#'
#'
#'
mspe_FH_Boot=function(m,p,X,beta,A,D,B1,B2,R)
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
  mspe.boot1<-matrix(0,R,m)
  mspe.boot2<-matrix(0,R,m)

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
    #MSPE estimation using the bootstrap method:
    ##############################################double Bootstrap method########################
    MSE.bot1<-matrix(0,B1,m)
    AMSE.bot2<-matrix(0,B1,m)

    for(b1 in 1:B1){

      random1<-rnorm(m,0,1)
      random2<-rnorm(m,0,1)
      yib.boot1<-rep(0,m)
      mu.boot1<-rep(0,m)
      for(i in 1:n){
        yib.boot1[i]<-as.numeric(X[i,]%*%BETA.hat)+sqrt(A.hat)*random1[i]+SD[i]*random2[i]
        mu.boot1[i]<-as.numeric(X[i,]%*%BETA.hat)+sqrt(A.hat)*random1[i]
      }
      for(i in (n+1):m){
        yib.boot1[i]<-as.numeric(X[i,]%*%BETA.hat)+sqrt(A.hat)*random1[i]+SD[i]*random2[i]
        mu.boot1[i]<-as.numeric(X[i,]%*%BETA.hat)+sqrt(A.hat)*random1[i]
      }
      #####We use the following estimators (based on the PR and EBLUE) approach:

      A.hat.boot1<-(t(yib.boot1)%*%PXI%*%yib.boot1-tr(PXI%*%D1))/(m-p)
      A.hat.boot1<-ifelse(A.hat.boot1<=0,0.0001,A.hat.boot1)


      sumx.tem<-matrix(0,p,p)
      sumy.tem<-matrix(0,p,1)
      tem1=0
      for(i in 1:m){

        sumx.tem<-sumx.tem+( X[i,]%*%t(X[i,] ) ) /drop(A.hat.boot1+SD[i]^2)
        sumy.tem<-sumy.tem+( X[i,]*yib.boot1[i] )/drop(A.hat.boot1+SD[i]^2)
        tem1<-tem1+(A.hat.boot1+D[i])^2
      }
      BETA.hat.boot1<-solve(sumx.tem)%*%sumy.tem
      ############################################MSPE:
      for(i in 1:m){
        MSE.bot1[b1,i]<-((A.hat.boot1*yib.boot1[i])/(A.hat.boot1+D[i]) +(D[i]*X[i,]%*%BETA.hat.boot1)/(A.hat.boot1+D[i])-mu.boot1[i])^2
      }
      #cat(b1)

      #################
      MSE.bot2<-matrix(0,B2,m)

      for(b2 in 1:B2){

        random1<-rnorm(m,0,1)
        random2<-rnorm(m,0,1)
        yib.boot2<-rep(0,m)
        mu.boot2<-rep(0,m)
        for(i in 1:n){
          yib.boot2[i]<-as.numeric(X[i,]%*%BETA.hat.boot1)+sqrt(A.hat.boot1)*random1[i]+SD[i]*random2[i]
          mu.boot2[i]<-as.numeric(X[i,]%*%BETA.hat.boot1)+sqrt(A.hat.boot1)*random1[i]
        }
        for(i in (n+1):m){
          yib.boot2[i]<-as.numeric(X[i,]%*%BETA.hat.boot1)+sqrt(A.hat.boot1)*random1[i]+SD[i]*random2[i]
          mu.boot2[i]<-as.numeric(X[i,]%*%BETA.hat.boot1)+sqrt(A.hat.boot1)*random1[i]
        }
        #####We use the following estimators (based on the PR and EBLUE) approach:

        A.hat.boot2<-(t(yib.boot2)%*%PXI%*%yib.boot2-tr(PXI%*%D1))/(m-p)
        A.hat.boot2<-ifelse(A.hat.boot2<=0,0.0001,A.hat.boot2)


        sumx.tem<-matrix(0,p,p)
        sumy.tem<-matrix(0,p,1)
        tem1=0
        for(i in 1:m){

          sumx.tem<-sumx.tem+( X[i,]%*%t(X[i,] ) ) /drop(A.hat.boot2+SD[i]^2)
          sumy.tem<-sumy.tem+( X[i,]*yib.boot2[i] )/drop(A.hat.boot2+SD[i]^2)
          tem1<-tem1+(A.hat.boot2+D[i])^2
        }
        BETA.hat.boot2<-solve(sumx.tem)%*%sumy.tem
        ############################################MSPE:
        for(i in 1:m){
          MSE.bot2[b2,i]<-((A.hat.boot2*yib.boot2[i])/(A.hat.boot2+D[i]) +(D[i]*X[i,]%*%BETA.hat.boot2)/(A.hat.boot2+D[i])-mu.boot2[i])^2
        }
        #cat(b2)
      }
      #################
      for(i in 1:m){AMSE.bot2[b1,i]<-mean(MSE.bot2[,i])}
    }
    ######################################################
    for(i in 1:m){
      AMSE.bot1.i<-mean(MSE.bot1[,i])
      AAMSE.bot2.i<-mean(AMSE.bot2[,i])

      # mspe.boot1[r,i]<-AMSE.bot1.i
      mspe.boot1[r,i]<-2*AMSE.bot1.i-AAMSE.bot2.i
      if(AMSE.bot1.i < AAMSE.bot2.i){mspe.boot1[r,i]<-AMSE.bot1.i*exp(-(AAMSE.bot2.i-AMSE.bot1.i)/AAMSE.bot2.i )}
      ifelse(mspe.boot1[r,i]>=0,mspe.boot1[r,i],0)

      mspe.boot2[r,i]<-(AMSE.bot1.i)^2/AAMSE.bot2.i
      ifelse(mspe.boot2[r,i]>=0,mspe.boot2[r,i],0)
      # mspe.boot4[r,i]<-AMSE.bot1.i+(1/m)*( cos(m*(AMSE.bot1.i-AAMSE.bot2.i))/sin(m*(AMSE.bot1.i-AAMSE.bot2.i)) )
      #  if(AMSE.bot1.i < AAMSE.bot2.i){mspe.boot4[r,i]<-(AMSE.bot1.i)^2/(AMSE.bot1.i+(1/m)*( cos(m*(AAMSE.bot2.i-AMSE.bot1.i))/sin(m*(AAMSE.bot2.i-AMSE.bot1.i)) ) )}
      #  ifelse(mspe.boot4[r,i]>=0,mspe.boot4[r,i],0)
    }
  }#loop R (simulation runs)
  ######################################################
  BETA.HAT<-apply(BETA.HAT,2,mean)
  A.HAT<-mean(A.HAT)
  MSPE.TRUE<-apply(MSPE,2,mean)
  mspe.boot1<-apply(mspe.boot1,2,mean)
  mspe.boot2<-apply(mspe.boot2,2,mean)

  Par<-c(BETA.HAT, A.HAT)
  MSPE.TRUE.Final<-c(MSPE.TRUE)
  mspe.boot1.Final<-c(mspe.boot1)
  mspe.boot2.Final<-c(mspe.boot2)

  RB.boot1<-(mspe.boot1.Final/MSPE.TRUE.Final)-1
  RB.boot2<-(mspe.boot2.Final/MSPE.TRUE.Final)-1

  list (Par=Par, MSPE.TRUE.Final=MSPE.TRUE.Final, mspe.boot1.Final=mspe.boot1.Final, mspe.boot2.Final=mspe.boot2.Final, RB.boot1=RB.boot1, RB.boot2=RB.boot2) #Result=Result)

}#loop function
