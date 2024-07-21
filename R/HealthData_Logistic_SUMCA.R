#' MSPE estimation in mixed logistic model (Health Insurance data) using SUMCA method.Calculate the mspe of mixed logistic model (Health Insurance data) using SUMCA method.
#'
#'
#' @param m number of domains
#' @param p number of complete model parameters
#' @param n.new sample size of each domain
#' @param y.new response variable
#' @param Xi covariates
#' @param cum.n.new Cummulative sum of n
#' @param yi.tem response variable for each individual
#' @param X.tem Individual level covariates
#' @param county.tem county
#' @param K number of Monte Carlo for the SUMCA method
#'
#' @return Par: return estimation of model parameters
#' @return Mu.hat: return prediction of domain parameters
#' @return mspe.Sumca: return mspe of small area (domain) predictor using the SUMCA method
#' @return sq.mspe.Sumca: return square root of mspe of small area predictor for non-zero domains using the SUMCA method
#' @export
#'
#' @importFrom lme4 glmer
#' @importFrom stats rnorm
#' @importFrom stats binomial
#' @importFrom stats integrate
#' @importFrom stats rbinom
#' @importFrom lme4 glmerControl
#' @importFrom lme4 .makeCC
#'
#'
#' @examples mspe_LOGISTIC_HealthData_SUMCA(20,3,c(2,1,2,2,1,2,3,1,1,3,1,3,2,3,
#' 3,1,2,1,3,3),c(3,4,2,2,3,3,4,3,4,1,4,1,3,5,4,7,1,3,1,2),
#' matrix(runif(60,0,1),nrow=20,byrow=TRUE),c(2,3,5,7,8,10,13,14,15
#' ,18,19,22,24,27,30,31,33,34,37,40),sample(c(0,1),replace=TRUE,40),
#' matrix(c(runif(40,7,10),runif(40,14,22),runif(40,2,4)),nrow=40,byrow=FALSE),
#' rep(1:20,each=2),10)
#

mspe_LOGISTIC_HealthData_SUMCA=function(m,p,n.new,y.new,Xi,cum.n.new,yi.tem,X.tem,county.tem,K)
{
  #=============================================
 # data=list(X.tem, yi.tem, m, n.new, county.tem, Xi, y.new, cum.n.new)

  ##################################################


  #source("MakingHealthDataReady.R")
#Data:
  #data=list(X.tem, yi.tem, m, n.new, county.tem, Xi, y.new, cum.n.new)
  ###################################
  #Binomial model
 res<-glmer(yi.tem ~ 1+X.tem+(1|county.tem),family="binomial"(link = "logit")) #, control = glmerControl(check.conv.grad     = .makeCC("warning", tol = 0.05, relTol = NULL),
      #                 check.conv.singular = .makeCC(action = "ignore",  tol = 0.05),
      #                check.conv.hess     = .makeCC(action = "warning", tol = 0.05)))

  res1<-summary(res)
  BETA.hat<-c(res1$coefficients[1,1],res1$coefficients[2,1],res1$coefficients[3,1],res1$coefficients[4,1])
 # BETA.hat = c(-4.18613543,  0.61847085,  0.00477239,  0.12482905)
  A.hat<-(res1$optinfo$val[1])^2
  A.hat<-ifelse(A.hat<=0.0001,0.001,A.hat)
 # A.hat
  #[1] 0.004412709
  Par<-c(BETA.hat, A.hat)
  #---------------------------------
  yi<-y.new
  ni<-n.new
  d=length(county.tem)
  #--------------------------------
  Mu.hat<-rep(0,m)
  for(i in 1:m){
    integrand1 <- function(x) {exp(yi[i]*(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat)-x^2/(2*A.hat))}
    integrand2 <- function(x) {(exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)/(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)) )*(exp(yi[i]*(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat)-x^2/(2*A.hat)))}
    Mu.hat[i]<-as.numeric(integrate(integrand2, lower = -Inf, upper = Inf)[1])/as.numeric(integrate(integrand1, lower = -Inf, upper = Inf)[1])
  }
  a1.hat<-Mu.hat
  #--------------------------------------------
  a2.hat<-rep(0,m)
  for(i in 1:m){
    integrand1 <- function(x) {exp(yi[i]*(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat)-x^2/(2*A.hat))}
    integrand2 <- function(x) {((exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)/(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)) )^2)*(exp(yi[i]*(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat)-x^2/(2*A.hat)))}
    a2.hat[i]<-as.numeric(integrate(integrand2, lower = -Inf, upper = Inf)[1])/as.numeric(integrate(integrand1, lower = -Inf, upper = Inf)[1])
  }

  a.hat<-a2.hat-(a1.hat)^2
  ############################################
  #MSPE estimation using the Sumca method:

  bias.MC<-matrix(0,K,m)

  for(k in 1:K){
    random1<-rnorm(m,0,1)
    yi.k<-rep(0,m)
    for(i in 1:m){
      tem.k<-BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+sqrt(A.hat)*random1[i]
      mu.k<-exp(tem.k)/(1+exp(tem.k))
      yi.k[i]<-rbinom(1,ni[i],mu.k)
                  }

    yi.tem.k<-rep(0,d)
    for(i in 1:m){
      if(yi.k[i]>0){
        a4<-ifelse(i==1,1,cum.n.new[i-1]+1)
        b4=a4+yi.k[i]-1
        yi.tem.k[a4:b4]<-rep(1,yi.k[i])
       }

      if(yi.k[i]==0){
        a4<-ifelse(i==1,1,cum.n.new[i-1]+1)
        yi.tem.k[a4:a4]<-0
      }
    }
    ###################################
    #tryCatch({
    yi.tem.k <- rbinom(n = length(county.tem), size = 1, prob = 0.5)
    res.k<-glmer(yi.tem.k~1+X.tem+(1|county.tem),family="binomial"(link = "logit") , control = glmerControl(check.conv.grad     = .makeCC("warning", tol = 0.05, relTol = NULL),
                          check.conv.singular = .makeCC(action = "ignore",  tol = 0.05)))
   #                       check.conv.hess     = .makeCC(action = "warning", tol = 0.05)))
   # }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    res1.k<-summary(res.k)
    BETA.hat.k<-c(res1.k$coefficients[1,1],res1.k$coefficients[2,1],res1.k$coefficients[3,1],res1.k$coefficients[4,1])
    A.hat.k<-(res1.k$optinfo$val[1])^2
    A.hat.k<-ifelse(A.hat.k<=0.0001,0.001,A.hat.k)
    ############################################
    a1.hat.k<-rep(0,m)
    for(i in 1:m) {
      tryCatch({
        integrand31.k <- function(x) {
          exp(yi.k[i] * (BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) -
                ni[i] * log(1 + exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x)) -
                0.5 * log(2 * pi * A.hat.k) - x^2 / (2 * A.hat.k))
        }

        integrand32.k <- function(x) {
          (exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) /
             (1 + exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x))) *
            (exp(yi.k[i] * (BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) -
                   ni[i] * log(1 + exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x)) -
                   0.5 * log(2 * pi * A.hat.k) - x^2 / (2 * A.hat.k)))
        }

        a1.hat.k[i] <- as.numeric(integrate(integrand32.k, lower = -Inf, upper = Inf)[1]) /
          as.numeric(integrate(integrand31.k, lower = -Inf, upper = Inf)[1])
      }, error = function(e) {
        print(paste("ERROR :", conditionMessage(e)))
      })
    }

    #--------------------------------------------
    a2.hat.k<-rep(0,m)
    for(i in 1:m) {
      tryCatch({
        integrand41.k <- function(x) {
          exp(yi.k[i] * (BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) -
                ni[i] * log(1 + exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x)) -
                0.5 * log(2 * pi * A.hat.k) - x^2 / (2 * A.hat.k))
        }

        integrand42.k <- function(x) {
          ((exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) /
              (1 + exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x)))^2) *
            (exp(yi.k[i] * (BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) -
                   ni[i] * log(1 + exp(BETA.hat.k[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x)) -
                   0.5 * log(2 * pi * A.hat.k) - x^2 / (2 * A.hat.k)))
        }

        a2.hat.k[i] <- as.numeric(integrate(integrand42.k, lower = -Inf, upper = Inf)[1]) /
          as.numeric(integrate(integrand41.k, lower = -Inf, upper = Inf)[1])
      }, error = function(e) {
        print(paste("ERROR:", conditionMessage(e)))
      })
    }


    a.hat.k<-a2.hat.k-(a1.hat.k)^2
    #===============================================
    a1.k.hat<-rep(0,m)
    for(i in 1:m){
      tryCatch({
        integrand51.k <- function(x) {
          exp(yi.k[i] * (BETA.hat[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) -
                ni[i] * log(1 + exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x)) -
                0.5 * log(2 * pi * A.hat) - x^2 / (2 * A.hat))
        }
        integrand52.k <- function(x) {
          (exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) /
             (1 + exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x))) *
            (exp(yi.k[i] * (BETA.hat[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x) -
                   ni[i] * log(1 + exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i,1] + BETA.hat.k[3] * Xi[i,2] + BETA.hat.k[4] * Xi[i,3] + x)) -
                   0.5 * log(2 * pi * A.hat) - x^2 / (2 * A.hat)))
        }
        a1.k.hat[i] <- as.numeric(integrate(integrand52.k, lower = -Inf, upper = Inf)[1]) /
          as.numeric(integrate(integrand51.k, lower = -Inf, upper = Inf)[1])
      }, error = function(e) {
        print(paste("ERROR:", conditionMessage(e)))
      })
    }

    #--------------------------------------------
    a2.k.hat<-rep(0,m)
    for(i in 1:m){
      tryCatch({
        integrand11.k <- function(x) {
          exp(yi.k[i] * (BETA.hat[1] + BETA.hat.k[2] * Xi[i, 1] + BETA.hat.k[3] * Xi[i, 2] + BETA.hat.k[4] * Xi[i, 3] + x) -
                ni[i] * log(1 + exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i, 1] + BETA.hat.k[3] * Xi[i, 2] + BETA.hat.k[4] * Xi[i, 3] + x)) -
                0.5 * log(2 * pi * A.hat) - x^2 / (2 * A.hat))
        }
        integrand21.k <- function(x) {
          ((exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i, 1] + BETA.hat.k[3] * Xi[i, 2] + BETA.hat.k[4] * Xi[i, 3] + x) /
              (1 + exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i, 1] + BETA.hat.k[3] * Xi[i, 2] + BETA.hat.k[4] * Xi[i, 3] + x)))^2) *
            (exp(yi.k[i] * (BETA.hat[1] + BETA.hat.k[2] * Xi[i, 1] + BETA.hat.k[3] * Xi[i, 2] + BETA.hat.k[4] * Xi[i, 3] + x) -
                   ni[i] * log(1 + exp(BETA.hat[1] + BETA.hat.k[2] * Xi[i, 1] + BETA.hat.k[3] * Xi[i, 2] + BETA.hat.k[4] * Xi[i, 3] + x)) -
                   0.5 * log(2 * pi * A.hat) - x^2 / (2 * A.hat)))
        }
        a2.k.hat[i] <- as.numeric(integrate(integrand21.k, lower = -Inf, upper = Inf)[1]) / as.numeric(integrate(integrand11.k, lower = -Inf, upper = Inf)[1])
      }, error = function(e) {
        print(paste("ERROR:", conditionMessage(e)))
      })
    }


    a.k.hat<-a2.k.hat-(a1.k.hat)^2
    #---------------------------------------------
    bias.MC[k,]<-a.k.hat-a.hat.k
   # cat(k)
  }#loop K

  bias.co<-rep(0,m)
  for(i in 1:m){
    bias.co[i]<-sum(bias.MC[,i],na.rm=TRUE)/(K-sum(is.na(bias.MC[,i])))
  }

  mspe.Sumca<-a.hat+bias.co
  #########################
  #square root of mspe for non-zero domains:
  sq.mspe.Sumca<-sqrt(mspe.Sumca[c(1,2,3,4,6,8,9,12,16,18,24,27,32,36,48,54,64,72,96)])

  list (Par=Par, Mu.hat=Mu.hat, mspe.Sumca=mspe.Sumca, sq.mspe.Sumca=sq.mspe.Sumca) #Result=Result)

}#loop function
