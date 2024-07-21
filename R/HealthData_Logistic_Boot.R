#' MSPE estimation in mixed logistic model (Health Insurance data) using bootstrap method. Calculate the mspe of mixed logistic model (Health Insurance data) using bootstrap method.
#'
#'
#' @param m number of domains
#' @param p number of complete model parameters
#' @param n.new sample size of each domain
#' @param y.new response variable
#' @param cum.n.new Cummulaticve sum of n
#' @param Xi covariates
#' @param yi.tem response variable for each individual
#' @param X.tem Individual level covariates
#' @param county.tem county
#' @param B number of bootstrap iterations
#'
#' @return Par: return estimation of model parameters
#' @return Mu.hat: return prediction of domain parameters
#' @return mspe.boot: return mspe of small area (domain) predictor using the bootstrap method
#' @return sq.mspe.boot: return square root of mspe of small area predictor for non-zero domains using the bootstrap method
#' @export
#'
#'
#' @importFrom lme4 glmer
#' @importFrom stats rnorm
#' @importFrom stats binomial
#' @importFrom stats integrate
#' @importFrom stats rbinom
#' @importFrom lme4 glmerControl
#' @importFrom lme4 .makeCC
#'
#' @examples mspe_LOGISTIC_HealthData_BOOT(20,3,c(2,1,2,2,1,2,3,1,1,3,1,3,2,3,3,
#' 1,2,1,3,3),c(3,4,2,2,3,3,4,3,4,1,4,1,3,5,4,7,1,3,1,2),c(2,3,5,7,8,10,13,14,15
#' ,18,19,22,24,27,30,31,33,34,37,40),
#' matrix(runif(60,0,1),nrow=20,byrow=TRUE),sample(c(0,1),replace=TRUE,40),
#' matrix(c(runif(40,7,10),runif(40,14,22),runif(40,2,4)),nrow=40,byrow=FALSE),
#' rep(1:20,each=2),10)
#'
#'
mspe_LOGISTIC_HealthData_BOOT=function(m,p,n.new,y.new,cum.n.new,Xi,yi.tem,X.tem,county.tem,B)
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
 # BETA.hat
 # [1] -4.18613543  0.61847085  0.00477239  0.12482905
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
  #a1.hat<-Mu.hat
  #--------------------------------------------
  #First-phase bootstrap method:
  mspe.boot1<-matrix(0,B,m)
  for(b in 1:B){
    random1<-rnorm(m,0,1)
    yib.boot<-rep(0,m)
    mu.boot<-rep(0,m)
    for(i in 1:m){
      tem.b<-BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+sqrt(A.hat)*random1[i]
      mu.boot[i]<-exp(tem.b)/(1+exp(tem.b))
      yib.boot[i]<-rbinom(1,ni[i],mu.boot[i])
    }

    yi.tem.b<-rep(0,d)
    for(i in 1:m){
      if(yib.boot[i]>0){
        a4<-ifelse(i==1,1,cum.n.new[i-1]+1)
        b4=a4+yib.boot[i]-1
        yi.tem.b[a4:b4]<-rep(1,yib.boot[i])
      }

      if(yib.boot[i]==0){
        a4<-ifelse(i==1,1,cum.n.new[i-1]+1)
        yi.tem.b[a4:a4]<-0
      }
    }
    ###################################
    yi.tem.b <- rbinom(n = length(county.tem), size = 1, prob = 0.5)

    res.b<-glmer(yi.tem.b~1+X.tem+(1|county.tem),family="binomial"(link = "logit"),control = glmerControl(check.conv.grad     = .makeCC("warning", tol = 0.05, relTol = NULL),
                 check.conv.singular = .makeCC(action = "ignore",  tol = 0.05)))

    res1.b<-summary(res.b)
    BETA.hat.b<-c(res1.b$coefficients[1,1],res1.b$coefficients[2,1],res1.b$coefficients[3,1],res1.b$coefficients[4,1])
    A.hat.b<-(res1.b$optinfo$val[1])^2
    A.hat.b<-ifelse(A.hat.b<=0.0001,0.001,A.hat.b)
    ############################################
    mu.hat.b<-rep(0,m)
    for(i in 1:m){
      tryCatch({
        integrand31.b <- function(x) {
          exp(yib.boot[i]*(BETA.hat.b[1]+BETA.hat.b[2]*Xi[i,1] +BETA.hat.b[3]*Xi[i,2] +BETA.hat.b[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat.b[1]+BETA.hat.b[2]*Xi[i,1] +BETA.hat.b[3]*Xi[i,2] +BETA.hat.b[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat.b)-x^2/(2*A.hat.b))
        }
        integrand32.b <- function(x) {
          (exp(BETA.hat.b[1]+BETA.hat.b[2]*Xi[i,1] +BETA.hat.b[3]*Xi[i,2] +BETA.hat.b[4]*Xi[i,3]+x)/(1+exp(BETA.hat.b[1]+BETA.hat.b[2]*Xi[i,1] +BETA.hat.b[3]*Xi[i,2] +BETA.hat.b[4]*Xi[i,3]+x)) )*(exp(yib.boot[i]*(BETA.hat.b[1]+BETA.hat.b[2]*Xi[i,1] +BETA.hat.b[3]*Xi[i,2] +BETA.hat.b[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat.b[1]+BETA.hat.b[2]*Xi[i,1] +BETA.hat.b[3]*Xi[i,2] +BETA.hat.b[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat.b)-x^2/(2*A.hat.b)))
        }
        mu.hat.b[i] <- as.numeric(integrate(integrand32.b, lower = -Inf, upper = Inf)[1]) / as.numeric(integrate(integrand31.b, lower = -Inf, upper = Inf)[1])
      }, error=function(e){
        print(paste("ERROR :", conditionMessage(e)))
      })
    }

    #--------------------------------------------
    mspe.boot1[b,]<-(mu.hat.b-mu.boot)^2
   # cat(b)
  }#loop of B
  mspe.boot<-apply(mspe.boot1,2,mean)
  #########################
  #square root of mspe for non-zero domains:
  sq.mspe.boot<-sqrt(mspe.boot[c(1,2,3,4,6,8,9,12,16,18,24,27,32,36,48,54,64,72,96)])

  list (Par=Par, Mu.hat=Mu.hat, mspe.boot=mspe.boot, sq.mspe.boot=sq.mspe.boot) #Result=Result)

}#loop function
