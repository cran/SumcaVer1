#' Model selection MSPE estimation in mixed logistic model using SUMCA method. Calculate the model selection mspe of mixed logistic model using SUMCA method.
#'
#' @param m number of small areas
#' @param p number of complete model parameters
#' @param ni sample size of each small area
#' @param X covariates for the complete model
#' @param beta regression coefficients of the complete model
#' @param A variance of area-specific random effects
#' @param R number of simulation runs
#' @param K number of Monte Carlo for the SUMCA method
#'
#' @return Par1: return estimation of model parameters of the complete model
#' @return Par2: return estimation of model parameters of the reduced model
#' @return MSPE: return empirical MSPE of small area predictor
#' @return mspe.Sumca: return mspe of small area predictor using the SUMCA method
#' @return RB.SUMCA: return relative bias (RB) of mspe of small area predictor using the SUMCA method
#' @return BIC: return BIC of the complete and reduced models
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
#' @examples mspe_MS_LOGISTIC_SUMCA(20,3,2,matrix(runif(60,0,1),nrow=20,byrow=TRUE),c(1,1,1),10,5,5)
#'
mspe_MS_LOGISTIC_SUMCA=function(m,p,ni,X,beta,A,K,R)
{

#=============================================
#Define the constants of the data and model

BETA.HAT<-matrix(0,R,p)
A.HAT<-rep(0,R)

p1=1
BETA.HAT.1<-matrix(0,R,p1)
A.HAT.1<-rep(0,R)

MSPE<-matrix(0,R,m)
mspe.Sumca<-matrix(0,R,m)

BIC.r<-matrix(0,R,2)

ai.hat<-matrix(0,R,m)
X1<-matrix(c(rep(1,m)),m,1)
#-----------------------------------
for(r in 1:R){

  random1<-rnorm(m,0,1)
  yib<-rep(0,m)
  mu<-rep(0,m)

  for(i in 1:m){
    tem<-as.numeric(X[i,]%*%beta)+sqrt(A)*random1[i]
    mu[i]<-exp(tem)/(1+exp(tem))
    yib[i]<-rbinom(1,ni,mu[i])
  }

  ######################################################################################
  #Binomial model (Model 2: complete model)

  county<-c(1:m)
  data.obs <- data.frame(cbind(yib, X, county, ni))
  #tryCatch({
  res<-glmer(cbind(yib, ni - yib) ~ 1+X[,-1]+(1|county),family="binomial"(link = "logit"), data=data.obs,
   control = glmerControl( check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),
              check.conv.singular = .makeCC(action = "ignore",  tol = 0.05)))
  #}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  r=2
  res1<-summary(res)
  BETA.hat<-c(res1$coefficients[1:p,1])
  BETA.HAT[r,]<-BETA.hat
  A.hat<-as.numeric(attr(res1$varcor$county,"stddev")^2)
  A.HAT[r]<-A.hat
  A.hat<-ifelse(A.hat<=0,0.0001,A.hat)

  BIC.r[r,2]<-as.numeric(res1$AICtab)[2]

  ######################################################################################
  #Binomial model (Model 1: reduced model)

  data.obs.1 <- data.frame(cbind(yib, county, ni))
  res.1<-glmer(cbind(yib, ni - yib) ~ 1+(1|county),family="binomial"(link = "logit"), data=data.obs.1,
             control = glmerControl(check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),
                                    check.conv.singular = .makeCC(action = "ignore",  tol = 0.05)))

  res1.1<-summary(res.1)
  BETA.hat.1<-c(res1.1$coefficients[1:p1,1])
  BETA.HAT.1[r,]<-BETA.hat.1
  A.hat.1<-as.numeric(attr(res1.1$varcor$county,"stddev")^2)
  A.HAT.1[r]<-A.hat.1
  A.hat.1<-ifelse(A.hat.1<=0,0.0001,A.hat.1)

  BIC.r[r,1]<-as.numeric(res1.1$AICtab)[2]

  num<-which.min(BIC.r[r,])  #WHICH MODEL IS SELECTED BASED ON BIC:
  #---------------------------------
  if(num==2){  ###################If Model 2 is selected:
    mu.hat<-rep(0,m)
    for(i in 1:m){
      integrand1 <- function(x) {exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))-x^2/(2*A.hat))}
      integrand2 <- function(x) {(exp(as.numeric(X[i,]%*%BETA.hat)+x)/(1+exp(as.numeric(X[i,]%*%BETA.hat)+x)) )*(exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))-x^2/(2*A.hat)))}
      mu.hat[i]<-as.numeric(integrate(integrand2, lower = -3.5, upper = 3.5)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
    }
    a1.hat<-mu.hat
    #--------------------------------------------
    a2.hat<-rep(0,m)
    for(i in 1:m){
      integrand1 <- function(x) {exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))-x^2/(2*A.hat))}
      integrand2 <- function(x) {((exp(as.numeric(X[i,]%*%BETA.hat)+x)/(1+exp(as.numeric(X[i,]%*%BETA.hat)+x)) )^2)*(exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))-x^2/(2*A.hat)))}
      a2.hat[i]<-as.numeric(integrate(integrand2, lower = -3.5, upper = 3.5)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
    }
    a.hat<-a2.hat-(a1.hat)^2

    ai.hat[r,]<-a.hat
    ############################################MSPE:
    for(i in 1:m){
      MSPE[r,i]<-(mu.hat[i]-mu[i])^2
    }
    ###############################################
    #MSPE estimation using the Sumca method:

    bias.MC <- matrix(0, K, m)

    for (k in 1:K) {
      random1 <- rnorm(m, 0, 1)
      yib.k <- rep(0, m)
      for (i in 1:m) {
        tem.k <- as.numeric(X[i,] %*% BETA.hat) + sqrt(A.hat) * random1[i]
        mu.k <- exp(tem.k) / (1 + exp(tem.k))
        yib.k[i] <- rbinom(1, ni, mu.k)
      }

      data.obs.k <- data.frame(cbind(yib.k, X, county, ni))
      res.k <- glmer(cbind(yib.k, ni - yib.k) ~ 1 + X[,-1] + (1 | county), family = "binomial"(link = "logit"), data = data.obs.k,
                     control = glmerControl(check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),
                                            check.conv.singular = .makeCC(action = "ignore", tol = 0.05)))

      res1.k <- summary(res.k)
      BETA.hat.k <- c(res1.k$coefficients[1:p, 1])
      A.hat.k <- as.numeric(attr(res1.k$varcor$county, "stddev")^2)
      A.hat.k <- ifelse(A.hat.k <= 0.0001, 0.001, A.hat.k)
      ############################################
      a1.hat.k <- rep(0, m)
      for (i in 1:m) {
        tryCatch({
          integrand31.k <- function(x) {exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k))}
          integrand32.k <- function(x) {(exp(as.numeric(X[i,] %*% BETA.hat.k) + x) / (1 + exp(as.numeric(X[i,] %*% BETA.hat.k) + x))) * (exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k)))}
          a1.hat.k[i] <- as.numeric(integrate(integrand32.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand31.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }
      #--------------------------------------------
      a2.hat.k <- rep(0, m)
      for (i in 1:m) {
        tryCatch({
          integrand41.k <- function(x) {exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k))}
          integrand42.k <- function(x) {((exp(as.numeric(X[i,] %*% BETA.hat.k) + x) / (1 + exp(as.numeric(X[i,] %*% BETA.hat.k) + x)))^2) * (exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k)))}
          a2.hat.k[i] <- as.numeric(integrate(integrand42.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand41.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }

      a.hat.k <- a2.hat.k - (a1.hat.k)^2
      #===============================================
      a1.k.hat <- rep(0, m)
      for (i in 1:m) {
        tryCatch({
          integrand51.k <- function(x) {exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat))}
          integrand52.k <- function(x) {(exp(as.numeric(X[i,] %*% BETA.hat) + x) / (1 + exp(as.numeric(X[i,] %*% BETA.hat) + x))) * (exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat)))}
          a1.k.hat[i] <- as.numeric(integrate(integrand52.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand51.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }
      #--------------------------------------------
      a2.k.hat <- rep(0, m)
      for (i in 1:m) {
        tryCatch({
          integrand11.k <- function(x) {exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat))}
          integrand21.k <- function(x) {((exp(as.numeric(X[i,] %*% BETA.hat) + x) / (1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)))^2) * (exp(yib.k[i] * (x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat)))}
          a2.k.hat[i] <- as.numeric(integrate(integrand21.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand11.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }

      a.k.hat <- a2.k.hat - (a1.k.hat)^2
      #---------------------------------------------
      bias.MC[k, ] <- a.k.hat - a.hat.k
    }
    #loop K

    bias.co<-apply(bias.MC,2,mean, na.rm = TRUE)
    mspe.Sumca[r,]<-a.hat+bias.co
  }#loop for Model 2 (complete model)
  ##################################################
  if(num==1){  ###################If Model 1 is selected:
    mu.hat<-rep(0,m)
    for(i in 1:m){
      integrand1 <- function(x) {exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))-x^2/(2*A.hat.1))}
      integrand2 <- function(x) {(exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)/(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)) )*(exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))-x^2/(2*A.hat.1)))}
      mu.hat[i]<-as.numeric(integrate(integrand2, lower = -3.5, upper = 3.5)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
    }
    a1.hat<-mu.hat
    #--------------------------------------------
    a2.hat<-rep(0,m)
    for(i in 1:m){
      integrand1 <- function(x) {exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))-x^2/(2*A.hat.1))}
      integrand2 <- function(x) {((exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)/(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)) )^2)*(exp(yib[i]*(x)-ni*log(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))-x^2/(2*A.hat.1)))}
      a2.hat[i]<-as.numeric(integrate(integrand2, lower = -3.5, upper = 3.5)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
    }
    a.hat<-a2.hat-(a1.hat)^2

    ai.hat[r,]<-a.hat
    ############################################MSPE:
    for(i in 1:m){
      MSPE[r,i]<-(mu.hat[i]-mu[i])^2
    }
    ###############################################
    #MSPE estimation using the Sumca method:
    bias.MC <- matrix(0, K, m)

    for(k in 1:K){
      random1 <- rnorm(m, 0, 1)
      yib.k <- rep(0, m)
      for(i in 1:m){
        tem.k <- as.numeric(X1[i,] %*% BETA.hat.1) + sqrt(A.hat.1) * random1[i]
        mu.k <- exp(tem.k) / (1 + exp(tem.k))
        yib.k[i] <- rbinom(1, ni, mu.k)
      }
      data.obs.k <- data.frame(cbind(yib.k, county, ni))
      res.k <- glmer(cbind(yib.k, ni - yib.k) ~ 1 + (1 | county),
                     family = "binomial"(link = "logit"),
                     data = data.obs.k,
                     control = glmerControl(
                       check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),
                       check.conv.singular = .makeCC(action = "ignore", tol = 0.05)
                     ))

      res1.k <- summary(res.k)
      BETA.hat.k <- c(res1.k$coefficients[1:p1, 1])
      A.hat.k <- as.numeric(attr(res1.k$varcor$county, "stddev")^2)
      A.hat.k <- ifelse(A.hat.k <= 0.0001, 0.001, A.hat.k)
      ############################################
      a1.hat.k <- rep(0, m)
      for(i in 1:m){
        tryCatch({
          integrand31.k <- function(x) {exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X1[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k))}
          integrand32.k <- function(x) {(exp(as.numeric(X1[i,] %*% BETA.hat.k) + x) / (1 + exp(as.numeric(X1[i,] %*% BETA.hat.k) + x))) * (exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X1[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k)))}
          a1.hat.k[i] <- as.numeric(integrate(integrand32.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand31.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }
      #--------------------------------------------
      a2.hat.k <- rep(0, m)
      for(i in 1:m){
        tryCatch({
          integrand41.k <- function(x) {exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X1[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k))}
          integrand42.k <- function(x) {((exp(as.numeric(X1[i,] %*% BETA.hat.k) + x) / (1 + exp(as.numeric(X1[i,] %*% BETA.hat.k) + x)))^2) * (exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X1[i,] %*% BETA.hat.k) + x)) - x^2 / (2 * A.hat.k)))}
          a2.hat.k[i] <- as.numeric(integrate(integrand42.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand41.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }

      a.hat.k <- a2.hat.k - (a1.hat.k)^2
      #===============================================
      a1.k.hat <- rep(0, m)
      for(i in 1:m){
        tryCatch({
          integrand51.k <- function(x) {exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat))}
          integrand52.k <- function(x) {(exp(as.numeric(X[i,] %*% BETA.hat) + x) / (1 + exp(as.numeric(X[i,] %*% BETA.hat) + x))) * (exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat)))}
          a1.k.hat[i] <- as.numeric(integrate(integrand52.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand51.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }
      #--------------------------------------------
      a2.k.hat <- rep(0, m)
      for(i in 1:m){
        tryCatch({
          integrand11.k <- function(x) {exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat))}
          integrand21.k <- function(x) {((exp(as.numeric(X[i,] %*% BETA.hat) + x) / (1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)))^2) * (exp(yib.k[i]*(x) - ni * log(1 + exp(as.numeric(X[i,] %*% BETA.hat) + x)) - x^2 / (2 * A.hat)))}
          a2.k.hat[i] <- as.numeric(integrate(integrand21.k, lower = -3.5, upper = 3.5)[1]) / as.numeric(integrate(integrand11.k, lower = -3.5, upper = 3.5)[1])
        }, error = function(e) {print(paste("ERROR :", conditionMessage(e)))})
      }

      a.k.hat <- a2.k.hat - (a1.k.hat)^2
      #---------------------------------------------
      bias.MC[k,] <- a.k.hat - a.hat.k
    }
    #loop K

    bias.co<-apply(bias.MC,2,mean, na.rm = TRUE)

    mspe.Sumca[r,]<-a.hat+bias.co
  }#loop for Model 1
  #####################################################
 }#loop R (simulation runs)
######################################################
BETA.HAT<-apply(BETA.HAT,2,mean)
A.HAT<-mean(A.HAT)
BETA.HAT.1<-apply(BETA.HAT.1,2,mean)
A.HAT.1<-mean(A.HAT.1)

MSPE<-apply(MSPE,2,mean)
mspe.Sumca<-apply(mspe.Sumca,2,mean)

Par<-c(BETA.HAT, A.HAT)
Par1<-c(BETA.HAT.1, A.HAT.1)

RB.SUMCA<-(mspe.Sumca/MSPE)-1

BIC<-apply(BIC.r,2, mean)

list (Par=Par, Par1=Par1, MSPE=MSPE, mspe.Sumca=mspe.Sumca, RB.SUMCA=RB.SUMCA, BIC=BIC) #Result=Result)

}#loop function

