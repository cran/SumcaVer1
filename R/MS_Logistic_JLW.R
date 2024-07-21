#' Model selection MSPE estimation in mixed logistic model using jackknife method.Calculate the model selection mspe of mixed logistic model using jackknife method.
#'
#' @param m number of small areas
#' @param p number of complete model parameters
#' @param ni sample size of each small area
#' @param X covariates for the complete model
#' @param beta regression coefficients of the complete model
#' @param A variance of area-specific random effects
#' @param R number of simulation runs
#'
#' @return Par1: return estimation of model parameters of the complete model
#' @return Par2: return estimation of model parameters of the reduced model
#' @return MSPE: return empirical MSPE of small area predictor
#' @return mspe.JLW: return mspe of small area predictor using the jackknife method
#' @return RB.JLW: return relative bias (RB) of mspe of small area predictor using the jackknife method
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
#' @examples mspe_MS_LOGISTIC_JLW(20,3,2,
#' matrix(runif(60,0,1),nrow=20,byrow=TRUE),c(1,1,1),10,2)
#'
mspe_MS_LOGISTIC_JLW=function(m,p,ni,X,beta,A,R)
{


#=============================================
#Define the constants of the data and model

BETA.HAT<-matrix(0,R,p)
A.HAT<-rep(0,R)
n=sum(ni)
p1=1
BETA.HAT.1<-matrix(0,R,p1)
A.HAT.1<-rep(0,R)

MSPE<-matrix(0,R,m)
mspe.JLW<-matrix(0,R,m)
Bi<-matrix(0,R,m)
Bi.j<-matrix(0,R,m)
mu2.hat.j<-matrix(0,R,m)

BIC.r<-matrix(0,R,2)

comb_with_replacement <- function(n, r){
  return( factorial(n) / (factorial(r) * factorial(n - r)) )
}
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
  res<-glmer(cbind(yib, ni - yib) ~ 1+X[,-1]+(1|county),family="binomial"(link = "logit"), data=data.obs) #, control = glmerControl(check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),

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
               control = glmerControl( check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),
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
    #-----------------------------
    for(i in 1:m){
      for(j in 1:ni){
        integrand1 <- function(x) {exp(j*(x)-ni*log(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))-x^2/(2*A.hat))}
        integrand2 <- function(x) {(exp(as.numeric(X[i,]%*%BETA.hat)+x)/(1+exp(as.numeric(X[i,]%*%BETA.hat)+x)) )*(exp(j*(x)-ni*log(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))-x^2/(2*A.hat)))}
        integrand21 <- function(x) {((exp(as.numeric(X[i,]%*%BETA.hat)+x)/(1+exp(as.numeric(X[i,]%*%BETA.hat)+x)) )^2)*(exp(j*(x)-ni*log(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))-x^2/(2*A.hat)))}
        a1.hat.j<-as.numeric(integrate(integrand2, lower = -3, upper = 3)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
        a2.hat.j<-as.numeric(integrate(integrand21, lower = -3, upper = 3)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])


        integrand11 <- function(x) {
          tem<-exp(as.numeric(X[i,]%*%BETA.hat)+x)/(1+exp(as.numeric(X[i,]%*%BETA.hat)+x))
          return(((tem^j)*(1-tem)^(ni-j))*exp(-0.5*log(2*pi*A.hat)-x^2/(2*A.hat)))
        }
        a.hat.k<-as.numeric(integrate(integrand11, lower = -3.4, upper = 3.4)[1])
        Bi[r,i]<-Bi[r,i]+(comb_with_replacement(ni,j))*(a2.hat.j-a1.hat.j^2)*a.hat.k #(term.MC/K)
      }
    }
    ############################################MSPE:
    for(i in 1:m){
      MSPE[r,i]<-(mu.hat[i]-mu[i])^2
    }
    ###############################################
    #MSPE estimation using the jackknife method:
    BETA.HAT.j<-matrix(0,m,p)
    A.HAT.j<-rep(0,m)
    for(j in 1:m){

      X.tem.j<-matrix(0,m-1,p)
      yi.j<-rep(0,(m-1))

      aa<-1
      for(i in 1:m){
        if (i !=j){
          X.tem.j[aa,]<-X[i,]
          yi.j[aa]<-yib[i]
          aa<-aa+1
        } }

      county.j<-c(1:(m-1))
      data.obs.j <- data.frame(cbind(yi.j, X.tem.j, county.j, ni))

      res.j<-glmer(cbind(yi.j, ni - yi.j) ~ 1+X.tem.j[,-1]+(1|county.j),family="binomial"(link = "logit"), data=data.obs.j,
               control = glmerControl( check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),
                              check.conv.singular = .makeCC(action = "ignore",  tol = 0.05)))

      res1.j<-summary(res.j)
      BETA.hat.j<-c(res1.j$coefficients[1:p,1])
      BETA.HAT.j[j,]<-BETA.hat.j
      A.hat.j<-as.numeric(attr(res1.j$varcor$county.j,"stddev")^2)
      A.hat.j<-ifelse(A.hat.j<=0.0001,0.001,A.hat.j)
      A.HAT.j[j]<-A.hat.j
    } #loop of j
    #---------------------------------
    for(z in 1:m){
      for(i in 1:m) {
        integrand1 <- function(x) {
          exp(yib[z] * (x) - ni * log(1 + exp(as.numeric(X[z,] %*% t(BETA.HAT.j[i,])) + x)) - x^2 / (2 * A.HAT.j[i]))
        }
        integrand2 <- function(x) {
          (exp(as.numeric(X[z,] %*% t(BETA.HAT.j[i,])) + x) / (1 + exp(as.numeric(X[z,] %*% t(BETA.HAT.j[i,])) + x))) *
            (exp(yib[z] * (x) - ni * log(1 + exp(as.numeric(X[z,] %*% t(BETA.HAT.j[i,])) + x)) - x^2 / (2 * A.HAT.j[i])))
        }
        tryCatch({
          mu2.hat.j[r,z] <- mu2.hat.j[r,z] + (as.numeric(integrate(integrand2, lower = -3, upper = 3)$value) /
                                                as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)$value) - mu.hat[z])^2
        }, error = function(e) {
          print(paste("ERROR:", conditionMessage(e)))
        })
      }

    }
    #--------------------------------------------
    for(z in 1:m){
      for(i in 1:m){
        tem<-0
        for(j in 1:ni) {
          integrand1 <- function(x) {
            exp(j * (x) - ni * log(1 + exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x)) - x^2 / (2 * A.HAT.j[i]))
          }
          integrand2 <- function(x) {
            (exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x) / (1 + exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x))) *
              (exp(j * (x) - ni * log(1 + exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x)) - x^2 / (2 * A.HAT.j[i])))
          }
          integrand21 <- function(x) {
            ((exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x) / (1 + exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x))) ^ 2) *
              (exp(j * (x) - ni * log(1 + exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x)) - x^2 / (2 * A.HAT.j[i])))
          }
          tryCatch({
            a1.hat.j <- as.numeric(integrate(integrand2, lower = -3, upper = 3)$value) / as.numeric(integrate(integrand1, lower = -Inf, upper = Inf)$value)
            a2.hat.j <- as.numeric(integrate(integrand21, lower = -3, upper = 3)$value) / as.numeric(integrate(integrand1, lower = -Inf, upper = Inf)$value)
          }, error = function(e) {
            print(paste("ERROR:", conditionMessage(e)))
          })

          integrand11.j <- function(x) {
            tem.j <- exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x) / (1 + exp(as.numeric(X[z, ] %*% t(BETA.HAT.j[i, ])) + x))
            return(((tem.j ^ j) * (1 - tem.j) ^ (ni - j)) * exp(-0.5 * log(2 * pi * A.HAT.j[i]) - x^2 / (2 * A.HAT.j[i])))
          }

          tryCatch({
            a.hat.k.j <- as.numeric(integrate(integrand11.j, lower = -3.4, upper = 3.4)$value)
          }, error = function(e) {
            print(paste("ERROR:", conditionMessage(e)))
          })

          tem <- tem + ((comb_with_replacement(ni, j)) * (a2.hat.j - a1.hat.j^2) * a.hat.k.j) #(term.MC/K))
        }
        #for j
        Bi.j[r,z]<-Bi.j[r,z]+ (tem-Bi[r,z])
      } #for i
    }#for z
    mspe.JLW[r,]<-Bi[r,]-((m-1)/m)*(Bi.j[r,]-mu2.hat.j[r,])
  }#loop for Model 2
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
    #-----------------------------
    for(i in 1:m){
      for(j in 1:ni){
        integrand1 <- function(x) {exp(j*(x)-ni*log(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))-x^2/(2*A.hat.1))}
        integrand2 <- function(x) {(exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)/(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)) )*(exp(j*(x)-ni*log(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))-x^2/(2*A.hat.1)))}
        integrand21 <- function(x) {((exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)/(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)) )^2)*(exp(j*(x)-ni*log(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))-x^2/(2*A.hat.1)))}
        a1.hat.j<-as.numeric(integrate(integrand2, lower = -3, upper = 3)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
        a2.hat.j<-as.numeric(integrate(integrand21, lower = -3, upper = 3)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])

        integrand11 <- function(x) {
          tem<-exp(as.numeric(X1[i,]%*%BETA.hat.1)+x)/(1+exp(as.numeric(X1[i,]%*%BETA.hat.1)+x))
          return(((tem^j)*(1-tem)^(ni-j))*exp(-0.5*log(2*pi*A.hat)-x^2/(2*A.hat)))
        }
        a.hat.k<-as.numeric(integrate(integrand11, lower = -3.4, upper = 3.4)[1])
        Bi[r,i]<-Bi[r,i]+(comb_with_replacement(ni,j))*(a2.hat.j-a1.hat.j^2)*a.hat.k #(term.MC/K)
      }
    }
    ############################################MSPE:
    for(i in 1:m){
      MSPE[r,i]<-(mu.hat[i]-mu[i])^2
    }
    ###############################################
    #MSPE estimation using the jackknife method:
    BETA.HAT.j<-matrix(0,m,p1)
    A.HAT.j<-rep(0,m)

    for(j in 1:m){
      X1.tem.j<-matrix(0,m-1,p1)
      yi.j<-rep(0,(m-1))

      aa<-1
      for(i in 1:m){
        if (i !=j){
          X1.tem.j[aa,]<-X1[i,]
          yi.j[aa]<-yib[i]
          aa<-aa+1
        } }

      county.j<-c(1:(m-1))
      data.obs.j <- data.frame(cbind(yi.j, county.j, ni))

      res.j<-glmer(cbind(yi.j, ni - yi.j) ~ 1+(1|county.j),family="binomial"(link = "logit"), data=data.obs.j,
             control = glmerControl( check.conv.grad = .makeCC("warning", tol = 0.05, relTol = NULL),
                              check.conv.singular = .makeCC(action = "ignore",  tol = 0.05)))


      res1.j<-summary(res.j)
      BETA.hat.j<-c(res1.j$coefficients[1:p1,1])
      BETA.HAT.j[j,]<-BETA.hat.j
      A.hat.j<-as.numeric(attr(res1.j$varcor$county.j,"stddev")^2)

      A.hat.j<-ifelse(A.hat.j<=0.0001,0.001,A.hat.j)
      A.HAT.j[j]<-A.hat.j
    } #loop of j
    #---------------------------------
    for(z in 1:m){
      for(i in 1:m){
        integrand1 <- function(x) {exp(yib[z]*(x)-ni*log(1+exp(as.numeric(X1[z,]%*%t(BETA.HAT.j[i,]))+x))-x^2/(2*A.HAT.j[i]))}
        integrand2 <- function(x) {(exp(as.numeric(X1[z,]%*%t(BETA.HAT.j[i,]))+x)/(1+exp(as.numeric(X1[z,]%*%t(BETA.HAT.j[i,]))+x)) )*(exp(yib[z]*(x)-ni*log(1+exp(as.numeric(X1[z,]%*%t(BETA.HAT.j[i,]))+x))-x^2/(2*A.HAT.j[i])))}
        mu2.hat.j[r,z]<-mu2.hat.j[r,z]+(as.numeric(integrate(integrand2, lower = -3, upper = 3)[1])/as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1]) -mu.hat[z])^2
      }
    }
    #--------------------------------------------
    for(z in 1:m){
      for(i in 1:m){
        tem<-0
        for(j in 1:ni){
          integrand1 <- function(x) {
            exp(j * x - ni * log(1 + exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x)) - x^2 / (2 * A.HAT.j[i]))
          }
          integrand2 <- function(x) {
            (exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x) / (1 + exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x))) *
              (exp(j * x - ni * log(1 + exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x)) - x^2 / (2 * A.HAT.j[i])))
          }
          integrand21 <- function(x) {
            ((exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x) / (1 + exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x)))^2) *
              (exp(j * x - ni * log(1 + exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x)) - x^2 / (2 * A.HAT.j[i])))
          }

          tryCatch({
            a1.hat.j <- as.numeric(integrate(integrand2, lower = -3, upper = 3)[1]) / as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
            a2.hat.j <- as.numeric(integrate(integrand21, lower = -3, upper = 3)[1]) / as.numeric(integrate(integrand1, lower = -3.5, upper = 3.5)[1])
          }, error = function(e) {
            print(paste("ERROR:", conditionMessage(e)))
          })

          integrand11.j <- function(x) {
            tem.j <- exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x) / (1 + exp(as.numeric(X1[z,] %*% t(BETA.HAT.j[i,])) + x))
            return(((tem.j^j) * (1 - tem.j)^(ni - j)) * exp(-0.5 * log(2 * pi * A.HAT.j[i]) - x^2 / (2 * A.HAT.j[i])))
          }

          tryCatch({
            a.hat.k.j <- as.numeric(integrate(integrand11.j, lower = -3.4, upper = 3.4)[1])
          }, error = function(e) {
            print(paste("ERROR:", conditionMessage(e)))
          })

          tem <- tem + ((comb_with_replacement(ni, j)) * (a2.hat.j - a1.hat.j^2) * a.hat.k.j) #(term.MC/K))
        }
        #for j
        Bi.j[r,z]<-Bi.j[r,z]+ (tem-Bi[r,z])
      } #for i
    }#for z
    mspe.JLW[r,]<-Bi[r,]-((m-1)/m)*(Bi.j[r,]-mu2.hat.j[r,])
  }#loop for Model 1
  #####################################################
 }#loop R (simulation runs)
######################################################
BETA.HAT<-apply(BETA.HAT,2,mean)
A.HAT<-mean(A.HAT)
BETA.HAT.1<-apply(BETA.HAT.1,2,mean)
A.HAT.1<-mean(A.HAT.1)

MSPE<-apply(MSPE,2,mean)
mspe.JLW<-apply(mspe.JLW,2,mean)

Par<-c(BETA.HAT, A.HAT)
Par1<-c(BETA.HAT.1, A.HAT.1)

RB.JLW<-(mspe.JLW/MSPE)-1

BIC<-apply(BIC.r,2, mean)

list (Par=Par, Par1=Par1, MSPE=MSPE, mspe.JLW=mspe.JLW, RB.JLW=RB.JLW, BIC=BIC) #Result=Result)

}#loop function
