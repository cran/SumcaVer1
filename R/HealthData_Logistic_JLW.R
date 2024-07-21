#' MSPE estimation in mixed logistic model (Health Insurance data) using jackknife method.Calculate the mspe of mixed logistic model (Health Insurance data) using jackknife method.
#'
#'
#' @param m number of domains
#' @param p number of complete model parameters
#' @param n.new sample size of each domain
#' @param y.new response variable
#' @param Xi covariates for each domain
#' @param yi.tem response variable for each individual
#' @param cum.n.new Cummulative sum of n
#' @param county.tem county
#' @param X.tem Individual level covariates
#'
#' @return Par: return estimation of model parameters
#' @return Mu.hat: return prediction of domain parameters
#' @return mspe.JLW: return mspe of small area (domain) predictor using the jackknife method
#' @return sq.mspe.JLW: return square root of mspe of small area predictor for non-zero domains using the jackknife method
#' @export
#'
#' @importFrom lme4 glmer
#' @importFrom stats binomial
#' @importFrom stats integrate
#' @importFrom lme4 glmerControl
#' @importFrom lme4 .makeCC
#'
#' @examples mspe_LOGISTIC_HealthData_JLW(20,3,c(2,1,2,2,1,2,3,1,1,3,1,3,2,3,3,
#' 1,2,1,3,3),c(3,4,2,2,3,3,4,3,4,1,4,1,3,5,4,7,1,3,1,2),
#' matrix(runif(60,0,1),nrow=20,byrow=TRUE),sample(c(0,1),replace=TRUE,40),
#' c(2,3,5,7,8,10,13,14,15,18,19,22,24,27,30,31,33,34,37,40),rep(1:20,each=2),
#' matrix(c(runif(40,7,10),runif(40,14,22),runif(40,2,4)),nrow=40,byrow=FALSE))
#'
#'
mspe_LOGISTIC_HealthData_JLW=function(m,p,n.new,y.new,Xi,yi.tem,cum.n.new,county.tem,X.tem)
{
  #=============================================
 # data=list(X.tem, yi.tem, m, n.new, county.tem, Xi, y.new, cum.n.new)
  ave.inc.tem = X.tem[,1]
  ave.educ.tem = X.tem[,2]
  ave.fm.size.tem = X.tem[,3]
  ##################################################
  Bi<-rep(0,m)
  Bi.j<-rep(0,m)
  mu2.hat.j<-rep(0,m)

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
  a1.hat<-Mu.hat
  #--------------------------------------------
  a2.hat<-rep(0,m)
  for(i in 1:m){
    integrand1 <- function(x) {exp(yi[i]*(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat)-x^2/(2*A.hat))}
    integrand2 <- function(x) {((exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)/(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)) )^2)*(exp(yi[i]*(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x)-ni[i]*log(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3]+x))-0.5*log(2*pi*A.hat)-x^2/(2*A.hat)))}
    a2.hat[i]<-as.numeric(integrate(integrand2, lower = -Inf, upper = Inf)[1])/as.numeric(integrate(integrand1, lower = -Inf, upper = Inf)[1])
  }
  #-------------------------------------------
  for(i in 1:m){
    integrand12 <- function(x) {
      tem<-exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2]+ BETA.hat[4]*Xi[i,3] +x)/(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3] +x))
      return((tem^2)*exp(-0.5*log(2*pi*A.hat)-x^2/(2*A.hat)))
    }
    a1.i2<-as.numeric(integrate(integrand12, lower = -3.4, upper = 3.4)[1])

    integrand11 <- function(x) {
      tem<-exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2]+ BETA.hat[4]*Xi[i,3] +x)/(1+exp(BETA.hat[1]+BETA.hat[2]*Xi[i,1] +BETA.hat[3]*Xi[i,2] +BETA.hat[4]*Xi[i,3] +x))
      return((tem)*exp(-0.5*log(2*pi*A.hat)-x^2/(2*A.hat)))
    }
    a1.i1<-as.numeric(integrate(integrand11, lower = -3.4, upper = 3.4)[1])

    Bi[i]<-a1.i2-(a1.i1)^2
  }
  Bi[is.na(Bi)] <- mean(Bi,na.rm=TRUE)
  ############################################
  #MSPE estimation using the jackknife method:

  p=4
  BETA.HAT.j<-matrix(0,m,p)
  A.HAT.j<-rep(0,m)
  # mu.hat.j<-rep(0,m)

  #print(r)
  for(j in 1:m){
    # tem<-length(ni[j])
    yi.tem.j<-rep(0,(d-ni[j])) #(m-1)*ni)
    county.tem.j<-rep(0,(d-ni[j])) #ni*(m-1))
    ave.fm.size.tem.j<-matrix(0,(d-ni[j])) #ni*(m-1))
    ave.educ.tem.j<-matrix(0,(d-ni[j])) #ni*(m-1))
    ave.inc.tem.j<-matrix(0,(d-ni[j])) #ni*(m-1))

    if(j==1) {yi.tem.j<-yi.tem[-(1:cum.n.new[j])]
    county.tem.j<-county.tem[-(1:cum.n.new[j])]
    ave.fm.size.tem.j<-ave.fm.size.tem[-(1:cum.n.new[j])]
    ave.educ.tem.j<-ave.educ.tem[-(1:cum.n.new[j])]
    ave.inc.tem.j<-ave.inc.tem[-(1:cum.n.new[j])]
    }
    if(j!=1){yi.tem.j<-yi.tem[-((cum.n.new[j-1]+1):cum.n.new[j])]
    county.tem.j<-county.tem[-((cum.n.new[j-1]+1):cum.n.new[j])]
    ave.fm.size.tem.j<-ave.fm.size.tem[-((cum.n.new[j-1]+1):cum.n.new[j])]
    ave.educ.tem.j<-ave.educ.tem[-((cum.n.new[j-1]+1):cum.n.new[j])]
    ave.inc.tem.j<-ave.inc.tem[-((cum.n.new[j-1]+1):cum.n.new[j])]
    }

    xi.tem.j<-cbind(ave.inc.tem.j, ave.educ.tem.j, ave.fm.size.tem.j)

    res.j<-glmer(yi.tem.j ~ 1+xi.tem.j+(1|county.tem.j),family="binomial"(link = "logit") ,
                 control = glmerControl(check.conv.grad     = .makeCC("warning", tol = 0.05, relTol = NULL),
                          check.conv.singular = .makeCC(action = "ignore",  tol = 0.05)))

    res1.j<-summary(res.j)
    BETA.hat.j<-c(res1.j$coefficients[1,1],res1.j$coefficients[2,1],res1.j$coefficients[3,1],res1.j$coefficients[4,1])
    BETA.HAT.j[j,]<-BETA.hat.j
    A.hat.j<-(res1.j$optinfo$val[1])^2
    A.hat.j<-ifelse(A.hat.j<=0.0001,0.001,A.hat.j)
    A.HAT.j[j]<-A.hat.j
  } #loop of j
  #---------------------------------
  for(z in 1:m){
    for(i in 1:m){
      integrand1 <- function(x) {exp(yi[z]*(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3]+x)-ni[z]*log(1+exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3]+x))-0.5*log(2*pi*A.HAT.j[i])-x^2/(2*A.HAT.j[i]))}
      integrand2 <- function(x) {(exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3]+x)/(1+exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3]+x)) )*(exp(yi[z]*(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3]+x)-ni[z]*log(1+exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3]+x))-0.5*log(2*pi*A.HAT.j[i])-x^2/(2*A.HAT.j[i])))}
      mu2.hat.j[z]<-mu2.hat.j[z]+(as.numeric(integrate(integrand2, lower = -3, upper = 3)[1])/as.numeric(integrate(integrand1, lower = -Inf, upper = Inf)[1]) -Mu.hat[z])^2
    }
  }
  #--------------------------------------------
  for(z in 1:m){
    for(i in 1:m){
      tem<-0
      integrand12.j <- function(x) {
        temj2<-exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3] +x)/(1+exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3] +x))
        return((temj2^2)*exp(-0.5*log(2*pi*A.HAT.j[i])-x^2/(2*A.HAT.j[i])))
      }
      a1.i2.j<-as.numeric(integrate(integrand12.j, lower = -3.4, upper = 3.4)[1])

      integrand11.j <- function(x) {
        temj1<-exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3] +x)/(1+exp(BETA.HAT.j[i,1]+BETA.HAT.j[i,2]*Xi[z,1]+BETA.HAT.j[i,3]*Xi[z,2]+BETA.HAT.j[i,4]*Xi[z,3] +x))
        return((temj1)*exp(-0.5*log(2*pi*A.HAT.j[i])-x^2/(2*A.HAT.j[i])))
      }
      a1.i1.j<-as.numeric(integrate(integrand11.j, lower = -3.4, upper = 3.4)[1])

      tem<-a1.i2.j-(a1.i1.j)^2
      Bi.j[z]<-Bi.j[z]+ (tem-Bi[z])
    }
  }
  Bi.j[is.na(Bi.j)] <- mean(Bi.j,na.rm=TRUE)
  #---------------------------
  mspe.JLW<-Bi-((m-1)/m)*(Bi.j-mu2.hat.j)
  #########################
  #square root of mspe for non-zero domains:
  mspe.JLW[16]<-Bi[16] #as we get a negative value for mspe.JLW[16]
  sq.mspe.JLW<-sqrt(mspe.JLW[c(1,2,3,4,6,8,9,12,16,18,24,27,32,36,48,54,64,72,96)])

  list (Par=Par, Mu.hat=Mu.hat, mspe.JLW=mspe.JLW, sq.mspe.JLW=sq.mspe.JLW) #Result=Result)

}#loop function
