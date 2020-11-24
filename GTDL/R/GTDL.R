#'@name GTDL
#'@title The Distribution GTDL
#'
#'@description Density, survival function, failure function and random generation for the GTDL distribution.
#'
#'@param t non-negative random variable representing the failure time and leave the snapshot failure rate, or danger.
#'@param alpha,gamma scalars.
#'@param lambda non-negative.
#'@param n number of observations. If length(n) > 1, the length is taken to be the number required.
#'
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}
#'
#'
#'
#'@references
#'
#' Mackenzie,G. ,(2016),Regression Models for Survival Data: The Generalized Time-Dependent Logistic Family
#'
#'@examples
#'t <- seq(0,20,by = 0.01)
#'lambda <- 1
#'alpha <- -0.05
#'gamma <- -1
#'denGTDL<- dGTDL(t,lambda,alpha,gamma,log = FALSE)
#'plot(x = t,y = denGTDL)

#'@rdname GTDL
#'@export

dGTDL<-function(t,lambda,alpha,gamma,log = FALSE){
  if(lambda<=0){
    stop("lambda is not greater than zero")
  }

  if(alpha==0){
    stop("alpha is no different than zero")
  }

  if(gamma==0){
    stop("gamma is no different than zero")
  }

  part1<-(lambda*exp(t*alpha+gamma)/(1+exp(t*alpha+gamma)))*((1+exp(t*alpha+gamma))/(1+exp(gamma)))^(-lambda/alpha)
  if(log == FALSE){
  return(part1)
  }
  if(log == TRUE){
  return(log(part1))
  }
}

#'@rdname GTDL
#'@export

hGTDL<-function(t,lambda,alpha,gamma,log = FALSE){

  if(lambda<=0){
    stop("lambda is not greater than zero")
  }

  if(gamma==0){
    stop("gamma is no different than zero")
  }

  return(lambda*exp(t*alpha+gamma)/(1+exp(t*alpha+gamma)))

}

#'@rdname GTDL
#'@export

sGTDL<- function(t,lambda,alpha,gamma,log = FALSE){
  return(dGTDL(t,lambda,alpha,gamma)/hGTDL(t,lambda,alpha,gamma))
}

#'@rdname GTDL
#'@export

rGTDL<-function(n,lambda,alpha,gamma){
  u<-runif(n)
  t<-(1/alpha)*(log((1+exp(gamma))*(1-u)^(-alpha/lambda)-1)-gamma)
  return(t)
}
#'@rdname GTDL
#'@export
rCENGTDL<-function(n,param,pcensura){
  lambda<-param[1]
  alpha<-param[2]
  gamma<-param[3]
  ncensura<-n*pcensura
  amostra<-data.frame(t=0,censura = 0)

  while(sum(amostra$censura==0)<ncensura+1){
    t1<-rGTDL(1,lambda,alpha,gamma)
    t2<-rGTDL(1,lambda,alpha,gamma)
    if(min(t1,t2)==t1){
      aux<-c(t1,1)
      amostra<-rbind(amostra,aux)
    }
    if(min(t1,t2)==t2){
      aux<-c(t2,0)
      amostra<-rbind(amostra,aux)
    }
  }
  while(sum(amostra$censura==1)!=n-ncensura){
    t1<-rGTDL(1,lambda,alpha,gamma)
    t2<-rGTDL(1,lambda,alpha,gamma)
    if(min(t1,t2)==t1){
      aux<-c(t1,1)
      amostra<-rbind(amostra,aux)
    }
  }




  return(amostra[-1,])
}




