#'@name MaxGTDL
#'@title Maximum estimation GTDL
#'
#'@description
#'
#'@param
#'@author Jalmar M. F. Carrasco \email{carrascojalmar@gmail.com}
#'@author Luciano S. Santos \email{lucianno0800@gmail.com}

require(GTDL)

MaxGTDL<-function(start,...){
  likeGTDL<-function(param){
    f<-sum(log((lambda*exp(t*alpha+gamma)/(1+exp(t*alpha+gamma)))*((1+exp(t*alpha+gamma))/(1+exp(gamma)))^(-lambda/alpha)))
    return(f)
  }
  aux<-maxLik(likeGTDL,start = start
              ,grad = NULL,hess = NULL)
  return(aux)
}
