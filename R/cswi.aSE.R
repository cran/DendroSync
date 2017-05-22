#' Within-group synchrony for homoscedastic compound symmetry model
#'
#' @description The function calculates the within-group synchrony (a^) and standard error (SE) for homoscedastic compound symmetry models (mCS).
#'
#' @usage cswi.aSE(model)
#' 
#' @param model a class "\code{lme}" compound symmetry model (\code{mCS}) produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{TRUE}.
#'
#' @details The function calculates within-group synchrony for homoscedastic compound symmetry model (mCS).
#' 
#' @return 
#' The function returns a \code{matrix} containing within-group synchrony and SE for each combinations of \code{varGroup} levels. This function is used internally in \code{\link{sync}}.  
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#'  
#' @references Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. \emph{Tree Physiology} 34: 819-838.
#' @references Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
#'  
#' @seealso \code{\link{sync}} for a clear description of synchrony evaluation.
#' 
#' @examples ## Calculate within-group homoscedastic synchrony and SE
#'  # for compound symmetry homocedastic model of conifersIP data:
#'  data(conifersIP)
#'  
#'  #Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN)
#'  # using taxonomic grouping criteria (ie. Species)
#'  ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conifersIP, homoscedastic = TRUE)
#'  summary(ModHm)
#'    
#'  #Obtain the compound symmetry model within-group synchrony and SE
#'  # for each varGroup stratum.
#'  cswi.aSE(ModHm$mCS)#compound symmetry
#'    
#' @export cswi.aSE  
#'
cswi.aSE <- function(model){
  if(class(model)[1] != "lme") {stop("a unstructure model (mCS) of nlme class is needed")
  }
  if (model$call[4] == "(~1 | vrTi/vrGr)()"){
    t0 <- as.numeric(VarCorr(model)[4,1])
    t1 <- as.numeric(VarCorr(model)[2,1])
    sig2 <- as.numeric(VarCorr(model)[5,1])
    a_vrGr <- (t0+t1)/(t0+t1+sig2)
    a_Group <- round(a_vrGr, 3)
    Ny <- (max(model$data$vrTi)-min(model$data$vrTi))-1
    SE_vrGr <- sqrt(1.96*a_vrGr^2*(1-a_vrGr)^2/Ny)
    SE_Group <- round(SE_vrGr, 3)
    }
    
    else {stop("a compound compound symmetry structure model (mCS) of nlme class is needed")}
    
  res <- cbind(a_Group, SE_Group)
  return(res)
}
    