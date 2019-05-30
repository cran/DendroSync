#' Within-group synchrony for homoscedastic models
#'
#' @description The function calculates for each varGroup stratum the within-group synchrony (a^) and standard error (SE). However, it only works for homoscedastic broad evaluation, narrow evaluation and unstructured models (mBE, mNE, mUN).
#'
#' @usage gen.aSE(model)
#' 
#' @param model a class "\code{lme}" model produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{TRUE}.
#'
#' @details The function calculates the within-group synchrony values for each varGroup stratum based on the metodology described in Shestakova et al. (2014) and Shestakova et al. (2016). Note that this function is designed to work only in 3 homoscedastic models (mBE, mNE, mUN).
#' 
#' @return 
#' The function returns a \code{matrix} containing the synchrony and SE for each level of \code{varGroup}. This function is used internally in \code{\link{sync}}.  
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @references Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. \emph{Tree Physiology} 34: 819-838.
#' @references Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
#' 
#' @seealso \code{\link{sync}} for a clear description of synchrony evaluation.
#' 
#' @examples ## Calculate within-group homoscedastic synchrony and SE for conifersIP data:
#'  data(conifersIP)
#'  
#'  #Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN) 
#'  # using taxonomic grouping criteria (ie. Species)
#'  ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                          data = conifersIP, homoscedastic = TRUE)
#'  summary(ModHm)
#'    
#'  #Obtain the within-group synchrony and SE for each varGroup stratum.
#'  gen.aSE(ModHm$mBE)#Broad evaluation
#'  gen.aSE(ModHm$mUN)#Unstructured
#'  
#' @import stats
#'  
#' @export gen.aSE  
#'  
#'
gen.aSE <- function(model){
  if(class(model)[1] != "lme") {stop("a model of nlme class is needed")
  }
  if (model$call[5] != "varIdent(form = ~1 | vrGr)()" && model$call[4] != "(~1 | vrTi/vrGr)()"){
    Sigma <- as.numeric(VarCorr(model)[, 1])
    Sigma <- as.list(na.omit(Sigma[-length(Sigma)]))}
  else {stop("a homoscedastic mBE, mNE or mUN model of nlme class is needed")}
  
  a_Spp <- sapply(Sigma, function(Sigma1) {
    Z <- (Sigma1)/(Sigma1+summary(model)$sigma^2) })
  
  ASpp <- as.list(a_Spp)
  SE_Spp <- sapply(ASpp, function(ASpp1) {
    Ny <- (max(model$data$vrTi)-min(model$data$vrTi))-1
    Z <- sqrt(1.96*(ASpp1)^2*(1-ASpp1)^2/Ny) }) 
  
  a_Group <- round(as.numeric(a_Spp), 3)
  SE_Group <- round(as.numeric(SE_Spp), 3)
  
  res <- cbind(a_Group, SE_Group)
  return(res)
}






