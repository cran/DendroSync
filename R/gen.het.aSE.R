#' Within-group synchrony for heteroscedastic mixed models
#'
#' @description The function calculates for each \code{varGroup} stratum the synchrony (a^) and standard error (SE), but only for heteroscedastic unstructured and narrow evaluation mixed models (mHeNE, mHeUN).
#'
#' @usage gen.het.aSE (model)
#' 
#' @param model a class "\code{lme}" model produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{FALSE}.
#'
#' @details The function calculates the within-group synchrony values for each \code{varGroup} stratum based on the metodology described in Shestakova et al. (2014) and Shestakova et al. (2016). Note that this function is designed to work only in heteroscedastic narrow evaluation and unstructured mixed models (ie. mHeNE and mHeUN).
#' 
#' @return 
#' The function returns a \code{matrix} containing  within-group synchrony and SE for each level of \code{varGroup}. This function is used internally in \code{\link{sync}}.  
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @references Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. \emph{Tree Physiology} 34: 819-838.
#' @references Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
#' 
#' @seealso \code{\link{sync}} for a clear description of synchrony evaluation.
#' 
#' @examples ## Calculate within-group heterosdecastic synchrony and SE for conifersIP data:
#'  data(conifersIP)
#'  
#'  #Fit the heteroscedastic set of varcov models (mBE, mHeNE, mHeCS, mHeUN)
#'  # using taxonomic grouping criteria (i.e. Species)
#'  ModHt <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conifersIP, homoscedastic = FALSE)
#'  summary(ModHt)
#'    
#'  #Obtain the heteroscedastic within-group synchrony and SE for each varGroup stratum.
#'  gen.het.aSE(ModHt$mHeNE)#Narrow evaluation model
#'  gen.het.aSE(ModHt$mHeUN)#Unstructured model
#'  
#' @import stats
#'    
#' @export gen.het.aSE 
#'
gen.het.aSE <- function(model){
  if(class(model)[1] != "lme") {stop("a model of nlme class is needed")
    }
  if(model$call[5] == "varIdent(form = ~1 | vrGr)()" && model$call[4] != "list(vrTi = pdCompSymm(~vrGr - 1))()"){
    Sigma <- as.numeric(VarCorr(model)[, 1])
    Sigma <- as.list(na.omit(Sigma[-length(Sigma)]))
    SpHet <- het.var(model)} 
  else {stop("a heteroscedastic mHeNE or mHeUN model of nlme class is needed")}
  
  a_Spp <- sapply(Sigma, function(Sigma1) {
    Z <- (Sigma1)/(Sigma1+SpHet) })
  
  ASpp <- as.list(diag(a_Spp))
  SE_Spp <-  sapply(ASpp, function(ASpp1) {
    Ny <- (max(model$data$vrTi)-min(model$data$vrTi))-1
    Z <- sqrt(1.96*(ASpp1)^2*(1-ASpp1)^2/Ny) }) 
  
  a_Group <- round(as.numeric(ASpp), 3)
  SE_Group <- round(as.numeric(SE_Spp), 3)
  
  res <- cbind(a_Group, SE_Group)
  return(res)
}






