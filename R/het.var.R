#' Variances per varGroup stratum for heteroscedastic models
#'
#' @description The function obtains the heteroscedastic variances for each \code{varGroup} elements for a selected model (mHeCS, mHeUN, mHeNE).
#'
#' @usage het.var(model)
#' 
#' @param model a class "\code{lme}" model produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{FALSE}.
#'
#' @details The function extracts the variances for each \code{varGroup} stratum using the within-group heteroscedastic structure of the fitted models (\code{varIdent} constant variance per stratum).
#'  Note that this function only works for heteroscedastic models: mHeCS, mHeNE,mHeUN.
#' 
#' @return 
#' The function returns a \code{numeric} vector containing the variance per each level of \code{varGroup}. They are used internally to calculate synchrony (\code{\link{sync}}).
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @examples ## Calculate within-group heteroscedastic variances for conifersIP data:
#'  data(conifersIP)
#'  
#'  #Fit the heteroscedastic set of models (mBE, mHeCS, mHeNE, mHeUN)
#'  # using taxonomic grouping criteria (i.e. Species)
#'  ModHt <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conifersIP, homoscedastic = FALSE)
#'  
#'  #Obtain the within-group variances for the model of interest
#'  het.var(ModHt$mHeCS)#Heterogeneous variant of compound symmetry model
#'  het.var(ModHt$mHeUN)#Heterogeneous unstructured model
#'  
#' @import stats
#'    
#' @export het.var  
#'        
het.var <- function(model){
  if(class(model)[1] != "lme") {stop("a model of nlme class is needed")
  }
  if (model$call[5] == "varIdent(form = ~1 | vrGr)()"){
    qaq <- coef(model$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE)} 
  else {stop("a nlme class model with iweights is needed (heteroscedastic)")}
    hetero.Spec <- (model$sigma^2)*(qaq^2)
    hetero.Spec <- hetero.Spec[order(names(hetero.Spec))]
  return(hetero.Spec)
  }







