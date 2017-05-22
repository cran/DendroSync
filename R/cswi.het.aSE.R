#' Within-group synchrony for heteroscedastic compound symmetry mixed models
#'
#' @description The function calculates for each \code{varGroup} stratum the synchrony (a^) and standard error (SE) for heteroscedastic compound symmetry models (mHeCS).
#'
#' @usage cswi.het.aSE (model)
#' 
#' @param model a class "\code{lme}" model produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{FALSE}.
#'
#' @details The function calculates the within-group synchrony values for each \code{varGroup} stratum based on the metodology described in Shestakova et al. (2014) and Shestakova et al. (2016). Note that this function is designed to work only in heteroscedastic compound symmetry mixed models (i.e. mHeCS).
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
#'  #Obtain the heteroscedastic compound symmetry 
#'  #within-group synchrony and SE for each varGroup stratum.
#'  cswi.het.aSE(ModHt$mHeCS)
#'  
#'  
#' @import stats
#'    
#' @export cswi.het.aSE
#'
cswi.het.aSE <- function(model){
  if(class(model)[1] != "lme") {stop("a heteroscedastic compound symmetry model of nlme class is needed (mHeCS")
    }
  if(model$call[5] == "varIdent(form = ~1 | vrGr)()" && model$call[4] == "list(vrTi = pdCompSymm(~vrGr - 1))()"){
    SpHet <- het.var(model)
    di1 <- length(SpHet)+1
    t0 <- as.numeric(getVarCov(model)[1,1])
    t1 <- as.numeric(getVarCov(model)[2,1])
    t2 <- t0-t1
    sig2 <- as.numeric(VarCorr(model)[di1,1])
    a_vrGr <- c()
    for (i in 1:length (SpHet)){
      a_vrGr[[i]] <- as.list(((t0)/(t0+SpHet[i])))
      }
    Ny <- (max(model$data$vrTi)-min(model$data$vrTi))-1
    a_Group <- as.vector(round(unlist(a_vrGr), 3))
    SE_vrGr <- sqrt(1.96*a_Group^2*(1-a_Group)^2/Ny)
    SE_Group <- round(SE_vrGr, 3)
  } 
  else {stop("a heteroscedastic mHeCS model of nlme class is needed")}
  
  res <- cbind(a_Group, SE_Group)
  return(res)
}






