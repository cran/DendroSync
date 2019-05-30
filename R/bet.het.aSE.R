#' Between-group synchrony for a heteroscedastic unstructured model
#'
#' @description The function calculates the between-group synchrony (a^) and standard error (SE) for a heteroscedastic unstructured model (mHeUN).
#'
#' @usage bet.het.aSE(model)
#' 
#' @param model a class "\code{lme}" unstructured model (\code{mHeUN}) produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{FALSE}.
#'
#' @details The function calculates between-group synchrony for a heteroscedastic unstructured model (mHeUN).
#' 
#' @return 
#' The function returns a \code{matrix} containing between-group synchrony and SE for each combination of \code{varGroup} levels. This function is used internally in \code{\link{sync}}.  
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @references Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. \emph{Tree Physiology} 34: 819-838.
#' @references Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
#'   
#' @seealso \code{\link{sync}} for a clear description of synchrony evaluation.
#'   
#' @examples ## Calculate between-group synchrony and SE 
#'  ##for a heteroscedastic unstructured model for conifersIP data:
#'  data(conifersIP)
#'  
#'  #Fit the heteroscedastic set of varcov models (mBE, mHeNE, mHeCS, mHeUN)
#'  # using taxonomic grouping criteria (ie. Species)
#'  ModHt <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conifersIP, homoscedastic = FALSE)
#'  
#'  #between-group synchrony and SE for each varGroup stratum combination
#'  # in heteroscedastic unstructured models.
#'  bet.het.aSE(ModHt$mHeUN)#Unstructured model
#'  
#' @export bet.het.aSE  
#'
bet.het.aSE  <- function(model){
  if(class(model)[1] != "lme") {stop("a unstructure heteroscedastic model (mHeUN) of nlme class is needed")
    }
  if (model$call[5] == "varIdent(form = ~1 | vrGr)()" && model$call[4] == "list(vrTi = pdSymm(~vrGr - 1))()"){
    Sigma <- as.numeric(VarCorr(model)[, 1])
    nc <- dim(getVarCov(model))[2]
    nsa <- rev(c(1:nc-1))
    Sigma1 <- rep((Sigma[1:nc]), nsa)
    Siga <- c()
    for (i in 1:nc){
      Siga[[i]] <- as.list((Sigma[i:nc]))
    }
    Sigma3 <- unlist(Siga[2:nc])
    Sovma <- getVarCov(model)
    Sovma1 <- Sovma[lower.tri(Sovma) == T]
    SpHet <- het.var(model)
    SpHet1 <- rep((SpHet[1:nc]), nsa)
    for (i in 1:nc){
      Siga[[i]] <- as.list((SpHet[i:nc]))
    }
    SpHet2 <- unlist(Siga[2:nc])
  } 
  else {stop("a unstructure heteroscedastic model (mHeUN) of nlme class is needed")}
  
  Ba_Spp <- sapply(Sovma1, function(Sovma2) {
    Z <- (Sovma2)/sqrt((Sigma1+SpHet1)*(Sigma3+SpHet2))
    } )
  
  Ba_Spp <- diag(as.matrix(Ba_Spp))
  
  BSpp <- as.list(Ba_Spp)
  BetSE_Spp <-  sapply(BSpp, function(BSpp1) {
    Ny <- (max(model$data$vrTi)-min(model$data$vrTi))-1
    Z <- sqrt(1.96*(BSpp1)^2*(1-BSpp1)^2/Ny)
    } ) 
  
  a_betw_Grp <- round(as.numeric(Ba_Spp), 3)
  SE_betw_Grp <- round(as.numeric(BetSE_Spp), 3)
  
  res <- cbind(a_betw_Grp, SE_betw_Grp)
  return(res)
}

