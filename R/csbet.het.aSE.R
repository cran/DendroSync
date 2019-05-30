#' Between-group synchrony for heteroscedastic compound symmetry model
#'
#' @description The function calculates the between-group synchrony (a^) and standard error (SE) for heteroscedastic compound symmetry model (mHeCS).
#'
#' @usage csbet.het.aSE(model)
#' 
#' @param model a class "\code{lme}" compound symmetry model (\code{mHeCS}) produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{FALSE}.
#'
#' @details The function calculates between-group synchrony for heteroscedastic compound symmetry models (mHeCS).
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
#'  ##for heteroscedastic compound symmetry model for conifersIP data:
#'  data(conifersIP)
#'  
#'  #Fit the heteroscedastic set of varcov models (mBE, mHeNE, mHeCS, mHeUN)
#'  # using taxonomic grouping criteria (ie. Species)
#'  ModHt <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conifersIP, homoscedastic = FALSE)
#'  
#'  #between-group synchrony and SE for each varGroup stratum combination
#'  # in heteroscedastic compound symmetry models.
#'  csbet.het.aSE(ModHt$mHeCS)
#'  
#' @export csbet.het.aSE  
#'
csbet.het.aSE  <- function(model){
  if(class(model)[1] != "lme") {stop("a heteroscedastic compound symmetry model of nlme class is needed (mHeCS)")
    }
  if(model$call[5] == "varIdent(form = ~1 | vrGr)()" && model$call[4] == "list(vrTi = pdCompSymm(~vrGr - 1))()"){
    SpHet <- het.var(model)
    t0 <- as.numeric(getVarCov(model)[1,1])
    t1 <- as.numeric(getVarCov(model)[2,1])
    t2 <- t0-t1
    nc <- length (SpHet)
    nsa <- rev(c(1:nc-1))
    SpHet1 <- rep(SpHet, nsa)
    Siga <- c()
    for (i in 1:nc){
      Siga[[i]] <- as.list((SpHet[i:nc]))
    }
    SpHet2 <- unlist(Siga[2:nc])
    Ba_H_vrGr <- sapply(SpHet1, function(SpHet1) {
      Z <- (t1)/sqrt((t1+t2+SpHet1)*(t1+t2+SpHet2))
    } )
    Ba_H_vrGr <- diag(as.matrix(Ba_H_vrGr)) 
    BvrGr <- as.list(Ba_H_vrGr)
    B_H_SE_vrGr <- sapply(BvrGr, function(BvrGr1) {
      Ny <- (max(model$data$vrTi)-min(model$data$vrTi))-1
      Z <- sqrt(1.96*(BvrGr1)^2*(1-BvrGr1)^2/Ny)
    } ) 
  } 
  else {stop("a compound symmetry heteroscedastic model (mHeCS) of nlme class is needed")}
  
  a_betw_Grp <- round(as.numeric(Ba_H_vrGr), 3)
  SE_betw_Grp <- round(as.numeric(B_H_SE_vrGr), 3)
  res <- cbind(a_betw_Grp, SE_betw_Grp)
  return(res)
}

