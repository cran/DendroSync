#' Between-group synchrony for a homoscedastic unstructured model
#'
#' @description The function calculates the between-group synchrony (a^) and standard error (SE) for a homoscedastic unstructured model or full model (mUN).
#'
#' @usage bet.aSE(model)
#' 
#' @param model a class "\code{lme}" unstructured model (\code{mUN}) produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{TRUE}.
#'
#' @details The function calculates between-group synchrony for a homoscedastic unstructured model (mUN).
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
#'  # for an unstructured homocedastic model of conifersIP data:
#'  data(conifersIP)
#'  
#'  #Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN)
#'  # using taxonomic grouping criteria (ie. Species)
#'  ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conifersIP, homoscedastic = TRUE)
#'  summary(ModHm)
#'    
#'  #Obtain the unstructured model between-group synchrony and SE
#'  # for each varGroup stratum.
#'  bet.aSE(ModHm$mUN)#Unstructure
#'    
#' @export bet.aSE  
#'
bet.aSE <- function(model){
    if(class(model)[1] != "lme") {stop("a unstructure model (mUN) of nlme class is needed")
        }
    if (model$call[4] == "list(vrTi = pdSymm(~vrGr - 1))()"){
        Sigma <- as.numeric(VarCorr(model)[, 1])
        nc <- dim(getVarCov(model))[2]
        nsa <- rev(c(1:nc-1))
        Sigma1 <- rep((Sigma[1:nc]), nsa)
        Siga <- c()
        for (i in 1:nc){
          Siga[[i]] <- as.list((Sigma[i:nc]))}
        
        Sigma3 <- unlist(Siga[2:nc])
        Sovma <- getVarCov(model)
        Sovma1 <- Sovma[lower.tri(Sovma) == T]
        }
  
    else {stop("a unstructure model (mUN) of nlme class is needed")}
  
    Ba_Spp <- sapply(Sovma1, function(Sovma2) {
      Z <- (Sovma2)/sqrt((Sigma1+summary(model)$sigma^2)*(Sigma3+summary(model)$sigma^2))
      })
      Ba_Spp <- diag(Ba_Spp) 
      BSpp <- as.list(Ba_Spp)

    BetSE_Spp <- sapply(BSpp, function(BSpp1) {
      Ny <- (max(model$data$vrTi)-min(model$data$vrTi))-1
      Z <- sqrt(1.96*(BSpp1)^2*(1-BSpp1)^2/Ny)
      }) 
  
    a_betw_Grp <- round(as.numeric(Ba_Spp), 3)
    SE_betw_Grp <- round(as.numeric(BetSE_Spp), 3)
    res <- cbind(a_betw_Grp, SE_betw_Grp)
    
    return(res)
  }

