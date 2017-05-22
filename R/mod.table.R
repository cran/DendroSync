#' Function to calculate goodness-of-fit statistics for variance-covariance models
#'
#' @description The function provides a table to compare fitted variance-covariance (VCOV) mixed models by AIC, AICc, BIC and LogLik. The restricted log-likelihood (LogLik) statistics for different models can be compared by Chi-square test, while Akaike information criterion (AIC), corrected Akaike information criterion (AICc) and Bayesian information criterion (BIC) are in the smaller-is-better form.
#' 
#' @usage mod.table(modelList)
#' 
#' @param modelList a \code{list} of variance-covariance (VCOV) mixed models of the type as produced by \code{\link{dendro.varcov}}.
#'
#' @details The function returns a table to compare the fitted variance-covariance (VCOV) mixed models to the same data based on information criteria. The smaller AIC, AICc or BIC, the better fit. Also, LogLik value is included. 
#'          AICc is calculated according to the formula AIC + 2*npar*(nobs/(nobs-npar-2)), where npar represents the number of parameters and nobs the number of observations in the fitted model. 
#'
#' @return The function returns a \code{data.frame} with rows corresponding to the objects and columns containing the following components:
#'    \item{n}{ the number of observations used in the model fit.}
#'    \item{df}{ the number of parameters in the fitted model.}
#'    \item{AIC}{ Akaike's Information Criterion of the fitted model.}
#'    \item{AICc}{ corrected Akaike's Information Criterion of the fitted model.}
#'    \item{BIC}{ Bayesian Information Criterion of the fitted model.}
#'    \item{LogLik}{log-likelihood of the fitted model}
#'  
#'  
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @references Hurvich, C.M. & Tsai, C.L. (1989). Regression and time series model selection in small samples. \emph{Biometrika} 76: 297-307.
#' 
#' @seealso \code{\link{AIC}}, \code{\link{BIC}}, \code{\link{logLik}}
#' 
#' @examples ## Compare homoscedastic variance-covariance models on Iberian Peninsula
#'  # conifer ring chronologies using taxonomic grouping criteria (i.e. Species).
#'  data(conifersIP)
#'  ModHmSp <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                           data = conifersIP, homoscedastic = TRUE)
#'  
#'  mod.table(ModHmSp)# a data.frame containing information criterion values
#' 
#'  ## Compare homoscedastic variance-covariance models on Iberian Peninsula conifers
#'  # ring chronologies using geographic criteria (ie. Region).
#'  ModHmGoe <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Region", 
#'                            data = conifersIP, homoscedastic = TRUE)
#'  
#'  mod.table(ModHmGoe)
#' 
#' @import stats
#' 
#' @export mod.table
#' 
#'  
mod.table <- function(modelList){
  
    if (all(class(modelList)!= "list")){modelList = list(modelList)}
    if (all(class(modelList$mBE) != "lme")){stop("'modelList' is no list output of function dendro.varcov")}
    output <- do.call(rbind, lapply(modelList, function(model) {
    output = data.frame(
      n = nobs(model)[1],
      df = attr(logLik(model), "df"),
      AIC = AIC(model),
      AICc = AIC(model) + (2 * (attr(logLik(model), "df")) * (attr(logLik(model), "df"))) / (nobs(model) - attr(logLik(model), "df") - 2),
      BIC = BIC(model),
      LogLik = (-2*logLik(model))
      )          
      }
      ))
      return(output)
  }
  
  
