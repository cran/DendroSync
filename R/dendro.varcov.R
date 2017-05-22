#' Fit variance-covariance mixed models on tree-ring chronologies
#'
#' @description The function calculates variance-covariance (VCOV) mixed models from a \code{data.frame} with tree-ring width index and years for each chronology following the methodology described in Shestakova et al. (2014). 
#'              The mixed models relate tree-ring width (\code{Y}) against specific names of tree-ring width chronologies (\code{A}), using years and grouping variable as random factors to characterize the strength of the common signal across the grouping variable.
#'              First, a linear mixed-effect model with null positive-definite matrix structure or broad evaluation is fitted and the subsequent models are consequently derived from it using the function \code{update}.
#'              When a \code{data.frame} with tree-ring width index has NAs the models are fitted with \code{\link{na.action}} = \code{na.omit}. Simultaneously, \code{\link{complete.cases}} is applied to guarantee that rows have no missing values across the entire \code{data.frame}.
#'
#' @usage dendro.varcov(formula, varTime = "", varGroup = "", data, 
#'                       homoscedastic = TRUE, null.mod = FALSE, all.mod = FALSE)
#' 
#' @param formula a model \code{formula} such as \code{Y ~ A}, where \code{Y} is usually tree-ring width and \code{A} will be a factor variable such as the specific names of tree-ring width chronologies (\code{\link{conifersIP}}) or \code{~1} for a null model.
#' @param varTime a \code{character} specifying the time variable to consider in calculating synchrony estimates. Models with varTime variable with less than 10 different time-points produce unreliable results. 
#' @param varGroup a \code{character} grouping variable. In dendrochronological studies different grouping strategies can be used. We used here two strategies following taxonomic (i.e. species) or geographic (i.e. region) criteria.
#' @param data a \code{data.frame} with tree-ring chronologies, years and grouping variables as columns.
#' @param homoscedastic \code{logical} if \code{TRUE} models do not included an optional \code{\link{varFunc}} object. If \code{FALSE} models will include a one-sided formula describing the within-group heteroscedasticity structure (\code{\link{varIdent}}). 
#' @param null.mod \code{logical} if \code{TRUE} only broad evaluation model will be fitted. Default \code{FALSE}.
#' @param all.mod \code{logical} if \code{TRUE} all homoscedastic and heteroscedastic model types will be fitted. Default \code{FALSE}.
#' 
#' @details The function fits a set of variance-covariance mixed models following Shestakova et al. (2014). A total of 7 different variance-covariance mixed models can be fitted: a null positive-definite matrix structure (mBE), and the homoscedastic and heteroscedastic versions of a diagonal positive-definite matrix structure (mNE, mHeNE), a positive-definite matrix with compound symmetry structure (mCS, mHeCS) and a general positive-definite matrix structure (mUN, mHeUN). Note that if null.mod is \code{TRUE} the function only fits broad evaluation model (mBE), this is set to \code{FALSE} by default. If all.mod is \code{TRUE} the function fits heteroscedastic and homoscesdastic versions of all models. This is set to \code{FALSE} by default, because for large-datasets it may take a long time to converge.
#'
#' @return The function returns a \code{list} containing the following components:
#'  \itemize{\item{for \code{null.mod = TRUE}:}}
#'  \item{mBE}{ an object of class "lme" representing the linear mixed-effects model fit of null positive-definite matrix structure or broad evaluation. See \code{\link{lmeObject}} for the components of the fit.}
#' 
#'  \itemize{\item{for \code{homoscedastic = TRUE}:}}
#'  \item{mNE}{ an object of class "lme" representing the linear mixed-effects model fit of a diagonal positive-definite matrix structure or narrow evaluation. See \code{\link{lmeObject}} for the components of the fit.} 
#'  \item{mCS}{ an object of class "lme" representing the linear mixed-effects model fit of a positive-definite matrix with compound symmetry structure. See \code{\link{lmeObject}} for the components of the fit.}
#'  \item{mUN}{ an object of class "lme" representing the linear mixed-effects model fit of a general positive-definite matrix structure or unstructured. See \code{\link{lmeObject}} for the components of the fit.} 
#' 
#'  \itemize{\item{for \code{homoscedastic = FALSE}:}}
#'  \item{mHeNE}{ an object of class "lme" representing the linear mixed-effects model fit of the heteroscedastic variant of a diagonal positive-definite matrix structure or narrow evaluation. See \code{\link{lmeObject}} for the components of the fit.}
#'  \item{mHeCS}{ an object of class "lme" representing the linear mixed-effects model fit of the heteroscedastic variant of a positive-definite matrix with compound symmetry structure. See \code{\link{lmeObject}} for the components of the fit.} 
#'  \item{mHeUN}{ an object of class "lme" representing the linear mixed-effects model fit of the heteroscedastic variant of a general positive-definite matrix structure or unstructured. See \code{\link{lmeObject}} for the components of the fit.} 
#' 
#'  \itemize{\item{for \code{all.mod = TRUE}:}}
#'  \item{}{ The function returns the homoscedastic and heteroscedastic versions of all fitted models.} 
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @references Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. \emph{Tree Physiology} 34: 819-838.
#' @references Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
#' 
#' @seealso \code{\link{lmeObject}}, \code{\link{na.action}}, \code{\link{complete.cases}}
#' 
#' @examples ## Calculate variance-covariance models on Iberian Peninsula conifers
#'  # chronologies using two different grouping strategies.
#'  # Tree-ring width chronologies are grouped according to taxonomic (i.e. Species)
#'  # or geographic (i.e. Region) criteria.
#'  #User-defined homoscedastic or heteroscedastic variances can be fitted.
#'  data(conifersIP)
#'  
#'  #Chop the data from 1960 to 1989.
#'  conif.30 <- conifersIP[conifersIP$Year>1959 & conifersIP$Year<1990,]
#'  summary(conif.30$Year)
#'  
#'  ##Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN)
#'  # using taxonomic grouping criteria (ie. Species)
#'  ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conif.30, homoscedastic = TRUE)
#'  
#'  summary(ModHm)# Class and length of list elements
#'  ModHm
#'  ModHm[2]#mNE fitted model results
#'  
#'  ##Fit the heteroscedastic set of varcov models (mBE, mHeNE, mHeCS, mHeUN) 
#'  # using geographic grouping criteria (ie. Region)
#'  ModHt <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Region", 
#'                         data = conif.30, homoscedastic = FALSE)
#'  
#'  summary(ModHt)# Class and length of list elements
#'  ModHt
#'  ModHt[3]#mHeCS fitted model results
#'        
#'                    
#' @import nlme
#' @import stats
#' 
#' @export dendro.varcov
#' 
#'  
dendro.varcov<-function(formula, varTime = "", varGroup = "", data = stop("A data.frame must be provided"),
                      homoscedastic = TRUE, null.mod = FALSE, all.mod = FALSE) {
  
  data$vrTi <- data[,varTime]
  data$vrGr <- data[,varGroup]
  if(data.class(data$vrTi) != "numeric")
    {stop("a numeric time variable must be provided (e.g years)")}
  
  if(length(unique(data$vrTi)) <= 10)
    {warning("for a precise evaluation of syncrony more than 10 different varTime points are needed")}
  
  if(data.class(data$vrGr) != "factor")
    {stop("a factor Grouping variable must be provided (e.g Species names)")}
  
  data <- data[complete.cases(data),]
  
  print("Please wait. I am fitting the models now :)")
  
    mBE <- lme(formula, random = ~1|vrTi, data = data, control = lmeControl(msMaxIter = 100, msVerbose = F), na.action = na.omit)
  if(null.mod){
    output <- list()
    output$mBE <- mBE
    return(invisible(output))
    }
  else{
    if(all.mod){
      mCS <- update(mBE, random = ~1|vrTi/vrGr, correlation = corCompSymm(form = ~1|vrTi/vrGr))
      mUN <- update(mBE, random = list(vrTi = pdSymm(~vrGr-1)))
      mNE <- update(mBE, random = list(vrTi = pdDiag(~vrGr-1)))
      mHeCS <- update(mBE, random = list(vrTi = pdCompSymm(~vrGr-1)), weights = varIdent(form = ~1|vrGr))
      mHeUN <- update(mBE, random = list(vrTi = pdSymm(~vrGr-1)), weights = varIdent(form = ~1|vrGr))
      mHeNE <- update(mBE, random = list(vrTi = pdDiag(~vrGr-1)), weights = varIdent(form = ~1|vrGr))
      output3 <- list()
      output3$mBE <- mBE
      output3$mNE <- mNE
      output3$mCS <- mCS
      output3$mUN <- mUN
      output3$mHeNE <- mHeNE
      output3$mHeCS <- mHeCS
      output3$mHeUN <- mHeUN
      return(invisible(output3)) 
    }
    else{
      if(homoscedastic){
        mCS <- update(mBE, random = ~1|vrTi/vrGr, correlation = corCompSymm(form = ~1|vrTi/vrGr))
        mUN <- update(mBE, random = list(vrTi = pdSymm(~vrGr-1)))
        mNE <- update(mBE, random = list(vrTi = pdDiag(~vrGr-1)))
        output <- list()
        output$mBE <- mBE
        output$mNE <- mNE
        output$mCS <- mCS
        output$mUN <- mUN
        return(invisible(output))
      }  
      else{
        mHeCS <- update(mBE, random= list(vrTi = pdCompSymm(~vrGr-1)), weights = varIdent(form = ~1|vrGr))
        mHeUN <- update(mBE, random= list(vrTi = pdSymm(~vrGr-1)), weights = varIdent(form = ~1|vrGr))
        mHeNE <- update(mBE, random= list(vrTi = pdDiag(~vrGr-1)), weights = varIdent(form = ~1|vrGr))
        output2 <- list()
        output2$mBE <- mBE
        output2$mHeNE <- mHeNE
        output2$mHeCS <- mHeCS
        output2$mHeUN <- mHeUN
        return(invisible(output2)) 
    }
    }
    }
  }      
      







  
  
  
  