#' Calculate temporal trends of synchrony
#' 
#' @description The function calculates temporal trends of spatial synchrony from a \code{data.frame} with tree-ring width chronologies using a moving window as described in Shestakova et al. (2016). This method splits the time variable (\code{varTime}) in 30 years windows plus a 5 years lag, and in each window the within- or between-group level (\code{varGroup}) synchronies are calculated. The function can also be used to find synchrony with similar time series \code{data.frame} from other fields.
#' 
#' @usage sync.trend (formula, varTime="", varGroup="", data,  window = 30, lag = 5, 
#'                     null.mod = TRUE, selection.method = c("AIC", "AICc", "BIC"), 
#'                     all.mod = FALSE, homoscedastic = TRUE, between.group = FALSE)
#' 
#' @param formula a \code{formula} a typical model formula such as \code{Y ~ A}, where \code{Y} is usually tree-ring width and \code{A} may be a grouping factor such as the specific names of tree-ring width chronologies (\code{\link{conifersIP}}).
#' @param varTime a \code{character} specifying the time variable to consider in calculating synchrony estimates. Models with less than 10 different time-points may produce unreliable results. 
#' @param varGroup a \code{character} grouping variable. In dendrochronological studies different grouping strategies can be used. We used here two strategies following taxonomic (i.e. species) or geographic (i.e. region) criteria.
#' @param data a \code{data.frame} with tree-ring chronologies, years and grouping variables as columns.
#' @param window an \code{integer} specifying the window size (i.e. number of years) to be used to calculate synchrony. Must be greater than 20 (>=20). Defaults to 20.
#' @param lag an \code{integer} specifying the lag that the window is moving (i.e. number of vrTirs moving window) to be used to calculate synchrony. Must be greater than 1 (>=1). Defaults to 5.
#' @param selection.method a \code{character} string of \code{"AIC"}, \code{"AICc"} or \code{"BIC"}, specifying the information criterion used for model selection.
#' @param null.mod a \code{logical} specifying if only the null model for general synchrony is fitted (broad evaluation, mBE). Default \code{TRUE}.
#' @param all.mod a \code{logical} specifying if all homoscedastic and heteroscedastic models should be fitted. Default \code{FALSE}.
#' @param homoscedastic a \code{logical} specifying if models should be an optional \code{varFunc} object or one-sided formula describing the within-group heteroscedasticity structure. Default \code{TRUE}
#' @param between.group a \code{logical} specifying if between-group synchrony is displayed instead of whitin-group synchrony. Default \code{FALSE}. 
#' 
#' @details The function fits by default (\code{"null.mod=T"}) the null model for general synchrony (broad evaluation, mBE) for a specified time window size and lag. If \code{"null.mod=F"} the function calculates \code{homoscedastic} or \code{heteroscedastic} versions of variance-covariance (VCOV) mixed models available (mBE, mNE, mCS, mUN, mHeNE, mHeCS, mHeUN; \code{\link{dendro.varcov}}) for each time window size and lag selected. In each window the best model is chosen based on the minimum information criterion selected between "AIC", "AICc" or "BIC".
#' When no \code{selection.method} is defined by default AIC is used. 
#' If \code{"all.mod=T"} the functions fits the homoscedastic and heteroscedastic versions of the 7 models (this is a higly time consuming process). 
#'
#' @return 
#' The function returns a \code{data.frame} containing the following components:
#' 
#' \itemize{\item{for \code{null.mod} \code{TRUE}:}}
#' \item{a_Group}{a column representing the within-group synchrony (mBE).}
#' \item{SE}{standard error of each observation.}
#' \item{Windlag}{a column representing the lag of the window used to split the time variable. A 0 value means that lag is 0, and then the defined time window starts from minimun varTime value.}
#' \item{varTime}{a column representing the \code{varTime} variable.}
#'  
#' \itemize{\item{for \code{null.mod} \code{FALSE}:}}
#' \item{Modname}{a column indicating the best model fit and the information criterion used.}
#' \item{GroupName}{a column indicating levels of the \code{varGroup} for each time-window selected.}
#' \item{a_Group}{a column indicating within-group synchrony for each \code{varGroup} level at time-window selected.} 
#' \item{a_betw_Grp}{a column indicating between-group synchrony for each \code{varGroup} level at time-window selected. Only if \code{between.group} is set to \code{TRUE}.} 
#' \item{SE}{standard error of each observation.}
#' \item{Windlag}{a column representing the lag of the window used to split the time variable. A 0 value means that lag is 0, and then the defined time window starts from minimun varTime value.}
#' \item{varTime}{a column representing the \code{varTime} variable window mean point.}
#'
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @references Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. Tree Physiology 34: 819-838.
#' @references Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
#'  
#' @examples ## Calculate  temporal trends of spatial synchrony for conifersIP data:
#'  data(conifersIP)
#'  
#'  ##Fit the null.model temporal trend (mBE) 
#'  #using taxonomic grouping criteria (i.e. Species)
#'  mBE.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                           data = conifersIP, null.mod = TRUE, window = 30, lag = 5)
#'  
#'  mBE.trend# it returns a data.frame
#'
#' \dontrun{ 
#'  ##Fit homoscedastic within-group trends (mBE, mNE, mCS, mUN) 
#'  # using geographic grouping criteria (i.e. Region)
#'  geo.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Region", 
#'                          data = conifersIP, window = 30, lag = 5, 
#'                          null.mod = FALSE, homoscedastic = TRUE)
#'                          
#'  geo.trend#a data.frame with varGroup syncrony for each time window.
#'  
#'  ##Fit heteroscedastic between-group trends (mBE, mHeNE, mHeCS, mHeUN) 
#'  #using geographic grouping criteria (i.e. Region) and BIC
#'  geo.het.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Region", 
#'                              data = conifersIP, window = 30, lag = 5, null.mod = FALSE, 
#'                              selection.method = c("BIC"), homoscedastic = FALSE, 
#'                              between.group = TRUE)
#'  geo.het.trend
#'  
#'  ##Fit homoscedastic and heterocedastic within-group trends
#'  # using taxonomic grouping criteria (i.e. Species) and BIC
#'  geo.tot.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                              data = conifersIP, window = 30, lag = 5, 
#'                              selection.method = c("BIC"), all.mod = TRUE)
#'  geo.tot.trend
#'  }
#' 
#' @import stats
#'    
#' @export sync.trend
#' 
#' 
#' 
sync.trend <- function(formula, varTime = "", varGroup = "", data = stop("A dataset must be provided"),
                        window = 30, lag = 5, null.mod = TRUE, selection.method = c("AIC", "AICc", "BIC"), all.mod = FALSE, homoscedastic = TRUE, between.group = FALSE){ 
              
         stopifnot(is.numeric(window), length(window) == 1, is.finite(window))
          if(window < 20) {stop("'window' must be > 20")}
         stopifnot(is.numeric(lag), length(lag) == 1, is.finite(lag))
          if(lag < 1) {stop("'lag' must be > 1")}

         data$vrTi <- data[,varTime]
         yrs.len <- (max(data$vrTi)-min(data$vrTi)+1)
         if(yrs.len < window) {stop("'varTime' must be longer than the window length")}
         
         tini <- min(data$vrTi)
         win.t <- (yrs.len-window)/lag
         wt <- floor(win.t)+1
         lwind <- c(1:wt-1)
         yr <-  sapply(lwind, function(S1) {Z <- tini+(window/2)+(lag*S1)} )
         
         mod.val <- list()
         
    if(null.mod){
        for(i in 1:wt) {
          chopval <- data[data$vrTi >= I(tini+(lag*(i-1))) & data$vrTi <= I(tini+window+(lag*(i-1))),]
          modHom <- dendro.varcov(formula, varTime = varTime, varGroup = varGroup, data = chopval, null.mod = T)
          a.mo <- rownames(mod.table(modHom))
          mod.val[[i]] <- sync(modHom, modname = a.mo)[1]
          }
      df <- do.call(rbind, lapply(mod.val, data.frame, stringsAsFactors = FALSE))[2:3]
      df$lagwindow <- lwind
      df$yr <- yr
      names(df) <- c("a_Group", "SE", "Windlag", varTime)
      output <- df
    }
         
    else{
        if(any(selection.method == c("AIC", "AICc", "BIC"))) {
            sel.met <- match.arg(selection.method, c("AIC", "AICc", "BIC"))}

        else{sel.met <- c("AIC")}
        
        if(between.group){
            for(i in 1:wt) {
              chopval <- data[data$vrTi >= I(tini+(lag*(i-1))) & data$vrTi <= I(tini+window+(lag*(i-1))),]
              modHom <- dendro.varcov(formula, varTime = varTime, varGroup = varGroup, data = chopval, all.mod = all.mod, homoscedastic = homoscedastic)
              mo.ta <- mod.table(modHom)
              w <- which(mo.ta[,sel.met] == min(mo.ta[,sel.met]))
              a.mo <- rownames(mo.ta[w,])
              mod.val[[i]] <- sync(modHom, modname = a.mo, trend.mBE = T)[2]
              }
            mv <- do.call(rbind, lapply(mod.val, data.frame, stringsAsFactors = FALSE))
            nl <- sum(1:(nlevels(data[,varGroup])-1))
            mv$lagwindow <- rep(1:wt-1, each = nl)
            mv$yr <- rep(yr, each = nl)
            ModName <- paste("Modname", sel.met, sep = "_")
            names(mv) <- c(ModName, "GroupName", "a_betw_Grp", "SE", "Windlag", varTime)
          output <- mv
          }
          else{
            for(i in 1:wt) {
              chopval <- data[data$vrTi >= I(tini+(lag*(i-1))) & data$vrTi <= I(tini+window+(lag*(i-1))),]
              modHom <- dendro.varcov(formula, varTime = varTime, varGroup = varGroup, data = chopval, all.mod = all.mod, homoscedastic = homoscedastic)
              mo.ta <- mod.table(modHom)
              w <- which(mo.ta[,sel.met] == min(mo.ta[,sel.met]))
              a.mo <- rownames(mo.ta[w,])
              mod.val[[i]] <- sync(modHom, modname = a.mo, trend.mBE = T)[1]
            }
            mv <- do.call(rbind, lapply(mod.val, data.frame, stringsAsFactors = FALSE))
            nl <- nlevels(data[,varGroup])
            mv$lagwindow <- rep(1:wt-1, each = nl)
            mv$yr <- rep(yr, each = nl)
            ModName <- paste("Modname", sel.met, sep = "_")
            names(mv) <- c(ModName, "GroupName", "a_Group", "SE", "Windlag", varTime)
            output <- mv
          }
          }
         class(output) <- c("sync.trend", "data.frame")
         return(invisible(output))
}


