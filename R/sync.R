#' Calculate within- and between-group synchrony
#' 
#' @description The function calculates spatial synchrony from a \code{list} of fitted mixed models with variance-covariance structures of the type as produced by \code{\link{dendro.varcov}}. Within- and between- \code{varGroup} level synchrony are calculated, quantifying the degree to which the values of N chronologies contain a common temporal signal. Different models allow for the estimation of intraclass correlations either at the intragroup or intergroup level. The underlying idea is to split the mean correlation estimated between all possible pairs of chronologies drawn from the whole dataset into: (i) a mean correlation between pairs of chronologies for every group; and (ii) a mean correlation between pairs of chronologies for pairs of groups.
#'  
#' @usage sync (modelList, modname = c("mBE", "mNE", "mCS", "mUN", "mHeNE",
#'                                     "mHeCS", "mHeUN"), trend.mBE = FALSE)
#' 
#' @param modelList a \code{list} of the type as produced by \code{\link{dendro.varcov}}.
#' @param modname a \code{character} string of \code{"mBE"}, \code{"mNE"}, \code{"mCS"}, \code{"mUN"}, \code{"mHeNE"}, \code{"mHeCS"} or \code{"mHeUN"}, specifying the variance-covariance structures selected for syncrony evaluation.
#' @param trend.mBE a \code{logical} specifying if a broad evaluation model (mBE) output for each grouping level is reported. This is a special mBE output to plot synchrony trends with \code{\link{sync.trend.plot}}. Default \code{FALSE}.
#' 
#' @details The function calculates the within- and between-group synchrony. For the more general (unstructured) model, the correlation of pairs of chronologies \code{i} and \code{i*} belonging to group \code{r} is:
#'          \deqn{rho(Wi,Wi*) = cov(Wi,Wi*)/sqrt(Var(Wi)*Var(Wi*)) = sigma^2yr/sigma^2yr+sigma^2e} 
#'          Where \code{Wi} is tree-ring width of \code{i}th chronology, \code{sigma^2yr} is a covariance between observations \code{Wi} and \code{Wi*} belonging to a group \code{r}, \code{sigma^2e} is a random deviation within the \code{r}th group.
#'          Conversely, the correlation of pairs of chronologies \code{i} and \code{i*} belonging to groups \code{r} and \code{r*} is:
#'          \deqn{rho(Wi,Wi*) = cov(Wi,Wi*)/sqrt(Var(Wi)*Var(Wi*)) =}
#'          \deqn{sigma^2yr*/sqrt((sigma^2yr+sigma^2e)+(sigma^2yr*+sigma^2e))} 
#'          Note that if no \code{modname} is provided a warning message appears indicating that synchrony will be only calculated for the first \code{modname} vector element, i.e. broad evaluation model (mBE). 
#' 
#' @return 
#' The function returns a \code{list} containing the following components:
#' 
#' \itemize{\item{for within-group synchrony:}}
#'  \item{Modname}{a column indicating the variance-covariance mixed models fit type:}
#'    \itemize{\item{mBE}{: null (or broad evaluation) structure.}}
#'    \itemize{\item{mNE}{: homoscedastic variant of banded main diagonal (or narrow evaluation) structure.}}
#'    \itemize{\item{mCS}{: homoscedastic variant of compound symmetry structure.}}
#'    \itemize{\item{mUN}{: homoscedastic variant of unstructured (or full) structure.}}
#'    \itemize{\item{mHeNE}{: heteroscedastic variant of banded main diagonal (or narrow evaluation) structure.}}
#'    \itemize{\item{mHeCS}{: heteroscedastic variant of compound symmetry structure.}}
#'    \itemize{\item{mHeUN}{: heteroscedastic variant of unstructured (or full) structure.}}
#'  \item{a_Group}{a column representing the within-group synchrony.}
#'  \item{SE_Group}{standard error of each observation.}
#'  
#' \itemize{\item{for between-group synchrony:}}
#' \item{Modname}{a column indicating the model fit type. See previous desription.}
#' \item{GroupName}{a column indicating between-group \code{varGroup} pairwise combinations r and r*.}
#' \item{a_betw_Grp}{a column indicating between-group \code{varGroup} synchrony.} 
#' \item{SE_betw_Grp}{standard error of each observation.}
#'
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @references Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. \emph{Tree Physiology} 34: 819-838.
#' @references Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
#' 
#' @seealso \code{\link{dendro.varcov}} for models details.
#' 
#' @examples ## Calculate synchrony for null.model (broad evaluation, mBE) and homoscedastic variant
#'  # of unstructured model (or full, mUN) for conifersIP data, 
#'  # and heteroscedastic variant for 1970-1999 period.
#'  data(conifersIP)
#'  
#'  ##Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN) 
#'  #using taxonomic grouping criteria (i.e. Species)
#'  ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                         data = conifersIP, homoscedastic = TRUE)
#'  
#'  summary(ModHm)# Class and length of list elements
#'  
#'  #Synchrony for mBE and mUN models
#'  sync(ModHm, modname = "mBE")
#'  sync(ModHm, modname = "mUN")
#'  
#'  ##Chop the data from 1970 to 1999.
#'  conif.30 <- conifersIP[conifersIP$Year>1969 & conifersIP$Year<2000,]
#'  summary(conif.30$Year)
#'  
#'  #Fit the heteroscedastic set of variance covariance mixed models (mBE, mHeNE, mHeCS, mHeUN)
#'  # using taxonomic grouping criteria (ie. Species)
#'  ModHt30 <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                           data = conif.30, homoscedastic = FALSE)
#'  sync(ModHt30, modname = "mBE")
#'  sync(ModHt30, modname = "mHeUN")
#'  
#' @import utils    
#' 
#' @export sync
#'
#' 

sync <- function(modelList, modname = c("mBE", "mNE", "mCS", "mUN", "mHeNE", "mHeCS", "mHeUN"), trend.mBE = FALSE){
   
  if (all(class(modelList) != "list")){modelList = list(modelList)} 
  if (all(class(modelList$mBE) != "lme")){stop("'modelList' is no list output of function dendro.varcov")}
  
  if(length(modname) > 1) {
    modname<-modname[1]
    warning("Synchrony will be evaluated only for the first element of modname vector")}
  
  if(any(modname == c("mBE", "mNE", "mCS", "mUN", "mHeNE", "mHeCS", "mHeUN"))) {}
  else{stop("One among the fitted variance-covariance structures must be selected: mBE, mNE, mCS, mUN, mHeNE, mHeCS, mHeUN")}
  
  if(any(modname == names(modelList))) {}
  else{stop("One among the fitted variance-covariance structures must be selected")}
  
  mod_mBE <- modelList$mBE
  GroupName <- levels(mod_mBE$data$vrGr)
  res <- data.frame(cbind(GroupName))
  mh <- (combn(levels(mod_mBE$data$vrGr), 2))
  GroupName2 <- paste(mh[1,], mh[2,], sep = "/")
  res2 <- data.frame(cbind(GroupName2))
  names(res2) <- "GroupName"
  Modname <- modname
  
  if(Modname == "mBE") {
      mod_mBE <- modelList$mBE
      ASE <- gen.aSE(mod_mBE)
      
    if(trend.mBE){
      res1_mBE <- cbind(Modname, res, ASE)
      a_betw_Grp = c(0)
      SE_betw_Grp = c(0)
      res2_mBE <- cbind(Modname, res2, a_betw_Grp, SE_betw_Grp)
                           
      output_mBE <- list(res1_mBE, res2_mBE)
      names(output_mBE)[[1]] <- "Within_Group_Synchrony"
      names(output_mBE)[[2]] <- "Between_Group_Synchrony"    
      class(output_mBE) <- "sync"
    }
      
    else{
      output_mBE <- data.frame(Modname)
      output_mBE$a_all <- ASE[,1]
      output_mBE$SE <- ASE[,2]
      output_mBE <- list(output_mBE)
      names(output_mBE)[[1]] <- "Within_Group_Synchrony"
      class(output_mBE) <- "sync"
      }
      return(output_mBE)
  }
  else{ 
    
  if(Modname == "mCS") {
      mod_mCS <- modelList$mCS
      ASE <- cswi.aSE(mod_mCS)
      res1_mCS <- cbind(Modname, res, ASE)
      
      BSE <- csbet.aSE(mod_mCS)
      res2_mCS <- cbind(Modname, res2, BSE)
      output <- list(res1_mCS, res2_mCS)
    }
  
  if(Modname == "mUN") {
      mod_mUN <- modelList$mUN
      ASE <- gen.aSE(mod_mUN)
      res1_mUN <- cbind(Modname, res, ASE)
        
      BSE <- bet.aSE(mod_mUN)
      res2_mUN <- cbind(Modname, res2, BSE)
      output <- list(res1_mUN, res2_mUN)
    }

  if(Modname == "mNE") {
      mod_mNE <- modelList$mNE
      ASE <- gen.aSE(mod_mNE)
      res1_mNE <- cbind(Modname, res, ASE)
    
      a_betw_Grp <- c(0)
      SE_betw_Grp <- c(0)
      res2_mNE <- cbind(Modname, res2, a_betw_Grp, SE_betw_Grp)
      output <- list(res1_mNE, res2_mNE)
    }

  if(Modname == "mHeCS") {
      mod_mCS_H<-modelList$mHeCS
      ASE <- cswi.het.aSE(mod_mCS_H)
      res1_mCS_H <- cbind(Modname, res, ASE)
      
      BSE <- csbet.het.aSE(mod_mCS_H)
      res2_mCS_H <- cbind(Modname, res2, BSE)
      output <- list(res1_mCS_H, res2_mCS_H)
    }
  
  if(Modname == "mHeUN") {
      mod_mUN_H <- modelList$mHeUN
      ASE <- gen.het.aSE(mod_mUN_H)
      res1_mUN_H <- cbind(Modname, res, ASE)
      
      BSE <- bet.het.aSE(mod_mUN_H)
      res2_mUN_H <- cbind(Modname, res2, BSE)
      output <- list(res1_mUN_H, res2_mUN_H)
    }  
  
  if(Modname == "mHeNE") {
      mod_mNE_H <- modelList$mHeNE
      ASE <- gen.het.aSE(mod_mNE_H)
      res1_mNE_H <- cbind(Modname, res, ASE)
      
      a_betw_Grp <- c(0)
      SE_betw_Grp <- c(0)
      res2_mNE_H <- cbind(Modname, res2, a_betw_Grp, SE_betw_Grp)
      output <- list(res1_mNE_H, res2_mNE_H)
    }
    
    names(output)[[1]] <- "Within_Group_Synchrony"
    names(output)[[2]] <- "Between_Group_Synchrony"
    class(output) <- "sync"      
    return(output)
  }
}


    
  
  
  
  