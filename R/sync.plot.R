#' Plot within- and between-group synchrony
#' 
#' @description The function creates dot plots of within- and between-group synchrony as produced by \code{\link{sync}} from a selected model produced by \code{\link{dendro.varcov}}.
#'              Note that broad evaluation model (mBE) can not be plotted since it produces only one value per model. 
#' 
#' @usage sync.plot (syncList)
#' 
#' @param syncList a \code{list} of the type as produced by \code{\link{sync}}.
#' 
#' @details The function makes a dot plots for within- and between-group synchrony for a user defined \code{varGroup} and \code{varTime} period in \code{\link{dendro.varcov}}. 
#'
#' @return 
#' Dotplot 
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @examples ## Plot homoscedastic narrow evaluation (mNE) and unstructured model (mUN)
#'  # synchronies for conifersIP data:
#'  data(conifersIP)
#'      
#'  ##Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN)
#'  # using geographic grouping criteria (ie. Region)
#'  ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Region", 
#'                         data = conifersIP, homoscedastic = TRUE)
#'  
#'  sync.plot(sync(ModHm, modname = "mNE"))
#'  sync.plot(sync(ModHm, modname = "mUN"))
#' 
#' 
#' @import ggplot2 
#' @import gridExtra
#' 
#'  
#' @export sync.plot
#' 
#' 
sync.plot <- function(syncList){
  
  stopifnot(is.list(syncList))
  if (class(syncList) != "sync") {
    stop("'syncList' is no a list output of function sync")
  }
  
  if(is.data.frame(syncList[1]) != FALSE) {
    stop("'syncList' is no a list output of function sync")
  }
  
  if(is.data.frame(syncList[2]) != FALSE) {
    stop("'syncList' is no a list output of function sync")
  }
  
  pd <- position_dodge(.2)
  aza1 <- do.call(rbind, lapply(syncList[1], data.frame, stringsAsFactors = FALSE))
  aza2 <- do.call(rbind, lapply(syncList[2], data.frame, stringsAsFactors = FALSE))

  if(dim(aza1)[1] == 1){stop("Broad evaluation plot has not sense (mBE)")}
  
  aexp <- expression(paste(bold("Within-group "),bolditalic(hat(a)["c"])))
  mdn <- aza1[2,1]
  bexp <- paste(mdn)
  
  aa1 <- aza1[[3]]-aza1[[4]]
  aa2 <- aza1[[3]]+aza1[[4]]
  am2 <- max(aa2)
  
  p1 <- ggplot(aza1, aes(x=aza1$GroupName,y=aza1$a_Group))+
    geom_errorbar(aes(ymin = aa1, ymax = aa2), width = 0.2, size = 0.7, position = pd, col = 4) +
    geom_point(shape = 16, size = 4, position = pd, col = 4) +
    labs(x = "Grouping variable", y = aexp)+
    expand_limits(y = 0)+
    scale_y_continuous()+
    ggtitle(bexp)+
    theme_bw()+
    theme(axis.title.y = element_text(vjust = 1.8),
          axis.title.x = element_text(vjust = -0.5),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=11),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")
          )
  
  aexp1 <- expression(paste(bold("Between-group "),bolditalic(hat(a)["c"])))
  ab1 <- aza2[[3]]-aza2[[4]]
  ab2 <- aza2[[3]]+aza2[[4]]
  
  cexp <- paste(mdn)
  
  p2 <- ggplot(aza2, aes(x = aza2$GroupName, y = aza2$a_betw_Grp))+
    geom_errorbar(aes(ymin = ab1, ymax = ab2), width = 0.2, size = 0.7, position = pd, col = "royalblue") +
    geom_point(shape = 16, size = 4, position = pd, col = "royalblue") +
    labs(x = "Grouping variable",y = aexp1)+
    expand_limits(y = c(0, am2))+
    scale_y_continuous()+
    ggtitle(cexp)+
    theme_bw()+
    theme(axis.title.y = element_text(vjust = 1.8),
          axis.title.x = element_text(vjust = -0.5),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=11),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")
          )
  
  grid.arrange(p1, p2, ncol = 2)
  }



