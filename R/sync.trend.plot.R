#' Plot temporal trends of synchrony
#' 
#' @description The function creates a line chart showing temporal trends of spatial synchrony from \code{data.frame} of the type as produced by \code{\link{sync.trend}}.
#' 
#' @usage sync.trend.plot (sync.trend.data)
#' 
#' @param sync.trend.data a \code{data.frame} of the type as produced by \code{\link{sync.trend}}.
#' 
#' @details The function makes a line chart showing synchrony trends across years from a data.frame produced by \code{\link{sync.trend}}. Within- or between- group synchrony and SE are indicated for a selected time window. If synchrony is defined using using \code{null.mod = TRUE} (\code{\link{sync.trend}}) only general synchrony is ploted. If synchrony is defined using using \code{null.mod = FALSE} (\code{\link{sync.trend}}) different synchronies for each group variable (\code{varGroup}) are fitted with different colours for each stratum.
#'
#' @return 
#' Line chart 
#' 
#' @author 
#' Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
#' 
#' @examples ## Calculate temporal trends of synchrony for conifersIP data:
#'  data(conifersIP)
#'  
#'  ##Fit the null.model temporal trend (mBE) using taxonomic grouping criteria (i.e. Species)
#'  mBE.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                           data = conifersIP, null.mod = TRUE, window = 30, lag = 5)
#'  
#'  mBE.trend# it returns a data.frame
#'  sync.trend.plot(mBE.trend)# Broad evaluation synchrony linechart
#'
#' \dontrun{ 
#'  ##Fit homoscedastic within-group trends (mBE, mNE, mCS, mUN) 
#'  # using geographic grouping criteria (i.e. Region)
#'  geo.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Region", 
#'                          data = conifersIP, window = 30, lag = 5, 
#'                          null.mod = FALSE, homoscedastic = TRUE)
#'  
#'  geo.trend#a data.frame with varGroup synchrony for each time window.
#'  sync.trend.plot(geo.trend)#Selected heteroscedastic between-group trends by AIC
#'  
#'  ##Fit heteroscedastic betwen-group trends (mBE, mHeNE, mHeCS, mHeUN) 
#'  # using geographic grouping criteria (i.e. Region) and AICc
#'  geo.het.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Region", 
#'                     data = conifersIP, window = 30, lag = 5, null.mod = FALSE, 
#'                     selection.method = c("AICc"), homoscedastic = FALSE, between.group = TRUE)
#'  
#'  geo.het.trend
#'  sync.trend.plot(geo.het.trend)#Selected heteroscedastic between-group trends by AICc
#'  
#'  ##Fit homoscedastic and heteroscedastic within-group trends 
#'  # using taxonomic grouping criteria (i.e. Species) and BIC
#'  geo.tot.trend <- sync.trend(TRW ~ Code, varTime = "Year", varGroup = "Species", 
#'                     data = conifersIP, window = 30, lag = 5, selection.method = c("BIC"),
#'                     all.mod = TRUE)
#'  geo.tot.trend
#'  #Selected homoscedastic and heteroscedastic within-group trends by BIC
#'  sync.trend.plot(geo.tot.trend)
#'  }
#'  
#'
#' @import ggplot2
#'
#' @export sync.trend.plot
#' 
#'  
#' 
sync.trend.plot <- function(sync.trend.data){
  
  stopifnot(is.data.frame(sync.trend.data))
  if(class(sync.trend.data)[1] != "sync.trend") {
    stop("'sync.trend.data' is no data.frame output of function sync.trend")
  }
   
    dat <- sync.trend.data
    aexp <- expression(paste(hat(a)["c"]))
    
    if(dim(dat)[2] == 4){
          ly <- length(unique(row(dat[4])))-1
          wyr <- max(dat[4])-min(dat[4])    
          yrby <- wyr/ly
          q <- ggplot(dat, aes(x = dat[4], y = dat$a_Group))
          q1 <- q + geom_point(size = 3) + geom_line(size = 1) + guides(colour = FALSE, fill = FALSE)+ 
              geom_ribbon(data = dat, aes(ymin = dat$a_Group - dat$SE, ymax = dat$a_Group + dat$SE), alpha = 0.2)+
              xlab("Windows") + ylab(aexp)+
              scale_x_continuous(limits = c(min(dat[4]), max(dat[4])), breaks = seq(min(dat[4]), max(dat[4]), by = yrby))+
              ggtitle("Broad evaluation synchrony (mBE)")+ theme(plot.title = element_text(hjust = 0.5))
        }
          
    else{
        gexp <- expression(paste(hat(a)["c"], " ", "synchrony"))
        ly <- dim(unique(dat[6]))[1]-1
        wyr <- max(dat[6])-min(dat[6])    
        yrby <- wyr/ly
        aa1 <- dat[3]-dat[4]
        aa2 <- dat[3]+dat[4]
        q <- ggplot(dat, aes(x = dat[6], y = dat[3], group = dat$GroupName, colour = factor(dat$GroupName)))
        q1 <- q+geom_point(size = 3)+
            stat_summary(fun.y = mean, geom = "line", aes(group = dat$GroupName))+
            guides(fill = guide_legend(title = "GroupName"), colour = FALSE)+
            scale_x_continuous(limits = c(min(dat[6]), max(dat[6])), breaks = seq(min(dat[6]), max(dat[6]), by = yrby))+
            scale_y_continuous(limits = c(min(aa1), max(aa2)))+
            geom_ribbon(data = dat, aes(ymin = aa1, ymax = aa2, fill = dat$GroupName), alpha = 0.4)+
            xlab("Windows")+ ylab(aexp)+ ggtitle(gexp)+ theme(plot.title = element_text(hjust = 0.5))
            
        }
    
    print(q1)
    
}
  