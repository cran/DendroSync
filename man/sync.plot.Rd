% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sync.plot.R
\name{sync.plot}
\alias{sync.plot}
\title{Plot within- and between-group synchrony}
\usage{
sync.plot (syncList)
}
\arguments{
\item{syncList}{a \code{list} of the type as produced by \code{\link{sync}}.}
}
\value{
Dotplot
}
\description{
The function creates dot plots of within- and between-group synchrony as produced by \code{\link{sync}} from a selected model produced by \code{\link{dendro.varcov}}.
             Note that broad evaluation model (mBE) can not be plotted since it produces only one value per model.
}
\details{
The function makes a dot plots for within- and between-group synchrony for a user defined \code{varGroup} and \code{varTime} period in \code{\link{dendro.varcov}}.
}
\examples{
## Plot homoscedastic narrow evaluation (mNE) and unstructured model (mUN)
 # synchronies for conifersIP data:
 data(conifersIP)
     
 ##Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN)
 # using geographic grouping criteria (ie. Region)
 ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Region", 
                        data = conifersIP, homoscedastic = TRUE)
 
 sync.plot(sync(ModHm, modname = "mNE"))
 sync.plot(sync(ModHm, modname = "mUN"))


}
\author{
Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
}
