% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen.het.aSE.R
\name{gen.het.aSE}
\alias{gen.het.aSE}
\title{Within-group synchrony for heteroscedastic mixed models}
\usage{
gen.het.aSE (model)
}
\arguments{
\item{model}{a class "\code{lme}" model produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{FALSE}.}
}
\value{
The function returns a \code{matrix} containing  within-group synchrony and SE for each level of \code{varGroup}. This function is used internally in \code{\link{sync}}.
}
\description{
The function calculates for each \code{varGroup} stratum the synchrony (a^) and standard error (SE), but only for heteroscedastic unstructured and narrow evaluation mixed models (mHeNE, mHeUN).
}
\details{
The function calculates the within-group synchrony values for each \code{varGroup} stratum based on the metodology described in Shestakova et al. (2014) and Shestakova et al. (2016). Note that this function is designed to work only in heteroscedastic narrow evaluation and unstructured mixed models (ie. mHeNE and mHeUN).
}
\examples{
## Calculate within-group heterosdecastic synchrony and SE for conifersIP data:
 data(conifersIP)
 
 #Fit the heteroscedastic set of varcov models (mBE, mHeNE, mHeCS, mHeUN)
 # using taxonomic grouping criteria (i.e. Species)
 ModHt <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
                        data = conifersIP, homoscedastic = FALSE)
 summary(ModHt)
   
 #Obtain the heteroscedastic within-group synchrony and SE for each varGroup stratum.
 gen.het.aSE(ModHt$mHeNE)#Narrow evaluation model
 gen.het.aSE(ModHt$mHeUN)#Unstructured model
 
}
\references{
Shestakova, T.A., Aguilera, M., Ferrio, J.P., Gutierrez, E. & Voltas, J. (2014). Unravelling spatiotemporal tree-ring signals in Mediterranean oaks: a variance-covariance modelling approach of carbon and oxygen isotope ratios. \emph{Tree Physiology} 34: 819-838.

Shestakova, T.A., Gutierrez, E., Kirdyanov, A.V., Camarero, J.J., Genova, M., Knorre, A.A., Linares, J.C., Resco de Dios, V., Sanchez-Salguero, R. & Voltas, J. (2016). Forests synchronize their growth in contrasting Eurasian regions in response to climate warming. \emph{Proceedings of the National Academy of Sciences of the United States of America} 113: 662-667.
}
\seealso{
\code{\link{sync}} for a clear description of synchrony evaluation.
}
\author{
Josu G. Alday, Tatiana A. Shestakova, Victor Resco de Dios, Jordi Voltas
}
