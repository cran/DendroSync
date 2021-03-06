% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bet.het.aSE.R
\name{bet.het.aSE}
\alias{bet.het.aSE}
\title{Between-group synchrony for a heteroscedastic unstructured model}
\usage{
bet.het.aSE(model)
}
\arguments{
\item{model}{a class "\code{lme}" unstructured model (\code{mHeUN}) produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{FALSE}.}
}
\value{
The function returns a \code{matrix} containing between-group synchrony and SE for each combination of \code{varGroup} levels. This function is used internally in \code{\link{sync}}.
}
\description{
The function calculates the between-group synchrony (a^) and standard error (SE) for a heteroscedastic unstructured model (mHeUN).
}
\details{
The function calculates between-group synchrony for a heteroscedastic unstructured model (mHeUN).
}
\examples{
## Calculate between-group synchrony and SE 
 ##for a heteroscedastic unstructured model for conifersIP data:
 data(conifersIP)
 
 #Fit the heteroscedastic set of varcov models (mBE, mHeNE, mHeCS, mHeUN)
 # using taxonomic grouping criteria (ie. Species)
 ModHt <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
                        data = conifersIP, homoscedastic = FALSE)
 
 #between-group synchrony and SE for each varGroup stratum combination
 # in heteroscedastic unstructured models.
 bet.het.aSE(ModHt$mHeUN)#Unstructured model
 
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
