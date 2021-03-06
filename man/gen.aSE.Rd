% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen.aSE.R
\name{gen.aSE}
\alias{gen.aSE}
\title{Within-group synchrony for homoscedastic models}
\usage{
gen.aSE(model)
}
\arguments{
\item{model}{a class "\code{lme}" model produced by \code{\link{dendro.varcov}} with \code{homoscedastic} equals \code{TRUE}.}
}
\value{
The function returns a \code{matrix} containing the synchrony and SE for each level of \code{varGroup}. This function is used internally in \code{\link{sync}}.
}
\description{
The function calculates for each varGroup stratum the within-group synchrony (a^) and standard error (SE). However, it only works for homoscedastic broad evaluation, narrow evaluation and unstructured models (mBE, mNE, mUN).
}
\details{
The function calculates the within-group synchrony values for each varGroup stratum based on the metodology described in Shestakova et al. (2014) and Shestakova et al. (2016). Note that this function is designed to work only in 3 homoscedastic models (mBE, mNE, mUN).
}
\examples{
## Calculate within-group homoscedastic synchrony and SE for conifersIP data:
 data(conifersIP)
 
 #Fit the homoscedastic set of varcov models (mBE, mNE, mCS, mUN) 
 # using taxonomic grouping criteria (ie. Species)
 ModHm <- dendro.varcov(TRW ~ Code, varTime = "Year", varGroup = "Species", 
                         data = conifersIP, homoscedastic = TRUE)
 summary(ModHm)
   
 #Obtain the within-group synchrony and SE for each varGroup stratum.
 gen.aSE(ModHm$mBE)#Broad evaluation
 gen.aSE(ModHm$mUN)#Unstructured
 
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
