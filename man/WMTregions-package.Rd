\name{WMTregions-package}
\alias{WMTregions-package}
\alias{WMTregions}
\docType{package}
\title{
Exact Calculation of the Weighted-Mean Trimmed Regions (WMTR)
}
\description{
This package realizes an exact algorithm (Bazovkin and Mosler, 201x) to compute the weighted-mean trimmed regions of a given data cloud in arbitrary dimension \eqn{d} as convex polytopes in \eqn{R^d}. Trimmed regions are a powerful tool of multivariate data analysis. They describe a probability distribution in Euclidean \eqn{d}-space regarding location, dispersion and shape, and they order multivariate data with respect to their centrality. In turn, weighted-mean trimmed regions possess attractive properties regarding continuity, subadditivity, and monotonicity. 
}
\details{
\tabular{ll}{
Package: \tab WMTregions\cr
Type: \tab Package\cr
Version: \tab 3.2\cr
Date: \tab 2011-12-12\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package contains functions calculating and representing the weighted-mean trimmed regions (WMTR). The data cloud is given in a text file. A function \code{\link{generTRsample}} can be used for generating a sample data file. A function \code{\link{WMTR}} is the main method for calculating WMT regions. Results are stored in data files and can be visualized by means of either \code{\link{showWMTR}} function that uses \pkg{rggobi} package or \code{\link{visualWMTR}} method producing \acronym{OpenGL} visualization. Apply \code{\link{loadWMTR}} to load the WMTR into an \R variable.  
}
\author{
Pavel Bazovkin, Karl Mosler

Maintainer: Pavel Bazovkin \email{econometrics@russia.ru} 
}

\seealso{
\code{\link{WMTR}}
\code{\link{visualWMTR}}
\code{\link{showWMTR}}
\code{\link{pointinTR}}
\code{\link{loadWMTR}}
\code{\link{generTRsample}}
}


\references{
Bazovkin, P. and Mosler, K. (2012). An Exact Algorithm for Weighted-Mean Trimmed Regions in Any Dimension. Journal of Statistical Software, \bold{47(13)}, pp. 1--29. URL \url{http://www.jstatsoft.org/v47/i13/}.

Dyckerhoff, R. and Mosler, K. (2011). Weighted-mean trimming of multivariate data. Journal of Multivariate analysis, \bold{102}, pp. 405--421.

Mosler, K. and Lange, T. and Bazovkin, P. (2009). Computing zonoid trimmed regions in dimension d > 2. Computational Statistics and Data Analysis, \bold{53}, pp. 2500--2510.

Mosler, K. (2002). \emph{Multivariate Dispersion, Central Regions and Depth: The Lift Zonoid Approach}. New York: Springer-Verlag.

}

\keyword{ package }
\keyword{ iplot }
\keyword{ programming }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ models }

\examples{

generTRsample("Cloud.dat")
WMTR("Cloud.dat")
visualWMTR()

}
