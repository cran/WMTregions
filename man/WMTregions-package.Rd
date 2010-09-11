\name{WMTregions-package}
\alias{WMTregions-package}
\alias{WMTregions}
\docType{package}
\title{
Exact calculation of the weighted-mean trimmed regions
}
\description{
This package realizes an exact algorithm (Bazovkin and Mosler, 2010) to compute the weighted-mean trimmed regions of a given data cloud in arbitrary dimension \eqn{d}{d} as convex polytopes in \eqn{R^d}{R^d}. Trimmed regions are a powerful tool of multivariate data analysis. They describe a probability distribution in Euclidean $d$-space regarding location, dispersion and shape, and they order multivariate data with respect to their centrality. In turn, weighted-mean trimmed regions possess attractive properties regarding continuity, subadditivity, and monotonicity. 
}
\details{
\tabular{ll}{
Package: \tab WMTregions\cr
Type: \tab Package\cr
Version: \tab 2.3.10\cr
Date: \tab 2010-09-11\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package contains functions calculating and representing weighted-mean trimmed regions (WMTR). The data cloud is given in a text file. A function \code{generTRsample} can be used for generating a sample data file. A function \code{WMTR} is the main method for calculating WMT regions. Results are stored in a data file and can be visualized by means of \code{showWMTR} function that uses \pkg{rggobi} package. Apply \code{loadWMTR} to load the WMTR to an R variable.  
}
\author{
Pavel Bazovkin, Karl Mosler

Maintainer: Pavel Bazovkin <econometrics@russia.ru> 
}
\references{
Pavel Bazovkin and Karl Mosler. An exact algorithm for weighted-mean trimmed regions in any dimension. 2010. Working paper.

Rainer Dyckerhoff and Karl Mosler. Weighted-mean trimming of multivariate data. 2010. To appear.

Karl Mosler, Tatjana Lange, and Pavel Bazovkin. Computing zonoid trimmed regions in dimension d > 2. Computational Statistics and Data Analysis, 53:2500-2510, 2009.

Karl Mosler. Multivariate Dispersion, Central Regions and Depth: The Lift Zonoid Approach. Springer, New York, 2002.

}

\keyword{ package }
\keyword{ iplot }
\keyword{ programming }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ models }

\examples{

\dontrun{generTRsample("Cloud.dat")}
\dontrun{WMTR("Cloud.dat")}
\dontrun{showWMTR()}

}
