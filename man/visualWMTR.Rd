\name{visualWMTR}
\alias{visualWMTR}
\title{
Visualizing the Calculated WMT Region for the Data in 3D
}
\description{
This function realizes the \eqn{3d}-visualization of the data and computational results of applying \code{\link{WMTR}} function. The parameter \code{fdir} is to be similar to that was used in \code{\link{WMTR}}.}
\usage{
visualWMTR(fdir = getwd())
}
\arguments{
  \item{fdir}{
A path to the directory where the output files of \code{\link{WMTR}} are located.
}
}
\author{
Pavel Bazovkin and Karl Mosler
}
\examples{
\dontrun{generTRsample("Cloud.dat")}     # Generate a sample data file in the working directory
\dontrun{WMTR("Cloud.dat")}              # Processing the file from the working directory
\dontrun{visualWMTR()}
}

\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }
