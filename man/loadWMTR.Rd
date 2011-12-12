\name{loadWMTR}
\alias{loadWMTR}
\title{
Loading a Calculated WMT Region From a File to a Matrix Object
}
\description{
This function loads the calculated WMT region into a matrix object in order to work with it as with a variable in \R. For example, for applying the function \code{\link{pointinTR}}.}
\usage{
loadWMTR(fname = "TRegion.dat", fdir = getwd())
}
\arguments{
  \item{fname}{
The name of the file with the calculated WMT region (the normal-intercept representation output file of \code{\link{WMTR}}) in the directory \code{fdir}.
}
  \item{fdir}{
A path to the directory where the output files of \code{\link{WMTR}} are located. If empty, the working directory is considered.
}

}

\value{
\item{\code{tregion}}{A matrix object containing the normal-intercept coordinates of the WMT region facets as its rows.}}
\author{
Pavel Bazovkin and Karl Mosler
}

\seealso{
\code{\link{WMTR}}
\code{\link{pointinTR}}
}


\examples{
generTRsample("Cloud.dat")
WMTR("Cloud.dat")
tregion <- loadWMTR( )
tregion
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

