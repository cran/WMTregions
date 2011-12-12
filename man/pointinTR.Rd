\name{pointinTR}
\alias{pointinTR}
\title{
Checking Whether a Point Is In a Specified Trimmed Region
}
\description{

Checks whether a point is in a specified trimmed region.

}
\usage{
pointinTR(dpoint, tregion)
}
\arguments{
  \item{dpoint}{
A vector containing the coordinates of the point to be checked.
}
  \item{tregion}{
A matrix object containing the WMT region in the normal-intercept representation.
}
}
\details{
\code{tregion} is normally produced by the \code{\link{loadWMTR}} function basing on a calculated WMT region.
}
\value{
\item{\code{IsIn}}{\code{TRUE} if the point is in the WMTR, otherwise - \code{FALSE}.}
}
\author{
Pavel Bazovkin and Karl Mosler
}

\seealso{
\code{\link{WMTR}}
\code{\link{loadWMTR}}
}


\examples{
generTRsample("Cloud.dat")
WMTR("Cloud.dat")
tregion <- loadWMTR( )
if(pointinTR( c(1,-2,3), tregion)) "In!" else "Out!"

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }


