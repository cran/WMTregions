\name{pointinTR}
\alias{pointinTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Checking whether a point is in a specified trimmed region
}
\description{

Checks whether a point is in a specified trimmed region

}
\usage{
pointinTR(dpoint, tregion)
}
\arguments{
  \item{dpoint}{
Floating point number vector - coordinates of the point to be tested
}
  \item{tregion}{
Matrix representing a WM trimmed region
}
}
\details{
}
\value{
Boolean: TRUE if the point is in the WMTR, otherwise - FALSE
}
\author{
Pavel Bazovkin and Karl Mosler
}

\seealso{
\code{\link{WMTR}}
}


\examples{
\dontrun{generTRsample("Cloud.dat")}
\dontrun{WMTR("Cloud.dat")}
\dontrun{tregion <- loadWMTR( )}
\dontrun{if(pointinTR( c(1,-2,3), tregion)) "In!" else "Out!"}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }


