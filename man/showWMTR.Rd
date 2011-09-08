\name{showWMTR}
\alias{showWMTR}
\title{
Demonstrating the Calculated WMT Region of Any Dimension
}
\description{
Demonstrates the calculated WMT region of any dimension by making multiple projections of its vertices into 3D. 
}
\usage{
showWMTR(fdir = getwd())
}
\arguments{
  \item{fdir}{
A path to the directory where the output files of \code{\link{WMTR}} are located. Must be similar to \code{fdir} that was used in \code{\link{WMTR}}.
}
}

\details{
The function visualizes a calculated WMT region as a convex polytope by representing its vertices in \pkg{rggobi}. The visualization is a series of projections into the space of dimension 3. In comparison with \code{\link{visualWMTR}}, \code{\link{showWMTR}} visualizes only vertices but, however, does it for the data of any dimension.
}
\author{
Pavel Bazovkin and Karl Mosler
}

\seealso{
\code{\link{WMTR}}
\code{\link{visualWMTR}}
}


\examples{
\dontrun{generTRsample("Cloud.dat")}
\dontrun{WMTR("Cloud.dat")}
\dontrun{showWMTR()}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

