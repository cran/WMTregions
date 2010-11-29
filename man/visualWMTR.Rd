\name{visualWMTR}
\alias{visualWMTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3d-visualization of the calculated weighted-mean trimmed region
}
\description{
The function uses OpenGL to visualize a WMTR for the data of dimension 3. 
}
\usage{
visualWMTR(fdir = getwd())
}
\arguments{
  \item{fdir}{
String argument describing a name of a directory where the computing results are located. If empty, the working directory is considered.
}
}
\details{
}
\value{
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
