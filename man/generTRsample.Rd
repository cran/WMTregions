\name{generTRsample}
\alias{generTRsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generating a sample data file for calculating WMTR
}
\description{
This function creates a uniformly distributed pseudorandom data cloud according to the given parameters.
}
\usage{
generTRsample(fname = "Cloud.dat", fdir = getwd(), dim = 3, num = 10, alpha = 0.221, trtype = "zonoid")
}
\arguments{
  \item{fname}{
A name of a target file
}
  \item{fdir}{
A directory of a target file (by default, the working directory)
}
  \item{dim}{
Dimension of the data cloud
}
  \item{num}{
Number of data points in the data cloud
}
  \item{alpha}{
Depth parameter \code{alpha} - floating point number in [0,1[, e.g. 0.221  
}
  \item{trtype}{
WMTR type - string value (choose from \var{zonoid}, \var{ECH}, \var{ECH*}, \var{geometrical})
}
}
\details{
}
\value{
The format of a data file:

\enumerate{
\item WMTR type - string value (choose from \var{zonoid}, \var{ECH}, \var{ECH*}, \var{geometrical})
\item Depth parameter \var{alpha} - floating point number in [0,1[, e.g. 0.221
\item Dimension \var{dim} - integer number > 2
\item Number of data points \var{num} - integer number > \var{dim}
\item Serially written coordinates of the data points - floating point numbers
}
}
\author{
Pavel Bazovkin and Karl Mosler
}

\seealso{
\code{\link{WMTR}}
}

\examples{
\dontrun{generTRsample("Cloud.dat", dim = 3, num = 10, alpha = 0.221, trtype = "zonoid")}
\dontrun{WMTR("Cloud.dat")}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

