\name{generTRsample}
\alias{generTRsample}
\title{
Generating Sample Data Cloud File For Applying \code{\link{WMTR}}
}
\description{
This function is an auxiliary one. It generates a random uniformly distributed data cloud of any size \code{num} and any dimension \code{dim} with the format appropriate for applying \code{\link{WMTR}}. }
\usage{
generTRsample(fname = "Cloud.dat", fdir = getwd(), dim = 3, num = 20, alpha = 0.05, trtype = "zonoid")
}
\arguments{
  \item{fname}{
A name of the target file.
}
  \item{fdir}{
A directory of a target file (by default, the working directory).
}
  \item{dim}{
Dimension of the data cloud.
}
  \item{num}{
Number of data points in the data cloud.
}
  \item{alpha}{
Depth parameter \code{alpha} - floating point number in \eqn{(0,\dots,1)}, e.g. \kbd{0.117}.  
}
  \item{trtype}{
The notion of the WMT region to be calculated.
}
}
\value{
\item{A file \code{fname}}{An input file for \code{\link{WMTR}}. See \code{\link{WMTR}} for the format details.}
}


\author{
Pavel Bazovkin and Karl Mosler
}

\seealso{
\code{\link{WMTR}}
}

\examples{
generTRsample("Cloud.dat", dim = 3, num = 20, alpha = 0.221, trtype = "zonoid")
WMTR("Cloud.dat")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

