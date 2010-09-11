\name{loadWMTR}
\alias{loadWMTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loading a calculated weighted-mean trimmed region from a text file to an R object
}
\description{
Having calculated a WMT region you can load it into a matrix for working with it. 
}
\usage{
loadWMTR(dim, fdir = getwd())
}
\arguments{
  \item{dim}{
Integer number representing affine dimension of the trimmed region. 
}
  \item{fdir}{
String argument describing a name of a directory where the computing results are located. If empty, the working directory is considered.
}

}
\details{

}

\value{
Matrix with each row representing a facet of the trimmed region
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
\dontrun{showWMTR()}
\dontrun{tregion <- loadWMTR( 3 )}
\dontrun{tregion}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

