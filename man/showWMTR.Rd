\name{showWMTR}
\alias{showWMTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualizing a calculated weighted-mean trimmed region
}
\description{
Visualizes a calculated WMTR as convex polytope by representing its vertices in the \pkg{rggobi} interactive graphics framework. The WMTR should be firstly calculated. 
}
\usage{
showWMTR(fdir = getwd())
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

\seealso{
\code{\link{WMTR}}
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

