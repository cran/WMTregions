\name{WMTR}
\alias{WMTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating WMTR from the empirical data given in a text file
}
\description{
The function calculates WMTR of a defined type and data from a specified text file.    
}
\usage{
WMTR(fname = "Cloud.dat", fdir = getwd())
}
\arguments{
  \item{fname}{
String argument describing a name of the data file.
}
  \item{fdir}{
String argument describing a name of a directory where the data file is located. If empty, the file is taken from the working directory.
}
}
\details{
The format of a data file:

1. WMTR type - string value (choose from \var{zonoid}, \var{ECH}, \var{ECH*}, \var{geometrical})

2. Depth parameter \var{alpha} - floating point number in [0,1[, e.g. 0.221

3. Dimension \var{dim} - integer number > 2

4. Number of data points \var{num} - integer number > \var{dim}

5. Serially written coordinates of the data points - \var{dim * var} floating point numbers

}
\value{
Computational results are stored in the working directory in a file "TRegion.dat". Format: each line represents a facet of the trimmed region, where first \var{dim} floating point numbers are coordinates of a normal to the facet and the last one - its absolute distance from the origin. 
}
\references{
Pavel Bazovkin and Karl Mosler. An exact algorithm for weighted-mean trimmed regions in any dimension. 2010. Working paper.
}
\author{
Pavel Bazovkin and Karl Mosler
}

\examples{
\dontrun{generTRsample("Cloud.dat")}     # Generate a sample data file in the working directory
\dontrun{WMTR("Cloud.dat")}              # Processing the file from the working directory
\dontrun{showWMTR()}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

