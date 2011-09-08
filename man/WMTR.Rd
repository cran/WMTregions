\name{WMTR}
\alias{WMTR}
\title{
Calculating WMT Regions of the Empirical Data
}
\description{
The function calculates a weighted-mean trimmed (WMT) region of a certain type and data from a specified text file.    
}
\usage{
WMTR(fname = "Cloud.dat", fdir = getwd(), bound = 0)
}
\arguments{
  \item{fname}{
The name of the data file.
}
  \item{fdir}{
A path to the directory where the input and output files should be located. The default value is the \R working directory.
}
  \item{bound}{
\kbd{1} if the upper boundary to be computed, \kbd{-1} if the lower boundary to be computed, \kbd{0} if the whole WMT region to be computed.
}
}
\details{
The format of the input data file:
\enumerate{
\item Type of the trimmed region (zonoid, ECH, ECH*, geometrically trimmed; given weight vector)

\emph{Format:} A text value from the following set: \samp{zonoid}, \samp{ECH}, \samp{ECH*}, \samp{geometrical}, \samp{general}.  \samp{general} is used for the case when the weights are given manually instead of being automatically generated basing on the WMT region type and the depth parameter. 

\item Depth parameter

\emph{Format:} A floating point number from the interval \eqn{[0,\dots,1)}.

\item Dimension

\emph{Format:} An integer number \eqn{d \ge 2}.

\item Number of points of the data cloud

\emph{Format:} An integer number \eqn{n > d}.

\item (\emph{If the type \samp{general} is selected)} The weight vector

\emph{Format:} \eqn{n} non-decreasing floating point numbers matching the requirements for the weight vector.

\item Coordinates of each point

\emph{Format:} \eqn{n} groups of \eqn{d} floating point numbers, each group representing the coordinates of a point from the data cloud. 
}

The points must be in the general position.
}
\value{
\item{\file{TRegion.dat}}{A file in the directory \code{fdir}. The calculated WMT region with facets represented by a normal and an intercept, namely by \eqn{d+1} numbers giving the equation of the hyperplane containing the facet. The first \eqn{d} of these numbers are coordinates of a normal to the facet, which is directed outward the WMT region. The last number defines the intercept.} 
\item{\file{TRegion_vertices.dat}}{A file in the directory \code{fdir}. The calculated WMT region with facets represented by a coordinates of its vertices. A vertex coordinates are grouped by brackets, while vertices are also grouped by brackets into a facet.} 
\item{\file{TRegion_bound.dat}}{An auxiliary file with a boundary of the calculated WMT region. The format is similar to \file{TRegion.dat}.}
\item{\file{TRegion_vertices_bound.dat}}{An auxiliary file with a boundary of the calculated WMT region. The format is similar to \file{TRegion_vertices.dat}.}

}
\references{
Bazovkin, P. and Mosler, K. An Exact Algorithm for Weighted-Mean Trimmed Regions in Any Dimension. 201x. Submitted to JSS.
}
\author{
Pavel Bazovkin and Karl Mosler
}

\seealso{
\code{\link{visualWMTR}}
}

\examples{
\dontrun{generTRsample("Cloud.dat")}     # Generate a sample data file in the working directory
\dontrun{WMTR("Cloud.dat")}              # Processing the file from the working directory
\dontrun{visualWMTR()}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iplot }
\keyword{ programming }
\keyword{ IO }
\keyword{ multivariate }
\keyword{ nonparametric }

