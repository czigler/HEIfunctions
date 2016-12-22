\name{plot.ps.spatmodel}
\alias{plot.ps.spatmodel}
\title{
Plotting function for propensity score spatial model
}
\description{
Plots output from a ps.spatmodel regression.
}
\usage{
plot.ps.spatmodel(psmod, plot = 1, group = NULL, contour = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psmod}{object of class ps.spatmodel}
  \item{plot}{number of plot you would like displayed, see details.}
  \item{group}{if group = NULL, plots for both groups are shown.  If group = 1, only intervention group is shown.  If group = 0, 
		   only no intervention group.}
  \item{contour}{if T, contours are displayed on plots.}
}
\details{
The different plot numbers correspond to:
1. Extrapolated surface from observed data only. \\
2. Estimated random effects. \\
3. Predicted surface.  
}
\value{
}
\references{
}
\author{
Lauren
}
\note{
}
\seealso{
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
