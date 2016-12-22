\name{plotps}
\alias{plotps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize propensity score.
}
\description{
Multiple plotting functions for helping visualize the propensity score.
}
\usage{
plotps(data, plot = 1, ps, a, basepollution = NULL, outcome = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{dataset}
  \item{plot}{number of plot, see details}
  \item{ps}{string with name of propensity score covariate, e.g. "ps"}
  \item{a}{string with name of intervention covariate, 1 intervention, 0 no intervention}
  \item{basepollution}{string with name of baseline pollution covariate (used with plots 3 and 4)}
  \item{outcome}{string with name of outcome covariate (used with option 4)}
}
\details{
1. Plot of propensity score (indexed by observation), intervention group in red.
2. Histogram of propensity score for intervention and no intervention groups. 
3. Plot of propensity score versus baseline pollution, intervention group in red.
4.  Plot of propensity score versus baseline pollution, highlighting where the missing outcome data is. 
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
