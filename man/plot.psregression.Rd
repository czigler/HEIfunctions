\name{plot.psregression}
\alias{plot.psregression}
\title{
Plotting function for Propensity Score Regression.
}
\description{
Makes multiple plots after fitting a regression model.
}
\usage{
plot.psregression(psmodel, data, a, ps, y = NULL, cov = NULL, plot = 1)
}
\arguments{
  \item{psmodel}{model of class psregression.}
  \item{data}{dataset used to fit the model}
  \item{a}{string variable with name of the intervention covariate, coded 0 no intervention, 1 intervention.}
  \item{ps}{string variable with name of the propensity score covariate}
  \item{y}{string variable with name of the outcome covariate}
  \item{cov}{string variable with name of additional covariate for plotting}
  \item{plot}{number of plot that will be displayed, see details}
}
\details{
Plot options: \\
1. Propensity score versus fitted values from PS regression, with no intervention in black, intervention in red.
2. Propensity score versus counterfactual fitted outcomes, with no intervention in black, intervention in red.
3. Plot of cov versus the fitted values, with no intervention in black, intervention in red.
4. Plot of original outcome y versus the fitted values of y, with no intervention in black, intervention in red.
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
