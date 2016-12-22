\name{psbayes}
\alias{psbayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian propensity score regression
}
\description{
Bayesian propensity score regression, no spatial.
}
\usage{
psbayes(formula, data, trt, family = "gaussian", approx = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Model formula within treatment group}
  \item{data}{dataset}
  \item{trt}{string denoting treatment variable}
  \item{family}{works for gaussian or poisson (only if approx = T).}
  \item{approx}{if approx = F, full MCMC.  if approx = T, Gelman's arm package is used.}
  \item{\dots}{pass other options to regression functions}
}
\details{
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
