\name{psregression}
\alias{psregression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propensity Score Regression
}
\description{
Frequentist propensity score regression
}
\usage{
psregression(formula, data, trt, nboot = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{model formula (should be a function of the propensity score)}
  \item{data}{dataset}
  \item{trt}{string variable denoting treatment indicator in dataset}
  \item{nboot}{number of bootstrap samples}
  \item{\dots}{other arguments to pass to glm}
}
\details{
known bug: I exclude any bootstrap samples that produce an error.  this happens when, for blocked PS model, bootstrap sample
is produced without all PS categories.
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
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
