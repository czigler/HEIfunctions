\name{mvpsmod}
\alias{mvpsmod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate linear spatial propensity score model.
}
\description{
Multivariate linear spatial propensity score model.  Jointly models the potential outcomes, rather than separate linear models.
}
\usage{
mvpsmod(formula, data, trt, coords, nsamp, nburn, thin, tuning = list(A = 0.1, psi = 0.2, theta = 1), prior = list(KIG = rep(0.01, 2), psi = rep(0.01, 2), theta1 = rep(0.6, 2), theta2 = rep(10, 2)), starting = list(B = rep(0, 2 * length(all.vars(formula))), A = NULL, psi = NULL, theta = NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{function of propensity score.  do not include trt.}
  \item{data}{dataset}
  \item{trt}{string with name of treatment variable; must be defined as 1 if intervention, 0 if no intervention}
  \item{coords}{spatial coordinates corresponding to data.}
  \item{nsamp}{total number of samples}
  \item{nburn}{number of samples to discard as burn in}
  \item{thin}{thinning of samples}
  \item{tuning}{tuning parameters}
  \item{prior}{priors}
  \item{starting}{starting values}
}
\details{
}
\value{
\item{samples}{MCMC samples for parameters}
\item{y0}{MCMC samples for y0 (if y0 known, it is constant)}
\item{y0}{MCMC samples for y1 (if y1 known, it is constant)}
\item{coords}{spatial coordinates corresponding to data}
\item{trt}{vector indicating whether intervention or no intervention}
}
\references{
}
\author{
Cory 
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
