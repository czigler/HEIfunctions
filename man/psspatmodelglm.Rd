\name{psspatmodelglm}
\alias{psspatmodelglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial propensity score model for Poisson outcomes.
}
\description{
Spatial propensity score model for Poisson outcomes.
}
\usage{
psspatmodelglm(formula, data, trt, coords, family = "poisson", na.rm = T, nsamp, nburn, tau.sq.IG, sigma.sq.IG, tau.sq.start = 1, sigma.sq.start = 1, w.start = 0, phi.tune = 0.01, sigma.sq.tune = 0.15, tau.sq.tune = 0.15, beta.tune = rep(0.1, 5), w.tune = 0.1, phiprior = NULL)
}
\arguments{
  \item{formula}{
formula
}
  \item{data}{
dataset
}
  \item{trt}{
string variable with name of intervention covariate, coded as 0 or 1.
}
  \item{coords}{
coordinates corresponding to data
}
  \item{family}{
right now, only works with family = "poisson"
}
  \item{na.rm}{
if na.rm=T, outcomes with missingness are deleted, and only counterfactuals are imputed.  if = F, missing outcomes are imputed.
}
  \item{nsamp}{
total number of samples
}
  \item{nburn}{
number of samples to discard as burnin
}
  \item{tau.sq.IG}{
prior for tau.sq
}
  \item{sigma.sq.IG}{
prior for sigma.sq
}
  \item{tau.sq.start}{
starting value for tau.sq
}
  \item{sigma.sq.start}{
starting value for sigma.sq
}
  \item{w.start}{
starting value for w
}
  \item{phi.tune}{
tuning parameter for phi
}
  \item{sigma.sq.tune}{
tuning parameter for sigma.sq
}
  \item{tau.sq.tune}{
tuning parameter for tau.sq
}
  \item{beta.tune}{
vector of tuning parameters for beta
}
  \item{w.tune}{
tuning parameter for w
}
  \item{phiprior}{
prior for phi
}
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


\seealso{
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
