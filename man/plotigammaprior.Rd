\name{plotigammaprior}
\alias{plotigammaprior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize inverse gamma prior
}
\description{
Simple plotting function to visualize the inverse gamma prior for variance parameters.
}
\usage{
plotigammaprior(vec, xlimit = c(0, 100), nsamp = 10000, sd = T)
}
\arguments{
  \item{vec}{vector of length 2, with parameters of inverse gamma}
  \item{xlimit}{limit of x-axis on the plot.}
  \item{nsamp}{number of samples for plot, default is usually okay.}
  \item{sd}{if sd=T, plots the standard deviation.  if sd=F, plots the variance.}
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
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
