\name{summary.spatpredmodel}
\alias{summary.spatpredmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize spatial prediction model
}
\description{
Summarizes the baseline spatial prediction model.
}
\usage{
summary.spatpredmodel(outmodel, covname = NULL, origdat, preddat, full = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outmodel}{spatial prediction model}
  \item{covname}{string with outcome name}
  \item{origdat}{original dataset}
  \item{preddat}{prediction dataset}
  \item{full}{if full = T, summary from spBayes original model is shown.  if full = F, only prediction summary is shown.}
}
\details{
}
\value{
Outputs table for easy input into LaTeX via xtable.
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
