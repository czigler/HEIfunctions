\name{mergepredictions}
\alias{mergepredictions}
\title{
Merge baseline predictions with original data
}
\description{
}
\usage{
mergepredictions(dat, preddat, cov, predcov, id, predid = NULL, name = NULL)
}
\arguments{
  \item{dat}{Original dataset}
  \item{preddat}{Dataset with predictions}
  \item{cov}{Name of covariate in original dataset}
  \item{predcov}{Name of covariate in predicted dataset}
  \item{id}{Name of unique identifier in original dataset that links the predicted and original data}
  \item{predid}{Name of unique identifier in predicted dataset that links the predicted and original data (uses id if not specified)}
  \item{name}{Name of new covariate that contains predicted and original data}
}
\details{
}
\value{
Returns data frame identical to original dataset, but with missing baseline filled in with predictions.
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
