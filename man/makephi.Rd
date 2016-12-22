\name{makephi}
\alias{makephi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Proposes a reasonable value for range parameter
}
\description{
Proposes a reasonable value for range parameter using method recommended in Banerjee et al. book.
}
\usage{
makephi(coords, scale = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{coordinates
}
  \item{scale}{denominator for fraction of distance at which you want correlation between two points to decay 
below 0.05 (e.g. 10 would correspond to 1/10 of the maximum distance between two points)}
}
\details{Silly, one line function but I always forget how the range parameter is parameterized in R, so I made this a function.
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
