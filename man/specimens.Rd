\name{specimens}
\alias{specimens}
\docType{data}
\title{
Fictitious data on specimens tested by two methods
}
\description{
stdtest has been done on everyone, and sampledtest has been done on a stratifed subsample of 275 out of 402 specimens (is NA on the other 127 specimens)
}
\usage{data(specimens)}
\format{
  A data frame with 402 observations on the following 3 variables.
  \describe{
    \item{\code{stratum}}{6 strata used for sampling}
    \item{\code{stdtest}}{standard test result available on all specimens}
    \item{\code{sampledtest}}{new test result available only on stratified subsample}
  }
}
\examples{
data(specimens)
}
\keyword{datasets}
