\name{fastLSA}
\alias{fastLSA}
\title{fastLSA: Efficent local simmilarity analysis}
\usage{
fastLSA(mat, d = 3, alpha = 0.05, minLSA = 0, rez = 1e+05,
  num_threads = 1)
}
\arguments{
  \item{mat}{An m-by-n matrix with the m columns
  representing observations accross n time-series column
  steps}

  \item{d}{A lead-lag window to search for out-of-phase,
  leading or lagged correlations. By default this is 3.}

  \item{alpha}{The significance level \eqn{\alpha} for each
  LSA value. fastLSA uses approximate p-value to
  efficiently calculate the upper bound. If this
  theoretical upper bound is less than the stated
  confidence level \eqn{(1-\alpha)}, the test is not
  reported. The default is \eqn{\alpha=0.05} or a
  confidence level of \eqn{(1-\alpha)=95%}.}

  \item{minLSA}{The minimum LSA value to report. The LSA
  similarity statistic has a range of \eqn{0.0-1.0}, with
  0.0 representing no correlation and 1.0 high correlation.
  By default all correlations are reported, i.e.
  \code{minLSA=0}.}

  \item{rez}{The p-value resolution. This parameter
  represents the quality of the p-value as calcuated by a
  Riemann integral. \code{rez} represents the number of
  equally spaced time-steps per standard devidation. By
  default this value is 100000.}

  \item{num_threads}{The number of threads to process the
  calculation. It is recommended that this number is never
  larger than the number of processors on a machine.}
}
\description{
Drives a parallel version of fastLSA
}
\details{
\code{fastLSA} has the capacity to use OpenMP for
parallelization. The optional parameter \code{num_threads}
represents the number thread processes to be started for a
given run.
}
\examples{
print('fastLSA')
}
\seealso{
Other fastpackage: \code{\link{cat_function}};
\code{\link{fastpackage}},
\code{\link{fastpackage-package}}; \code{\link{mutlinksC}};
\code{\link{subdiagC}}
}
\keyword{fastLSA}

