\name{subdiagC}
\alias{subdiagC}
\title{subdiagC Sub-diagonal Function}
\usage{
subdiagC(m, k)
}
\arguments{
  \item{m}{The input square matrix or data frame}

  \item{k}{The index of the sub-diagonal to return (i.e., 0
  for the diagonal, 1 for the first subdiagonal, 2 for the
  second, etc.)}
}
\description{
Given a square matrix or data frame, and the paramter k,
this function returns the kth subdiagonal vector of the
matrix. Implemented in C for extra speed.
}
\examples{
m <- matrix(1:25,5,5) # create 5x5 matrix
k = 2
subdiagC(m,k)
}
\seealso{
Other fastpackage: \code{\link{cat_function}};
\code{\link{fastLSA}}; \code{\link{fastpackage}},
\code{\link{fastpackage-package}}; \code{\link{mutlinksC}}
}
\keyword{demo}
\keyword{matrix,}
\keyword{square}
\keyword{sub-diagonal,}

