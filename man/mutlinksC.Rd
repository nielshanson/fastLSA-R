\name{mutlinksC}
\alias{mutlinksC}
\title{mutlinksC: calculate the average mutual links in an adjacency network}
\usage{
mutlinksC(adj, num_threads = 1)
}
\arguments{
  \item{adj}{A square adjecency matrix representing links
  in a network}

  \item{num_threads}{The number of threads to process the
  calculation with. It is recommended that this number is
  never larger than the number of processors on a machine.}
}
\description{
Network connetctions can be presented as a square adjacency
matrix m where element \deqn{m(i,j) = 1} if there is a link
between node \eqn{i} and node \eqn{j}, and \eqn{0}
otherwise. Given a square adjacency matrix m, mutlinksC
will compute and return the average mutual links between
all nodes in the network.
}
\details{
\code{mutlinksC} has the capacity to use OpenMP for
parallelization. The optional parameter \code{num_threads}
represents the number thread processes to be started for a
given run.
}
\examples{
# simple run
n <- 100 # square matrix size
m <- matrix(sample(0:1,n^2,replace=T),nrow=n)  # create adjacency matrix of 0 and 1s
mutlinksC(m)
# test threads
for (num_thread in 1:4) {
     t = system.time(mutlinksC(m,num_thread))
     print(t)
}
}
\seealso{
Other fastpackage: \code{\link{cat_function}};
\code{\link{fastLSA}}; \code{\link{fastpackage}},
\code{\link{fastpackage-package}}; \code{\link{subdiagC}}
}
\keyword{OpenMP}
\keyword{adjacency}
\keyword{matrix,}
\keyword{multilinks,}

