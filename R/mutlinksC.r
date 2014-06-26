#' mutlinksC: calculate the average mutual links in an adjacency network
#'
#' Network connetctions can be presented as a square adjacency matrix m where element \deqn{m(i,j) = 1} 
#' if there is a link between node \eqn{i} and node \eqn{j}, and \eqn{0} otherwise. Given a square adjacency matrix m, 
#' mutlinksC will compute and return the average mutual links between all nodes in the network.
#' 
#' \code{mutlinksC} has the capacity to use OpenMP for parallelization. The optional parameter \code{num_threads} represents the number
#' thread processes to be started for a given run.
#' @param adj A square adjecency matrix representing links in a network
#' @param num_threads The number of threads to process the calculation with. It is recommended that this number
#' is never larger than the number of processors on a machine.
#' @keywords multilinks, adjacency matrix, OpenMP
#' @family fastpackage
#' @useDynLib fastpackage
#' @export
#' @examples
#' # simple run
#' n <- 100 # square matrix size
#' m <- matrix(sample(0:1,n^2,replace=T),nrow=n)  # create adjacency matrix of 0 and 1s
#' mutlinksC(m)
#' # test threads
#' for (num_thread in 1:4) {
#'      t = system.time(mutlinksC(m,num_thread))
#'      print(t)
#' }
mutlinksC <- function(adj, num_threads=1) {
  # check to see if symetric and matrix
  if (!is.loaded("mutlinks")) {
    stop("Function mutlinks not loaded from C")
  } else if (dim(adj)[1] != dim(adj)[2]) {
    stop("Input must be square.")
  } else if (num_threads <= 0 || num_threads != floor(num_threads)) {
    stop("The Number of Threads parameter 'num_threads' must be strictly positive and an integer.")
  }
  
  # check to see if data frame and convert to matrix
  if(!is.matrix(adj)){
    adj = as.matrix(adj)
  }
  
  # check to see if all values of matrix are strictly 0 or 1
  if ( (sum(adj %in% c(0,1)) != (dim(adj)[1]*dim(adj)[2])) ) {
    warning("Input adjacency matrix has entries other than 0 or 1. Changing all non-zero entries to 1.")
    adj[adj != 0] = 1
  }
  
  # setting OpenMP parameter OMP_NUM_THREADS
  # Sys.setenv(OMP_NUM_THREADS=num_threads)
  
  n = dim(adj)[1] # dimention of matrix for C
  z <- .C("mutlinks",
          as.integer(adj),
          as.integer(n),
          result=double(1), 
          threads=as.integer(num_threads))
  
  return(z$result)
}