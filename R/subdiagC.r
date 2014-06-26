#' subdiagC Sub-diagonal Function
#'
#' Given a square matrix or data frame, and the paramter k, this function returns
#' the kth subdiagonal vector of the matrix. Implemented in C for extra speed.
#' @param m The input square matrix or data frame
#' @param k The index of the sub-diagonal to return (i.e., 0 for the diagonal, 1 for
#' the first subdiagonal, 2 for the second, etc.)
#' @keywords sub-diagonal, square matrix, demo
#' @family fastpackage
#' @useDynLib fastpackage
#' @export
#' @examples
#' m <- matrix(1:25,5,5) # create 5x5 matrix
#' k = 2
#' subdiagC(m,k)
subdiagC <- function (m, k){
  # check to see if symetric and matrix
  if (!is.loaded("subdiag")) {
    stop("ERROR: Function subdiag not loaded from C")
  } else if (dim(m)[1] != dim(m)[2]) {
    stop("ERROR: Input must be square.")
  } else if (length(k) != 1 || k != floor(k) || k < 0 || k >= dim(m)) {
    stop("Subdiagonal index k must be an integer > 0 and less than the size of matrix m.")
  }
  # check to see if data frame and convert to matrix
  if(!is.matrix(m)){
    m = as.matrix(m)
  }
  
  # call C: subdiag.cc
  temp = .C("subdiag", 
     as.double(m), 
     as.integer(dim(m)[1]), 
     as.integer(k), 
     result = double(dim(m)[1]-k) )
  # return sub-diag
  return(temp$result)
}