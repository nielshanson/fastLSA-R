#' fastLSA: Efficent local simmilarity analysis
#'
#' Drives a parallel version of fastLSA 
#' 
#' \code{fastLSA} has the capacity to use OpenMP for parallelization. The optional parameter \code{num_threads} represents the number
#' thread processes to be started for a given run.
#' @param mat An m-by-n matrix with the m columns representing observations accross n time-series column steps
#' @param d A lead-lag window to search for out-of-phase, leading or lagged correlations. By default this is 3.
#' @param alpha The significance level \eqn{\alpha} for each LSA value. fastLSA uses approximate p-value
#' to efficiently calculate the upper bound. If this theoretical upper bound is less than the stated confidence level \eqn{(1-\alpha)}, 
#' the test is not reported. The default is \eqn{\alpha=0.05} or a confidence level of \eqn{(1-\alpha)=95%}.
#' @param minLSA The minimum LSA value to report. The LSA similarity statistic has a range of \eqn{0.0-1.0}, with
#' 0.0 representing no correlation and 1.0 high correlation. By default all correlations are reported, i.e. \code{minLSA=0}.
#' @param rez The p-value resolution. This parameter represents the quality of the p-value as calcuated by
#' a Riemann integral. \code{rez} represents the number of equally spaced time-steps per standard devidation.
#' By default this value is 100000.
#' @param num_threads The number of threads to process the calculation. It is recommended that this number
#' is never larger than the number of processors on a machine.
#' @keywords fastLSA
#' @family fastpackage
#' @useDynLib fastpackage
#' @export
#' @examples
#' print('fastLSA')
fastLSA <- function(mat, d=3, alpha=0.05, minLSA=0, rez=100000, num_threads=1) {
  # check matrix
  if (!is.matrix(mat)) {
    mat = as.matrix(mat)
  }
  m = dim(mat)[1] # num rows (variables)
  n = dim(mat)[2] # num columns (time-steps)
  
  # check parameters
  # d
  if(d < 0 || floor(d) != d) {
    stop("Offset parameter d must be non-negative integer ")
  } else if (d > n) {
    stop("d can not be greater than the number of columns of n")
  }
  # alpha
  if(alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1")
  }
  # minLSA
  if(minLSA < 0 || minLSA > 1) {
    stop("minLSA must be between 0 and 1")
  }
  # rez
  if(rez < 1 || floor(d) != d) {
    stop("rez must be a strictly positive integer")
  }
  # num_threads
  if(num_threads <= 0 || num_threads != floor(num_threads)) {
    stop("num_threads must be a strictly positive integer")
  }
  
  # call C code
  z <- .C("lsarun",
          as.double(t(mat)),
          as.integer(m),
          as.integer(n),
          as.integer(d),
          as.double(alpha),
          as.double(minLSA),
          as.integer(rez),
          lsa_result=double(m*(m-1)),
          pval_result=double(m*(m-1)),
          threads=integer(1))
  return(z$lsa_result)
}