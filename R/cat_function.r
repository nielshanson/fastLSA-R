#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @family fastpackage
#' @export
#' @examples
#' cat_function()
#' cat_function(love=FALSE)
cat_function <- function(love=TRUE) {
  if (love==TRUE) {
    print ("I love cats!")
  } else {
    print("I don't love cats!")
  }
}