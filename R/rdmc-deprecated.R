# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************

#' Deprecated functions in package \pkg{rdmc}
#' 
#' The following functions are deprecated and may be removed as soon as the 
#' next version.
#' 
#' @name rdmc-deprecated
#' 
#' @keywords internal

NULL


#' @rdname rdmc-deprecated
#' 
#' @section Replacement functions:
#' Instead of \code{get_X()}, use \code{\link{get_completed}()} or 
#' \code{\link{get_imputed}()}.
#' 
#' @export

get_X <- function(object, ...) {
  .Deprecated(
    msg = c("'get_X' is deprecated.\n", 
            "Use 'get_completed' or 'get_imputed' instead.\n",
            "See help(\"Deprecated\") and help(\"rdmc-deprecated\").")
  )
  UseMethod("get_completed")
}
