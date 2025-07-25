# ************************************
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ************************************


#' Extract the completed (imputed) data matrix
#' 
#' Extract the completed (i.e., imputed) data matrix from an object returned by 
#' a matrix completion algorithm. 
#' 
#' @param object  an object returned by a matrix completion algorithm.
#' @param which  an integer specifying the index of the regularization 
#' parameter for which to extract the completed data matrix.
#' @param discretized  a logical indicating if the completed data matrix with 
#' or without the discretization step should be extracted. The default is 
#' \code{TRUE} if the discretization step was performed and \code{FALSE} 
#' otherwise.
#' @param \dots  additional arguments to be passed down to methods.
#' 
#' @return  The completed (i.e., imputed) data matrix.
#' 
#' @note Matrix completion and imputation are synonymous terms used in 
#' different streams of the literature, hence \code{get_imputed()} is an 
#' alias for \code{get_completed()} with the same functionality.
#' 
#' @author Andreas Alfons and Aurore Archimbaud
#' 
#' @seealso \code{\link{rdmc_tune}()}, \code{\link{soft_impute_tune}()}, 
#' \code{\link{median_impute}()}, \code{\link{mode_impute}()}
#' 
#' @examples
#' # toy example derived from MovieLens 100K dataset
#' data("MovieLensToy")
#' # robust discrete matrix completion with hyperparameter tuning
#' set.seed(20250723)
#' fit <- rdmc_tune(MovieLensToy, 
#'                  lambda = fraction_grid(nb_lambda = 6),
#'                  splits = holdout_control(R = 5))
#' # extract completed matrix with optimal regularization parameter
#' X_hat <- get_completed(fit)
#' head(X_hat)
#' 
#' # for more examples, see the help files of other functions for 
#' # matrix completion and imputation methods
#' 
#' @keywords utilities
#' 
#' @export

get_completed <- function(object, ...) UseMethod("get_completed")


#' @rdname get_completed
#' @export
get_completed.rdmc_tuned <- function(object, ...) {
  object$fit$X
}


#' @rdname get_completed
#' @export
get_completed.rdmc <- function(object, which, ...) {
  out <- object$X
  if (length(object$lambda) > 1L) {
    if (missing(which)) {
      stop("multiple regularization parameters used; use argument 'which' ", 
           "to specify the index for which to extract the completed matrix")
    } else if (!is.numeric(which) || length(which) != 1L) {
      stop("argument 'which' must be a single integer")
    }
    out <- out[[which]]
  }
  out
}

#' @rdname get_completed
#' @export
get_completed.soft_impute_tuned <- function(object, discretized = NULL, ...) {
  have_discretized <- !is.null(object$fit$X_discretized)
  if (is.null(discretized)) discretized <- have_discretized
  else {
    discretized <- isTRUE(discretized)
    if (discretized & !have_discretized) 
      stop("You requested the discretized completed matrix, but this step ", 
           "must be explicitly enabled in function soft_impute_tune(). Please ",
           "run soft_impute_tune() with 'discretize = TRUE' to perform this ", 
           "step.")
  }
  if (discretized) object$fit$X_discretized
  else object$fit$X
}

#' @rdname get_completed
#' @export
get_completed.soft_impute <- function(object, which, discretized = NULL, ...) {
  have_discretized <- !is.null(object$X_discretized)
  if (is.null(discretized)) discretized <- have_discretized
  else {
    discretized <- isTRUE(discretized)
    if (discretized & !have_discretized) 
      stop("You requested the discretized completed matrix, but this step ", 
           "must be explicitly enabled in function soft_impute(). Please ",
           "run soft_impute() with 'discretize = TRUE' to perform this step.")
  }
  if (discretized) out <- object$X_discretized
  else out <- object$X
  if (length(object$lambda) > 1L) {
    if (missing(which)) {
      stop("multiple regularization parameters used; use argument 'which' ", 
           "to specify the index for which to extract the completed matrix")
    } else if (!is.numeric(which) || length(which) != 1L) {
      stop("argument 'which' must be a single integer")
    }
    out <- out[[which]]
  }
  out
}

#' @rdname get_completed
#' @export
get_completed.median_impute <- function(object, discretized = NULL, ...) {
  have_discretized <- !is.null(object$X_discretized)
  if (is.null(discretized)) discretized <- have_discretized
  else {
    discretized <- isTRUE(discretized)
    if (discretized & !have_discretized) 
      stop("You requested the discretized completed matrix, but this step ", 
           "must be explicitly enabled in function median_impute(). Please ",
           "run median_impute() with 'discretize = TRUE' to perform this step.")
  }
  if (discretized) object$X_discretized
  else object$X
}

#' @rdname get_completed
#' @export
get_completed.mode_impute <- function(object, ...) {
  object$X
}


#' @rdname get_completed
#' @export
get_imputed <- get_completed


#' Extract the optimal value of the regularization parameter
#' 
#' Extract the optimal value of the regularization parameter from an object 
#' returned by a matrix completion algorithm.
#' 
#' @param object  an object returned by a matrix completion algorithm with a 
#' regularization parameter.
#' @param relative  logical; in case the values of the regularization parameter 
#' were given relative to a certain reference value computed from the data at 
#' hand, this allows to return the optimal value before or after multiplication 
#' with that reference value.
#' @param \dots  additional arguments to be passed down to methods.
#' 
#' @return  The optimal value of the regularization parameter.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{rdmc_tune}()}, \code{\link{soft_impute_tune}()}
#' 
#' @examples
#' # toy example derived from MovieLens 100K dataset
#' data("MovieLensToy")
#' # robust discrete matrix completion with hyperparameter tuning
#' set.seed(20250723)
#' fit <- rdmc_tune(MovieLensToy, 
#'                  lambda = fraction_grid(nb_lambda = 6),
#'                  splits = holdout_control(R = 5))
#' # extract optimal value of regularization parameter
#' get_lambda(fit)
#' 
#' # for more examples, see the help files of other functions for 
#' # matrix completion and imputation methods
#' 
#' @keywords utilities
#' 
#' @export

get_lambda <- function(object, ...) UseMethod("get_lambda")


#' @rdname get_lambda
#' @export
get_lambda.rdmc_tuned <- function(object, ...) {
  get_lambda(object$fit, ...)
}

#' @rdname get_lambda
#' @export
get_lambda.rdmc <- function(object, relative = TRUE, ...) {
  lambda <- object$lambda
  if (isTRUE(relative)) lambda 
  else lambda * object$d_max
}

#' @rdname get_lambda
#' @export
get_lambda.soft_impute_tuned <- function(object, ...) {
  get_lambda(object$fit, ...)
}

#' @rdname get_lambda
#' @export
get_lambda.soft_impute <- function(object, relative = TRUE, ...) {
  lambda <- object$lambda
  if (isTRUE(relative)) lambda 
  else lambda * object$lambda0
}


#' Extract the number of iterations
#' 
#' Extract the number of iterations from an object returned by a matrix 
#' completion algorithm.
#' 
#' @param object  an object returned by an iterative matrix completion 
#' algorithm.
#' @param \dots  currently ignored.
#' 
#' @return  The number of iterations performed in the iterative algorithm.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{rdmc_tune}()}
#' 
#' @examples
#' # toy example derived from MovieLens 100K dataset
#' data("MovieLensToy")
#' # robust discrete matrix completion with hyperparameter tuning
#' set.seed(20250723)
#' fit <- rdmc_tune(MovieLensToy, 
#'                  lambda = fraction_grid(nb_lambda = 6),
#'                  splits = holdout_control(R = 5))
#' # extract number of iterations with optimal regularization parameter
#' get_nb_iter(fit)
#' 
#' @keywords utilities
#' 
#' @export

get_nb_iter <- function(object, ...) UseMethod("get_nb_iter")


#' @rdname get_nb_iter
#' @export
get_nb_iter.rdmc_tuned <- function(object, ...) {
  object$fit$nb_iter
}

#' @rdname get_nb_iter
#' @export
get_nb_iter.rdmc <- function(object, ...) {
  object$nb_iter
}
