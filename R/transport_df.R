
#' @title Optimal Transport Plan
#' @description Calculate the optimal transport plan.
#' @details 
#' In case of two-dimensional grids, the pivot measure is used to calculate the optimal transport plan.
#' 
#' For one-dimensional optimal transport, we assume that the optimal transport plan is the monotone plan. For example, this is the case
#' for costs of the form \eqn{c(x, y) = | x - y |^p} for some \eqn{p \geq 1}.  
#' In this case, the north-west-corner algorithm is used. 
#' @return 
#' a data frame representing the optimal transport plan. It has columns `from`, `to` and `mass` that specify
#' between which points of the two grids how much mass is transported. 
#' 
#' In the two-dimensional case, a point is given by the index in column-mayor format and the data frame is actually stored in the element
#' `df` of an object of class `"otgridtransport"`.
#' @examples
#' ## one-dimensional example
#' set.seed(1)
#' wa <- rep(1/5, 5)
#' wb <- runif(6)
#' wb <- wb / sum(wb)
#' transport_df(wa, wb)
#' 
#' ## two-dimensional example
#' x <- otgrid(cbind(0:1, 1:0))
#' y <- otgrid(cbind(1:0, 0:1))
#' 
#' # first calculate pivot manually
#' pm <- pivot_measure(x, y)
#' pm <- transport_df(pm)
#' 
#' # or just
#' pm2 <- transport_df(x, y)
#' @seealso pivot measure [`pivot_measure`]
#' @export
transport_df <- \(x, ...) {
	UseMethod("transport_df")
}

#' @param x a vector of weights, an object of class `"otgridtransport"` or `"otgrid"`, in the latter case `...` must be the arguments of [`pivot_measure`].
#' @rdname transport_df
#' @method transport_df otgridtransport
#' @export
transport_df.otgridtransport <- \(x, ...) {
	x$df <- transportPlan(x$from$mass, x$to$mass, x$pivot)
	x
}

#' @param ... further arguments (for [`pivot_measure`] if `x` is an object of class `"otgrid"`).
#' @rdname transport_df
#' @method transport_df otgrid
#' @export
transport_df.otgrid <- \(x, ...) {
	transport_df(pivot_measure(x, ...))
}
