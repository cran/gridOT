

#' @title Optimal Transport Cost
#' @description Calculate the optimal transport cost.
#' @details 
#' In case of two-dimensional grids, the pivot measure is used to calculate the optimal transport cost.
#' 
#' For one-dimensional optimal transport, the cost function is given by \eqn{c(x, y) = | x - y |^p}. In this case, the north-west-corner
#' algorithm is used. 
#' @return the optimal transport cost or, in case of two-dimensional case, an object of class `"otgridtransport"` with element `cost` that contains it.
#' @examples
#' ## one-dimensional example
#' set.seed(1)
#' a <- 1:5
#' wa <- rep(1/5, 5)
#' b <- 1:6
#' wb <- runif(6)
#' wb <- wb / sum(wb)
#' transport_cost(a, b, wa, wb, p = 1)
#' 
#' ## two-dimensional example
#' x <- otgrid(cbind(0:1, 1:0))
#' y <- otgrid(cbind(1:0, 0:1))
#' 
#' # first calculate pivot manually
#' pm <- pivot_measure(x, y)
#' pm <- transport_cost(pm)
#' print(pm$cost)
#' 
#' # or just
#' pm2 <- transport_cost(x, y)
#' print(pm2$cost)
#' 
#' # or from transport plan and cost matrix
#' costm <- transport_costmat(pm)
#' tp <- transport_df(pm)
#' print(transport_cost(tp$df, costm))
#' @seealso pivot measure [`pivot_measure`]
#' @export
transport_cost <- \(x, ...) {
	UseMethod("transport_cost")
}

#' @param x a vector of points; a data frame with columns `from`, `to` and `mass` specifying the optimal transport plan; 
#' an object of class `"otgridtransport"` or `"otgrid"`, in the latter case `...` must be 
#' the arguments of [`pivot_measure`].
#' @param threshold small value that indicates when a value is considered to be zero.
#' @rdname transport_cost
#' @method transport_cost otgridtransport
#' @export
transport_cost.otgridtransport <- \(x, threshold = 1e-15, ...) {
	x$cost <- transportCost(x$from$x, x$from$y, x$from$mass, x$to$x, x$to$y, x$to$mass, x$p.1, x$p.2, x$pivot, threshold)
	x
}

#' @param ... further arguments (for [`pivot_measure`] if `x` is an object of class `"otgrid"`).
#' @rdname transport_cost
#' @method transport_cost otgrid
#' @export
transport_cost.otgrid <- \(x, ...) {
	transport_cost(pivot_measure(x, ...))
}

#' @param costm cost matrix of the transport
#' @rdname transport_cost
#' @method transport_cost data.frame
#' @export
transport_cost.data.frame <- \(x, costm, ...) {
	transportCostFromPlan(x$from, x$to, x$mass, costm)
}
