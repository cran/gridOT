
#' @title Cost Matrix for two-dimensional Optimal Transport
#' @description Calculate the cost matrix for the optimal transport between two-dimensional grids with respect to the cost function 
#' \deqn{c(x, y) = | x_1 - y_1 |^{p_1} + | x_2 - y_2 |^{p_2}.}
#' @param ... further arguments (currently unused).
#' @return a matrix giving the pairwise costs in column-mayor format. 
#' @examples
#' x <- otgrid(cbind(0:1, 1:0))
#' y <- otgrid(cbind(1:0, 0:1))
#' 
#' transport_costmat(x, y, p.1 = 1, p.2 = 3)
#' @export
transport_costmat <- \(x, ...) {
	UseMethod("transport_costmat")
}

#' @param x an object of class `"otgridtransport"` or `"otgrid"`, in the latter case it gives the object the mass is to be transported from.
#' @param y an object of class `"otgrid"` the mass is to be transported from.
#' @param p.1 the first power \eqn{\geq 1} of the cost.
#' @param p.2 the second power \eqn{\geq 1} of the cost.
#' @rdname transport_costmat
#' @method transport_costmat otgrid
#' @export
transport_costmat.otgrid <- \(x, y, p.1 = 2, p.2 = p.1, ...) {

	stopifnot(inherits(y, "otgrid"))
	stopifnot(is.numeric(p.1) && is.numeric(p.2) && isTRUE(p.1 >= 1) && isTRUE(p.2 >= 1))
	
	costMatrix(x$x, x$y, y$x, y$x, p.1, p.2)
}

#' @rdname transport_costmat
#' @method transport_costmat otgridtransport
#' @export
transport_costmat.otgridtransport <- \(x, ...) {
	transport_costmat(x$from, x$to, p.1 = x$p.1, p.2 = x$p.2)
}
