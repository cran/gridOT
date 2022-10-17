
#' @title Pivot Measure
#' @description Calculate the pivot measure of the optimal transport between two-dimensional grids.
#' @param x an object of class `"otgrid"` the mass is to be transported from or 
#' a data frame with columns `from`, `to` and `mass` specifying the optimal transport plan.
#' @param y an object of class `"otgrid"` the mass is to be transported to.
#' @param p.1 the first power \eqn{\geq 1} of the cost.
#' @param p.2 the second power \eqn{\geq 1} of the cost.
#' @param ... further arguments (currently unused).
#' @return an object of class `"otgridtransport"` that contains the following elements:
#' \tabular{ll}{
#'   `from`  \tab an object of class `"otgrid"` the mass is transported from. \cr
#'   `to`    \tab an object of class `"otgrid"` the mass is transported to. \cr
#'   `p.1`   \tab the first power \eqn{\geq 1} of the cost function. \cr
#'   `p.2`   \tab the second power \eqn{\geq 1} of the cost function. \cr
#'   `pivot` \tab a matrix representing the pivot measure.  
#' }
#' If the Frank-Wolfe algorithm is used to approximate the pivot measure, the element `conv` indicates whether or not
#' we can ensured that the error is less or equal to the given precision `tol`.
#' 
#' If `return.it = TRUE`, then it also contains the vectors `costs` and `dualgaps` of costs and dual gaps in each iteration. 
#' 
#' Also note that the functions `print` and `plot` are available for objects of class `"otgridtransport"`. 
#' @details 
#' Denote with \eqn{X_1 \times X_2} and \eqn{Y_1 \times Y_2} the two-dimensional grids that `x` and `y` lie on, respectively. 
#' The pivot measure is a measure on \eqn{Y_1 \times X_2} that specifies the whole transport between `x` and `y` given that the cost
#' is of the separable form \eqn{c(x, y) = | x_1 - y_1 |^{p_1} + | x_2 - y_2 |^{p_2}}. 
#' 
#' The pivot measure is approximated using the Frank-Wolfe algorithm. The algorithm starts with an initial guess (`start.pivot`), e.g., the
#' independent coupling (`"independent"`) or the north-west-corner rule (`"northwestcorner"`). Then, in each iteration step, the dual solutions of
#' multiple one-dimensional transport problems are calculated and combined to give the gradient. To ensure differentiability, the dual solutions
#' must be unique. There are three different calculators for the dual solutions, of which the last two ensure uniqueness:
#' \tabular{lll}{
#' 	 name                  \tab description \tab `dual.params` \cr
#'   `"discrete"`          \tab the distributions are basically unchanged, due to rounding errors the support points
#'                            are moved by a small positive amount.  \tab `right.margin = 1e-15` \cr
#'   `"epsilon-discrete"`  \tab additionally, the support is made connected by uniformly adding \eqn{\varepsilon}
#'                            mass. \tab + `eps = 1e-8` \cr
#'   `"epsilon-histogram"` \tab additionally, each point mass is distributed uniformly in a bin around said point \tab + `width = 1e-8`
#' }
#' 
#' A pivot measure can also be computed from an optimal transport plan.
#' 
#' Finally, the pivot measure can be used to calculate the full transport plan and cost. 
#' @examples 
#' x <- otgrid(rbind(c(0, 0, 0, 0, 0, 0, 0, 0, 0),
#'	                 c(0, 0, 1, 0, 0, 0, 0, 0, 0),
#'	                 c(0, 1, 2, 1, 0, 0, 0, 0, 0),
#'	                 c(0, 0, 1, 0, 0, 0, 0, 0, 0), 
#'                   0, 0, 0, 0, 0))
#' y <- otgrid(rbind(0, 0, 0, 0,
#'	                 c(0, 0, 0, 0, 0, 0, 0, 0, 0),
#'	                 c(0, 0, 0, 0, 0, 0, 1, 0, 0),
#'	                 c(0, 0, 0, 0, 0, 1, 2, 1, 0),
#'	                 c(0, 0, 0, 0, 0, 0, 1, 0, 0), 0))
#' 
#' pm <- pivot_measure(x, y, p.1 = 1, p.2 = 2, dual.method = "epsilon-discrete")
#' 
#' print(pm)
#' plot(pm)
#' 
#' # use pivot measure to calculate cost and plan
#' pm <- transport_cost(pm)
#' pm <- transport_df(pm)
#' print(pm)
#' 
#' # calculate pivot from plan
#' pm2 <- pivot_measure(pm$df, x, y)
#' plot(pm2)
#' @seealso transport plan [`transport_df`], transport cost [`transport_cost`], two-dimensional grid [`otgrid`], 
#' plot [`plot.otgridtransport`]
#' @export 
pivot_measure <- \(x, ...) {
	UseMethod("pivot_measure")
}

start_pivot <- \(x, y, method = c("independent", "northwestcorner")) {
	
	mu.2 <- x$marg.y
	nu.1 <- y$marg.x
	
	if (is.character(method)) {
		method <- match.arg(method)
		
		if (method == "independent") {
			tcrossprod(nu.1, mu.2)
		}
		else {
			northWestCorner(nu.1, mu.2, 1e-15)
		}
	}
	else {
		if (is.matrix(method)) {
			xi <- method
		}
		else {
			xi <- method(nu.1, mu.2)
		}
		
		stopifnot(is.matrix(xi), nrow(xi) == y$n && ncol(xi) == x$m,
		          all.equal(nu.1, rowSums(xi)), all.equal(mu.2, colSums(xi)))
		
		xi
	}
}

dual_params <- \(x, y, method, params) {
	stopifnot("dual.params is not a list" = is.list(params))
	
	if (is.null(params$right.margin)) {
		params$right.margin <- 1e-15
	}
	else {
		stopifnot(params$right.margin >= 0)
	}
	
	if (method == "epsilon-histogram") {
		
		min.width <- min(diff(x$x), diff(y$y))
		
		if (is.null(params$width)) {
			params$width <- min(min.width, 1e-8)
		}
		else {
			stopifnot(min.width >= params$width)
		}
	}
	
	if (method != "discrete") {
		if (is.null(params$eps)) {
			params$eps <- 1e-8
		}
		else {
			stopifnot(0 < params$eps && params$eps < 1)
		}
	}
	
	params
}

#' @param max.it the maximum number of iterations.
#' @param tol the desired accuracy.
#' @param threads number of threads to use.
#' @param start.pivot the start pivot to use: matrix representing an actual coupling or `"independent"` or `"northwestcorner"`.
#' @param dual.method the name of the dual calculators to use: `"discrete"`, `"epsilon-discrete"` or `"epsilon-histogram"`.
#' @param dual.params list of parameters for the dual calculators.
#' @param return.it logical value specifying whether or not the costs and dual gaps in each iteration are to be returned .
#' @rdname pivot_measure
#' @method pivot_measure otgrid
#' @export
pivot_measure.otgrid <- \(x, y, p.1 = 2, p.2 = p.1, max.it = 100, tol = 1e-4, threads = 1, 
                          start.pivot = c("independent", "northwestcorner"), 
                          dual.method = c("discrete", "epsilon-discrete", "epsilon-histogram"),
                          dual.params = list(),
                          return.it = FALSE, ...) {
	
	dual.method <- match.arg(dual.method)
	
	res <- otgridtransport(x, y, p.1, p.2)
	
	x.0 <- without_zeros(x)
	y.0 <- without_zeros(y)

	start.pivot <- start_pivot(x.0, y.0, start.pivot)
	
	dual.params <- dual_params(x.0, y.0, dual.method, dual.params)

	fw <- switch(
		dual.method,
		"discrete" = {
			frankWolfeDiscrete(x.0$x, x.0$y, x.0$mass, y.0$x, y.0$y, y.0$mass, 
			                   p.1, p.2, start.pivot, 
			                   max.it, tol, threads, return.it, 
			                   dual.params$right.margin)
		},
		"epsilon-discrete" = {
			frankWolfeEpsilonDiscrete(x.0$x, x.0$y, x.0$mass, y.0$x, y.0$y, y.0$mass, 
			                          p.1, p.2, start.pivot, 
			                          max.it, tol, threads, return.it, 
			                          dual.params$right.margin, dual.params$eps)
		},
		"epsilon-histogram" =  {
			frankWolfeEpsilonHistogram(x.0$x, x.0$y, x.0$mass, y.0$x, y.0$y, y.0$mass, 
			                           p.1, p.2, start.pivot, 
			                           max.it, tol, threads, return.it, 
			                           dual.params$right.margin, dual.params$eps, dual.params$width)
		}
	)
	
	res$pivot <- emb_mat(y$n, x$m, fw$pivot, y.0$ind.x, x.0$ind.y)
	res$conv <- fw$conv
	
	if (return.it) {
		res$costs <- fw$costs
		res$dualgaps <- fw$dualgaps
	}
	
	res
}

#' @param a an object of class `"otgrid"` the mass is to be transported from.
#' @param b an object of class `"otgrid"` the mass is to be transported to.
#' @rdname pivot_measure
#' @method pivot_measure data.frame
#' @export
pivot_measure.data.frame <- \(x, a, b, p.1 = 2, p.2 = p.1, ...) {

	stopifnot(c("from", "to", "mass") %in% colnames(x))
	
	res <- otgridtransport(a, b, p.1, p.2)

	res$pivot <- pivotMeasure(x$from, x$to, x$mass, a$n, a$m, b$n)

	res
}
