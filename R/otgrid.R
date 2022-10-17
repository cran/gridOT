
unitrect.points <- \(n, m, inner = TRUE) {
	
	h <- 1 / n
	
	if (inner) {
		s <- h / 2
	}
	else {
		s <- 0
		n <- n + 1
		m <- m + 1
	}

	list(x = seq(s, by = h, len = n), y = seq(s, by = h, len = m))
}

#' @title Two-dimensional Grid with Mass
#' @description Create an object that represents a probability measure that is supported on a two-dimensional grid.
#' @param mass matrix of non-negative weights.
#' @param x vector of support points of the first marginal.
#' @param y vector of support points of the second marginal.
#' @param sorted logical value specifying whether or not the support points are sorted.
#' @param normalize logical value specifying whether or not the total mass should be rescaled to 1.
#' @param remove.zeros logical value specifying whether or not marginals with no mass should be removed from the grid.
#' @return object of class `"otgrid"`. It contains the following elements:
#' \tabular{ll}{
#'   `x`     \tab vector of support points of the first marginal  \cr
#'   `y`     \tab vector of support points of the second marginal \cr
#'   `n`     \tab number of support points of the first marginal  \cr
#'   `m`     \tab number of support points of the second marginal \cr
#'   `mass`  \tab matrix of non-negative weights                  \cr
#'   `total` \tab total mass
#' }
#' 
#' Also note that the functions `print` and `plot` are available for objects of class `"otgrid"`.
#' @details 
#' If `x` and `y` are not specified, then a equispaced unit grid is chosen.
#' @seealso plot [`plot.otgrid`]
#' @examples
#' x <- otgrid(cbind(1:2, 0, 3:4), remove.zeros = TRUE)
#' 
#' print(x) # note that it's only 2 x 2 
#' plot(x)
#' @export
otgrid <- \(mass, x = NULL, y = NULL, sorted = FALSE, normalize = FALSE, remove.zeros = FALSE) {

	stopifnot(is.null(x) == is.null(y))

	sorted <- sorted || (is.null(y) && is.null(x))

	if (is.null(x)) {
		tmp <- unitrect.points(nrow(mass), ncol(mass))
		x <- tmp$x
		y <- tmp$y
	}
	
	stopifnot(is.numeric(x) && is.numeric(y) && is.numeric(mass))
	stopifnot("dimensions do not match" = nrow(mass) == length(x) && ncol(mass) == length(y),
	          "negative mass"           = all(mass >= 0))
	stopifnot(anyDuplicated(x) == 0 && anyDuplicated(y) == 0)

	if (!sorted) {
		i <- order(x)
		j <- order(y)
		
		x <- x[i]
		y <- y[j]
		
		mass <- mass[i, j, drop = FALSE]
	}
	
	if (remove.zeros) {
		a <- rowSums(mass) > 0
		b <- colSums(mass) > 0

		x <- x[a]
		y <- y[b]

		mass <- mass[a, b, drop = FALSE]
	}
	
	n <- length(x)
	m <- length(y)
	
	stopifnot("empty grid" = n > 0 && m > 0)
	
	total <- sum(mass)

	stopifnot(total > 0)

	if (normalize) {
		mass <- mass / total
		total <- 1
	}

	structure(
		list(
			x     = x,
			y     = y,
			n     = n,
			m     = m,
			mass  = mass,
			total = total
		), class = "otgrid"
	)
}

setstr <- \(z) {
	n <- length(z)
	if (n == 1) {
		sprintf("{%f}", z)
	}
	else if (n == 2) {
		sprintf("{%f,%f}", z[1], z[2])
	}
	else {
		sprintf("{%f,...,%f}", z[1], z[n])
	}
}

#' @export
print.otgrid <- \(x, ...) {
	
	cat(sprintf("%d x %d grid on \n", x$n, x$m))
	
	cat(setstr(x$x))
	cat(" x ")
	cat(setstr(x$y), "\n")
	
	cat("with total mass", x$total)
	
	invisible()
}

#' @title Plots for two-dimensional Optimal Transport
#' @description Plot two-dimensional grids or visualize the pivot measure.
#' @details 
#' For objects of class `"otgrid"`, the grid is plotted as a greyscale image.
#' 
#' For objects of class `"otgridtransport"`, the pivot measure is visualized as follows: the two grids of the transport are plotted 
#' in the left upper (`from`) and right lower (`to`) corner. The pivot measure is in the left lower corner such that the marginals match. 
#' @param x an object of class `"otgrid"` or `"otgridtransport"`
#' @param num.col number of colors (grey) to use.
#' @param useRaster parameter passed to the [`image`] function.
#' @param ... further arguments (currently unused).
#' @return No return value, called for side effects.
#' @seealso pivot measure [`pivot_measure`], grid [`otgrid`]
#' @rdname plot-otgrid
#' @export
plot.otgrid <- \(x, num.col = 256, useRaster = TRUE, ...) {
	
	n <- x$n
	m <- x$m
	
	a <- unitrect.points(n, m, inner = FALSE)

	mass <- (x$mass / x$total)[n:1, , drop = FALSE] |> t()
	
	zlim <- range(mass)
	col <- seq(0, 1, len = num.col) |> grey()
	
	plot.new()
	plot.window(xlim = c(0, a$y[m + 1]), ylim = c(0, a$x[n + 1]))
	
	image(a$y, a$x, mass, zlim = zlim, col = col, add = TRUE, useRaster = useRaster) 
	
	invisible()
}
