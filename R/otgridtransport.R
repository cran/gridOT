
otgridtransport <- \(x, y, p.1, p.2) {
	
	stopifnot(inherits(x, "otgrid") && inherits(y, "otgrid"))
	stopifnot(all.equal(x$total, y$total))
	stopifnot(is.numeric(p.1) && is.numeric(p.2) && isTRUE(p.1 >= 1) && isTRUE(p.2 >= 1))
	
	structure(
		list(
			from = x,
			to   = y,
			p.1  = p.1,
			p.2  = p.2
		), class = "otgridtransport"
	)
}

#' @export
print.otgridtransport <- \(x, ...) {
	
	cat("transport of", x$from$total, "mass between\n")
	cat(x$from$n, "x", x$from$m, "grid", setstr(x$from$x), "x", setstr(x$from$y))
	cat("\nand\n")
	cat(x$to$n, "x", x$to$m, "grid", setstr(x$to$x), "x", setstr(x$to$y))
	cat("\nwith respect to the cost |x_1 - y_1|^", x$p.1, " + |x_2 - y_2|^", x$p.2, "\n", sep = "")
	
	cat("with", x$to$n, "x", x$from$m, "pivot measure on", setstr(x$to$x), "x", setstr(x$from$y), "\n")
	
	if ("df" %in% names(x)) {
		cat("and transport plan of size", nrow(x$df), "\n")
	}
	if ("cost" %in% names(x)) {
		cat("and cost", x$cost, "\n")
	}
	
	invisible()
}

#' @param back.col color of the background.
#' @rdname plot-otgrid
#' @export
plot.otgridtransport <- \(x, num.col = 256, useRaster = TRUE, back.col = "lightblue", ...) {
	
	n.1 <- x$from$n
	n.2 <- x$from$m
	m.1 <- x$to$n
	m.2 <- x$to$m
	
	# rectangle for mu
	tmp <- unitrect.points(n.1, n.2, inner = FALSE)
	x.2 <- tmp$x	
	y.1 <- tmp$y

	# rectangle for xi
	tmp <- unitrect.points(m.1, m.2, inner = FALSE)
	x.1 <- tmp$x
	y.2 <- tmp$y + y.1[n.2 + 1]
	x.2 <- x.2 + x.1[m.1 + 1]

	# normalize
	tot <- x$from$total
	mu <- x$from$mass / tot
	nu <- x$to$mass / tot
	xi <- x$pivot / tot
	
	zlim <- range(mu, nu, xi)
	col <- seq(0, 1, len = num.col) |> grey()
	
	xm <- y.2[m.2 + 1]
	ym <- x.2[n.1 + 1]

	# x2 mu
	# x1 xi nu
	#    y1 y2

	plot.new()
	plot.window(xlim = c(0, xm), ylim = c(0, ym))
	
	rect(0, 0, xm, ym, col = back.col, border = NA)
	
	image(y.1, x.2, t(mu[n.1:1, , drop = FALSE]), zlim = zlim, col = col, add = TRUE, useRaster = useRaster) 
	image(y.1, x.1, t(xi[m.1:1, , drop = FALSE]), zlim = zlim, col = col, add = TRUE, useRaster = useRaster)
	image(y.2, x.1, t(nu[m.1:1, , drop = FALSE]), zlim = zlim, col = col, add = TRUE, useRaster = useRaster)
	
	invisible()
}
