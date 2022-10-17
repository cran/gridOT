
emb_vec <- \(n, vec, i.x, x) {
	
	if (n == length(vec)) {
		return(vec)
	}
	
	y <- numeric(n)
	y[-i.x] <- vec
	y[i.x] <- x	
	
	y
}

emb_mat <- \(n, m, mat, i.x, i.y) {
	
	n.0 <- nrow(mat)
	m.0 <- ncol(mat)
	
	if (n.0 == n && m.0 == m) {
		return(mat)
	}
	
	z <- matrix(0, n, m)
	
	if (n.0 == n && m.0 != m) {
		z[, -i.y] <- mat
	}
	else if (n.0 != n && m.0 == m) {
		z[-i.x, ] <- mat
	}
	else {
		z[-i.x, -i.y] <- mat		
	}
	
	z
}

# remove points where marginal weight is zero
without_zeros <- \(x) {

	marg.x <- rowSums(x$mass)
	marg.y <- colSums(x$mass) 

	mask.x <- marg.x > 0
	mask.y <- marg.y > 0

	x.0 <- x$x[mask.x]
	y.0 <- x$y[mask.y]

	n.0 <- length(x.0)
	m.0 <- length(y.0)
	
	stopifnot(n.0 > 0 && m.0 > 0)
	
	list(
		x      = x.0,
		y      = y.0,
		n      = n.0,
		m      = m.0,
		ind.x  = which(!mask.x),
		ind.y  = which(!mask.y),
		mass   = x$mass[mask.x, mask.y, drop = FALSE],
		marg.x = marg.x[mask.x],
		marg.y = marg.y[mask.y]
	)
}
