
check.weights <- \(wx, wy) {
	stopifnot(length(wx) > 0 && length(wy) > 0)
	stopifnot(all.equal(sum(wx), sum(wy)))
	stopifnot(all(wx >= 0) && all(wy >= 0))
}

#' @title North-west-corner Rule
#' @description Calculate the transport plan obtained by the north-west-corner rule.
#' @param wx first weight vector.
#' @param wy second weight vector.
#' @param threshold small value that indicates when a value is considered to be zero.
#' @return a matrix representing the transport plan obtained by the north-west-corner rule.
#' @examples 
#' set.seed(1)
#' wx <- rep(1/5, 5)
#' wy <- runif(6)
#' wy <- wy / sum(wy)
#' north_west_corner(wx, wy)
#' @export
north_west_corner <- \(wx, wy, threshold = 1e-15) {
	
	check.weights(wx, wy)
	
	northWestCorner(wx, wy, threshold = threshold)
}

#' @param x first vector of points.
#' @param y second vector of points.
#' @param wx weight vector of the first vector of points.
#' @param wy weight vector of the second vector of points.
#' @param p the power \eqn{\geq 1} of the cost function.
#' @param sorted logical value indicating whether or not `a` and `b` are sorted.
#' @rdname transport_cost
#' @method transport_cost numeric
#' @export
transport_cost.numeric <- \(x, y, wx, wy, p = 1, sorted = FALSE, threshold = 1e-15, ...) {
	
	n <- length(x)
	m <- length(y)
	
	stopifnot(p >= 1)
	stopifnot(n > 0 && n == length(wx))
	stopifnot(m > 0 && m == length(wy))
	check.weights(wx, wy)

	if (!sorted) {
		i <- order(x)
		j <- order(y)
		
		x <- x[i]
		y <- y[j]
		
		wx <- wx[i]
		wy <- wy[j]
	}
	
	transportCost1d(x, y, wx, wy, p, threshold = threshold)
}

#' @param x first weight vector.
#' @param y second weight vector.
#' @param threshold small value that indicates when a value is considered to be zero.
#' @rdname transport_df
#' @method transport_df numeric
#' @export
transport_df.numeric <- \(x, y, threshold = 1e-15, ...) {
	
	check.weights(x, y)
	
	transportPlan1d(x, y, threshold = threshold)
}

#' @title Dual Solution of one-dimensional Optimal Transport
#' @description Calculate an optimal dual pair for the optimal transport between discrete one-dimensional measures.
#' @details 
#' The pair \eqn{f, g} is an optimal dual pair if the optimal transport distance between the two distributions
#' with respect to the cost function \eqn{c(x, y) = | x - y |^p} is given by 
#' \deqn{\langle f, w_a \rangle + \langle g, w_b \rangle}
#' and the condition \eqn{f_i + g_j \leq | a_i - b_j |^p} holds.
#' @param a first vector of points.
#' @param b second vector of points.
#' @param wa weight vector of the first vector of points.
#' @param wb weight vector of the second vector of points.
#' @param p the power \eqn{\geq 1} of the cost function.
#' @param right.margin small amount the points are moved by.
#' @param sorted logical value indicating whether or not `a` and `b` are sorted.
#' @return a list containing the dual vectors `pot.a` and `pot.b`.
#' @examples 
#' set.seed(1)
#' a <- 1:5
#' wa <- rep(1/5, 5)
#' b <- 1:6
#' wb <- runif(6)
#' wb <- wb / sum(wb)
#' 
#' d <- dual1d(a, b, wa, wb, p = 1)
#' 
#' dc <- sum(d$pot.a * wa) + sum(d$pot.b * wb)
#' print(all.equal(dc, transport_cost(a, b, wa, wb, p = 1)))
#' @export
dual1d <- \(a, b, wa, wb, p = 1, right.margin = 1e-15, sorted = FALSE) {

	n <- length(a)
	m <- length(b)

	stopifnot(p >= 1)
	stopifnot(n > 0 && n == length(wa))
	stopifnot(m > 0 && m == length(wb))
	check.weights(wa, wb)

	if (!sorted) {
		ord.a <- order(a)
		ord.b <- order(b)
		
		a <- a[ord.a]
		wa <- wa[ord.a]

		b <- b[ord.b]
		wb <- wb[ord.b]
	}

	c.a <- cumsum(wa)
	c.b <- cumsum(wb)

	q.a <- b[findInterval(c.a[-n], c.b + right.margin, left.open = TRUE) + 1]
	pot.a <- c(0, cumsum(abs(q.a - a[-1])^p - abs(q.a - a[-n])^p))

	i <- findInterval(b, q.a, left.open = TRUE) + 1
	pot.b <- abs(a[i] - b)^p - pot.a[i] 

	if (!sorted) {
		pot.a <- pot.a[order(ord.a)]
		pot.b <- pot.b[order(ord.b)]
	}

	list(pot.a = pot.a, pot.b = pot.b)
}
