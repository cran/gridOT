
library(hexSticker)
library(tidyverse)
library(gggenes)
library(showtext)

#View(font_files() -> tmp.fonts)
font_add("Consolas", "consola.ttf") 
showtext_auto()

border.col <- "#004182" 
fill.col <- "#118DF0"
line.size <- 0.5

plot_grid <- \(x, y, mu, col) {
	
	xr <- range(x)
	yr <- range(y)
	
	h <- lapply(x, \(x) geom_line(aes(x = x, y = yr), size = line.size, col = border.col))
	v <- lapply(y, \(y) geom_line(aes(x = xr, y = y), size = line.size, col = border.col))
	
	df <- expand.grid(x = x, y = y) |>
		as_tibble() |>
		add_column(z = c(mu)) |>
		filter(z > 0)
	
	p <- geom_point(aes(x = x, y = y, size = z), df, col = col)
	
	c(h, v, p)
}

normalize <- \(mu) mu / sum(mu)

some_zero <- \(mu, k) {
	n <- nrow(mu)
	m <- ncol(mu)
	
	i <- sample(1:n, k)
	j <- sample(1:m, k)
	
	for (l in 1:k) {
		mu[i[l], j[l]] <- 0
	}
	
	mu
}

rmass <- \(n, m, k) runif(n * m) |> matrix(n, m) |> some_zero(k) |> normalize()

bx <- sqrt(3) / 2

set.seed(2)

ul <- list(
	x   = 1 - bx / 2 + c(-0.25, -0.05, 0.05, 0.2),
	y   = c(0.6, 0.75, 0.85, 1),
	mu  = rmass(4, 4, 2),
	col = "#FFC045"
)

dr <- list(
	x   = 1 + bx / 2 + c(-0.2, -0.12, 0.02, 0.1, 0.25),
	y   = c(0.5, 0.65, 0.85, 0.95, 1.1),
	mu  =  rmass(5, 5, 3),
	col = "#BBE1FA"
)

ggplot() + theme_sticker(legend.position = "none") + scale_size_continuous(range = c(0, 1.75)) +
	geom_hexagon(size = 1.25, color = border.col, fill = fill.col) +
	geom_pkgname("gridOT", x = 1, y = 1.4, family = "Consolas", 
				 size = 30, color = "#142850") +
	do.call(plot_grid, ul) + do.call(plot_grid, dr) +
	geom_gene_arrow(aes(xmin = 0.9, xmax = 1.1, y = 0.8), size = line.size, col = border.col,
	                arrowhead_height = unit(2, "mm"), fill = fill.col,
					arrow_body_height = unit(1, "mm"),
					arrowhead_width = unit(1.5, "mm"))

save_sticker("logo.png")
