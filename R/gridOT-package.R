
#' @description \if{html}{\figure{logo.png}{options: style='float: right' alt='logo' width='120'}}
#' Can be used for optimal transport between two-dimensional grids with respect to separable cost functions of l^p form. It utilizes the 
#' Frank-Wolfe algorithm to approximate so-called pivot measures: one-dimensional transport plans that fully describe the full transport, see G.
#' Auricchio (2021) \href{https://arxiv.org/abs/2105.07278}{arXiv:2105.07278}. For these, it offers methods for visualization and to extract the corresponding transport plans and costs. 
#' Additionally, related functions for one-dimensional optimal transport are available.
#' @references G. Auricchio (2021). On the Pythagorean Structure of the Optimal Transport for Separable Cost Functions. arXiv preprint \href{https://arxiv.org/abs/2105.07278}{arXiv:2105.07278}.
#' @keywords internal
"_PACKAGE"

#' @importFrom Rcpp evalCpp
#' @importFrom grDevices grey
#' @importFrom graphics image plot.new plot.window rect
#' @useDynLib gridOT, .registration = TRUE
NULL
