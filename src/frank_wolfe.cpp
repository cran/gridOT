#include<RcppArmadillo.h>
#include"frank_wolfe_generic.h"
#include"construct_duals.h"

// [[Rcpp::export]]
Rcpp::List frankWolfeDiscrete(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, 
                              const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, 
                              const double p1, const double p2, const arma::mat& startPivot,
                              const int maxIt, const double tol, const int threads, const bool returnIt,
                              const double rightMargin)
{
	return frankWolfeGeneric(x1, x2, mu, y1, y2, nu, p1, p2, startPivot, maxIt, tol, threads, returnIt,
	                         constructDiscrete, constructDiscrete, rightMargin);
}

// [[Rcpp::export]]
Rcpp::List frankWolfeEpsilonDiscrete(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, 
                                     const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, 
                                     const double p1, const double p2, const arma::mat& startPivot,
                                     const int maxIt, const double tol, const int threads, const bool returnIt,
                                     const double rightMargin, const double eps)
{
	return frankWolfeGeneric(x1, x2, mu, y1, y2, nu, p1, p2, startPivot, maxIt, tol, threads, returnIt,
	                         constructEpsilonDiscrete, constructEpsilonDiscrete, rightMargin, eps);
}

// [[Rcpp::export]]
Rcpp::List frankWolfeEpsilonHistogram(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, 
                                      const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, 
                                      const double p1, const double p2, const arma::mat& startPivot,
                                      const int maxIt, const double tol, const int threads, const bool returnIt,
                                      const double rightMargin, const double eps, const double width)
{
	return frankWolfeGeneric(x1, x2, mu, y1, y2, nu, p1, p2, startPivot, maxIt, tol, threads, returnIt,
	                         constructEpsilonHistogram, constructEpsilonHistogram, rightMargin, eps, width);
}
