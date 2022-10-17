#include<RcppArmadillo.h>
#include"transport_cost_1d.h"
#include"transport_cost.h"

// [[Rcpp::export]]
double transportCost(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, 
                     const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, 
                     const double p1, const double p2, const arma::mat& xi,
                     const double threshold)
{
	const int m1 = nu.n_rows;
	const int m2 = nu.n_cols;
	const int n1 = mu.n_rows;
	const int n2 = mu.n_cols;

	double cost = 0;

	for (int r = 0; r < m1; ++r)
	{
		cost += transportCost1d(y2, nu.row(r), m2, x2, xi.row(r), n2, p2, threshold);
	}
	
	for (int j = 0; j < n2; ++j)
	{
		cost += transportCost1d(x1, mu.col(j), n1, y1, xi.col(j), m1, p1, threshold);
	}

	return cost;
}

// [[Rcpp::export]]
double transportCostFromPlan(const Rcpp::IntegerVector& from, const Rcpp::IntegerVector& to,
                             const Rcpp::NumericVector& mass, const Rcpp::NumericMatrix& cost)
{
	double s = 0;

	const int n = from.size();
	for (int i = 0; i < n; ++i)
	{
		s += mass[i] * cost(from[i] - 1, to[i] - 1);
	}

	return s;
}
