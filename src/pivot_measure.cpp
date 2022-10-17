#include<Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix pivotMeasure(const Rcpp::IntegerVector& from, const Rcpp::IntegerVector& to, 
                                 const Rcpp::NumericVector& mass, const int n1, const int n2, const int m1)
{
	Rcpp::NumericMatrix xi(m1, n2);

	const int N = from.size();

	for (int k = 0; k < N; ++k)
	{
		int j = (from[k] - 1) / n1;
		int r = to[k] - ((to[k] - 1) / m1) * m1 - 1;

		xi(r, j) += mass[k];
	}

	return xi;
}
