#include<Rcpp.h>
#include<cmath>

// [[Rcpp::export]]
Rcpp::NumericMatrix costMatrix(const Rcpp::NumericVector& x1, const Rcpp::NumericVector& x2,
                               const Rcpp::NumericVector& y1, const Rcpp::NumericVector& y2, const double p1, const double p2)
{
	const int n1 = x1.size();
	const int n2 = x2.size();
	const int m1 = y1.size();
	const int m2 = y2.size();

	Rcpp::NumericMatrix cost(n1 * n2, m1 * m2);
	
	for (int i = 0; i < n1; ++i)
	{
		for (int r = 0; r < m1; ++r) 
		{
			const double dist1 = std::pow(std::fabs(x1[i] - y1[r]), p1);

			for (int j = 0; j < n2; ++j) 
			{
				for (int s = 0; s < m2; ++s) 
				{
					cost(i + j * n1, r + s * m1) = dist1 + std::pow(std::fabs(x2[j] - y2[s]), p2);
				}
			}
		}
	}
	
	return cost;
}
