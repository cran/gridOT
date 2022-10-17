#include<Rcpp.h>
#include<algorithm>

// [[Rcpp::export]]
Rcpp::NumericMatrix northWestCorner(const Rcpp::NumericVector& mu, const Rcpp::NumericVector& nu, const double threshold)
{
	const int n = mu.size();
	const int m = nu.size();

	Rcpp::NumericMatrix xi(n, m);

	int i = 0;
	int j = 0;

	double r = mu[i];
	double s = nu[j];
	double t = 0;

	while (true)
	{
		while (r <= threshold)
		{
			if (++i < n)
			{
				r = mu[i];
			}
			else
			{
				return xi;
			}
		}

		while (s <= threshold) 
		{
			if (++j < m)
			{
				s = nu[j];
			}
			else
			{
				return xi;
			}
		}

		t = std::min(r, s);
		
		xi(i, j) = t;

		r -= t;
		s -= t;
	}
}
