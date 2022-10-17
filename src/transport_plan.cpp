#include<RcppArmadillo.h>
#include<vector>
#include"transport_plan_1d.h"

// [[Rcpp::export]]
Rcpp::DataFrame transportPlan(const arma::mat& mu, const arma::mat& nu, const arma::mat& xi, const double threshold = 1e-15)
{
	const int m1 = nu.n_rows;
	const int m2 = nu.n_cols;
	const int n1 = mu.n_rows;
	const int n2 = mu.n_cols;
	
	std::vector<int> from;
	std::vector<int> to;
	std::vector<double> mass;
	
	// reserve how much?
	// from.reserve();
	// to.reserve();
	// mass.reserve();

	std::vector<std::vector<Transport>> plans2;
	plans2.reserve(n2);

	// precalculate the plans for the columns
	for (int j = 0; j < n2; ++j)
	{
		plans2.push_back(transportPlan1d(xi.col(j), m1, mu.col(j), n1, threshold));
	}

	for (int r = 0; r < m1; ++r)
	{
		std::vector<Transport> plan1 = transportPlan1d(xi.row(r), n2, nu.row(r), m2, threshold);
		const int n = plan1.size();
		
		for (int a = 0; a < n; ++a)
		{
			Transport* transp1 = &plan1[a];

			int j = transp1->from;
			int s = transp1->to;

			if (xi(r, j) <= 0)
			{	
				continue;
			}

			std::vector<Transport>* plan2 = &plans2[j];
			const int m = plan2->size();

			Transport* transp2;

			// find entry point for r
			int b = 0;
			while (b < m && (transp2 = &(*plan2)[b])->from < r)
			{
				++b;
			}

			// as long as we have transports from r
			while (b < m && (transp2 = &(*plan2)[b])->from == r)
			{
				int i = transp2->to;

				from.push_back(i + j * n1 + 1);
				to.push_back(r + s * m1 + 1);

				mass.push_back(transp1->mass * transp2->mass / xi(r, j));

				++b;
			}
		}
	}
	
    return Rcpp::DataFrame::create(Rcpp::Named("from") = from, Rcpp::Named("to") = to, Rcpp::Named("mass") = mass);
}
