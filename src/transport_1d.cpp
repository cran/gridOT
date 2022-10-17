#include<Rcpp.h>
#include"transport_cost_1d.h"
#include"transport_plan_1d.h"

// [[Rcpp::export(transportPlan1d)]]
Rcpp::DataFrame RcppTransportPlan1d(const Rcpp::NumericVector& wa, const Rcpp::NumericVector& wb, const double threshold = 1e-15)
{
	std::vector<Transport> plan = transportPlan1d(wa, wa.size(), wb, wb.size(), threshold);

	const int n = plan.size();

	Rcpp::IntegerVector from(n);
	Rcpp::IntegerVector to(n);
	Rcpp::NumericVector mass(n);

	for (int i = 0; i < n; ++i)
	{
		from[i] = plan[i].from + 1;
		to[i] = plan[i].to + 1;
		mass[i] = plan[i].mass;
	}
	
	return Rcpp::DataFrame::create(Rcpp::Named("from") = from, Rcpp::Named("to") = to, Rcpp::Named("mass") = mass);
}

// [[Rcpp::export(transportCost1d)]]
double RcppTransportCost1d(const Rcpp::NumericVector& a, const Rcpp::NumericVector& b, 
                           const Rcpp::NumericVector& wa, const Rcpp::NumericVector& wb,
                           const double p, const double threshold = 1e-15)
{
	return transportCost1d(a, wa, wa.size(), b, wb, wb.size(), p, threshold);
}
