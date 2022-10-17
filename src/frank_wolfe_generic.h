#ifndef GRIDOT_FRANK_WOLFE_GENERIC_H_INCLUDED
#define GRIDOT_FRANK_WOLFE_GENERIC_H_INCLUDED

#include<RcppArmadillo.h>
#include<memory>
#include<vector>
#include"network_simplex_simple.h"
#include"dual.h"
#include"transport_cost.h"

// function that constructs a dual calculator
template<typename Vector, typename... Args>
using DualConstructor = std::unique_ptr<Dual>(*)(const arma::vec& x, const arma::vec& y, Vector&& weights, 
                                                 const double marginalWeight, const double p, Args... args);

// performs the Frank-Wolfe algorithm with the specified dual calculators
// @param x1 first marginal support points of the first measure
// @param x2 second marginal support points of the first measure
// @param mu first measure
// @param y1 first marginal support points of the second measure
// @param y2 second marginal support points of the second measure
// @param nu second measure
// @param p1 first power of the cost
// @param p2 second power of the cost
// @param startPivot starting measure for the algorithm
// @param maxIt maximum number of iterations
// @param tol stopping criterion for the dual gap
// @param threads number of threads the algorithm uses
// @param returnIt whether or not the cost and dualgap in each iteration are returned
template<typename... DualArgs>
Rcpp::List frankWolfeGeneric(const arma::vec& x1, const arma::vec& x2, const arma::mat& mu, 
                             const arma::vec& y1, const arma::vec& y2, const arma::mat& nu, 
                             const double p1, const double p2, const arma::mat& startPivot,
                             const int maxIt, const double tol, const int threads, const bool returnIt,
                             const DualConstructor<arma::rowvec, DualArgs...>& dualConstructor1,
                             const DualConstructor<arma::colvec, DualArgs...>& dualConstructor2,
                             DualArgs... dualArgs)
{
#ifdef _OPENMP	
	omp_set_num_threads(threads);
#endif

	const int n2 = x2.n_elem;
	const int m1 = y1.n_elem;
	
	arma::colvec nu1 = arma::sum(nu, 1);
	arma::rowvec mu2 = arma::sum(mu, 0);
	
	std::vector<std::unique_ptr<Dual>> duals1;
	std::vector<std::unique_ptr<Dual>> duals2;
	duals1.reserve(m1);
	duals2.reserve(n2);
	
	// construct the dual calculators
	for (int r = 0; r < m1; ++r)
	{
		duals1.push_back(dualConstructor1(x2, y2, nu.row(r), nu1[r], p2, dualArgs...));
	}
	
	for (int j = 0; j < n2; ++j)
	{
		duals2.push_back(dualConstructor2(y1, x1, mu.col(j), mu2[j], p1, dualArgs...));
	}

	lemon::FullBipartiteDigraph graph(m1, n2);
	lemon::NetworkSimplexSimple<lemon::FullBipartiteDigraph, double, double, long long> 
	net(graph, true, m1 + n2, m1 * n2);
	
	// set the weights
	std::vector<double> weights1(m1);
	std::vector<double> weights2(n2);
	for (int r = 0; r < m1; ++r)
	{
		weights1[graph.nodeFromId(r)] = nu1[r];
	}
	for (int j = 0; j < n2; ++j)
	{
		weights2[graph.nodeFromId(j)] = -mu2[j];
	}
	net.supplyMap(weights1.data(), m1, weights2.data(), n2);
	
	arma::mat cost(m1, n2, arma::fill::none);
	arma::mat s(m1, n2, arma::fill::none);
	
	arma::mat xi = startPivot;

	double dualgap;
	bool conv = false;

	std::vector<double> costs;
	std::vector<double> gaps;
	if (returnIt)
	{
	    costs.reserve(maxIt);
	    gaps.reserve(maxIt);
	}

	for (int step = 1; step <= maxIt; ++step) 
	{
		// calculate the dual solutions
		#pragma omp parallel for
		for (int j = 0; j < n2; ++j)
		{
			duals2[j]->calculateDual(xi.begin_col(j), xi.end_col(j), cost.begin_col(j));
		}

		#pragma omp parallel for
		for (int r = 0; r < m1; ++r)
		{
			arma::rowvec dual1(n2, arma::fill::none);

			duals1[r]->calculateDual(xi.begin_row(r), xi.end_row(r), dual1.begin());
			
			cost.row(r) += dual1;
		}
		
		// set new costmatrix
		int64_t costId = 0;
		for (int r = 0; r < m1; ++r)
		{
			for (int j = 0; j < n2; ++j)
			{
				net.setCost(graph.arcFromId(costId++), cost(r, j));
			}
		}
		
		net.run();
		
		for (int64_t r = 0; r < m1; ++r) 
		{
			for (int64_t j = 0; j < n2; ++j)
			{
				// just store difference, not the plan itself
				s(r, j) = net.flow(graph.arcFromId(r * n2 + j)) - xi(r, j);
			}
		}
		
		// take abs to make sure its not negative
		dualgap = std::fabs(arma::accu(cost % s));
		
		xi += 2.0 / (step + 1.0) * s;
		
		if (returnIt)
		{
			costs.push_back(transportCost(x1, x2, mu, y1, y2, nu, p1, p2, xi));
		    gaps.push_back(dualgap);
		}

		if (dualgap <= tol)
		{
			conv = true;

			break;
		}
	}
	
	if (returnIt)
	{
		return Rcpp::List::create(Rcpp::Named("pivot") = xi, Rcpp::Named("conv") = conv, 
		                          Rcpp::Named("costs") = costs, Rcpp::Named("dualgaps") = gaps);
	}

	return Rcpp::List::create(Rcpp::Named("pivot") = xi, Rcpp::Named("conv") = conv);
}

#endif
