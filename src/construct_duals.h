#ifndef GRIDOT_CONSTRUCT_DUALS_H_INCLUDED
#define GRIDOT_CONSTRUCT_DUALS_H_INCLUDED

#include<RcppArmadillo.h>
#include<vector>
#include<memory>
#include"dual.h"
#include"step_function_dual.h"
#include"linear_interpolation_dual.h"

// methods for constructing the three different types of dual solution calculators:
// discrete, epsilon-discrete and epsilon-histogram	

template<typename Vector>
std::unique_ptr<Dual> constructDiscrete(const arma::vec& x, const arma::vec& y, Vector&& weights, 
                                        const double marginalWeight, const double p, const double rightMargin)
{
	const int n = weights.n_elem;
	
	std::vector<double> y0;
	std::vector<double> breaks;
	y0.reserve(n);
	breaks.reserve(n);
	
	double sum = rightMargin;
	for (int i = 0; i < n; ++i)
	{
		// only points with mass
		if (weights[i] > 0)
		{
		  y0.push_back(y[i]);
		  
		  sum += weights[i];
		  
		  breaks.push_back(sum);
		}
	}
	
	return std::unique_ptr<Dual>(new StepFunctionDual(x, y0, breaks, p));
}

template<typename Vector>
std::unique_ptr<Dual> constructEpsilonDiscrete(const arma::vec& x, const arma::vec& y, Vector&& weights, 
                                               const double marginalWeight, const double p, 
                                               const double rightMargin, const double eps)
{
	// remove points with no mass
	arma::uvec ind = arma::find(weights > 0);
	arma::vec weights0 = weights.elem(ind);
	arma::vec y1 = y.elem(ind);
	
	// edge case
	if (ind.n_elem == 1) 
	{
		std::vector<double> y2 = { y1[0] };
		std::vector<double> breaks2 = { 2.0 };
		
		return std::unique_ptr<Dual>(new StepFunctionDual(x, y2, breaks2, p));
	}
	
	const int n = weights0.n_elem - 1;
	
	// additional mass everywhere
	const double epsWeight = eps / (y1.back() - y1[0]) * marginalWeight;
	
	std::vector<double> y0(2 * n);
	std::vector<double> breaks(2 * n);
	
	const double invEps = 1 - eps;
	double sum = rightMargin;
	
	int k = 0;
	for (int i = 0; i < n; ++i)
	{
		// point mass
		sum += weights0[i] * invEps;
		
		y0[k] = y1[i];
		breaks[k++] = sum;
		
		// integrated additional eps mass
		sum += (y1[i + 1] - y1[i]) * epsWeight;
		
		y0[k] = y1[i + 1];
		breaks[k++] = sum;
	}
	
	return std::unique_ptr<Dual>(new LinearInterpolationDual(x, y0, breaks, p));
}

template<typename Vector>
std::unique_ptr<Dual> constructEpsilonHistogram(const arma::vec& x, const arma::vec& y, Vector&& weights, 
                                                const double marginalWeight, const double p, const double rightMargin,
                                                const double eps, const double width)
{
	// remove points with no mass
	arma::uvec ind = arma::find(weights > 0);
	arma::vec weights0 = weights.elem(ind);
	arma::vec y1 = y.elem(ind);
	
	// edge case
	if (ind.n_elem == 1) 
	{
		std::vector<double> y2 = { y1[0] };
		std::vector<double> breaks2 = { 2.0 };

		return std::unique_ptr<Dual>(new StepFunctionDual(x, y2, breaks2, p));
	}
	
	const int n = weights0.n_elem;
	
	// additional mass everywhere
	const double epsWeight = eps / (y1.back() - y1[0] + width) * marginalWeight;
	
	std::vector<double> y0(2 * n);
	std::vector<double> breaks(2 * n);

	const double h = width / 2;
	const double invEps = 1 - eps;	

	// extra mass per point bin
	const double epsExtra = epsWeight * width;
	double sum = rightMargin;
	
	int k = 0;
	
	// left bin point
	y0[k] = y1[0] - h;
	breaks[k++] = sum; 
	
	// mass in bin
	sum += weights0[0] * invEps + epsExtra;
	
	// right bin point
	y0[k] = y1[0] + h;
	breaks[k++] = sum; 
	
	for (int i = 1; i < n; ++i)
	{
		// mass between point bins
		sum += (y1[i] - y1[i - 1]) * epsWeight;
		
		// left bin point
		y0[k] = y1[i] - h;
		breaks[k++] = sum;
		
		// mass in bin
		sum += weights0[i] * invEps + epsExtra;
		
		// right bin point
		y0[k] = y1[i] + h;
		breaks[k++] = sum;
	}
	
	return std::unique_ptr<Dual>(new LinearInterpolationDual(x, y0, breaks, p));
}

#endif
