#ifndef GRIDOT_STEP_FUNCTION_DUAL_H_INCLUDED
#define GRIDOT_STEP_FUNCTION_DUAL_H_INCLUDED

#include<vector>
#include<cmath>
#include<RcppArmadillo.h>
#include"dual.h"

// dual calculator that uses constant interpolation between breakpoints
class StepFunctionDual : public Dual
{
	protected:
		// support points of the first measure
		const arma::vec& x;
		// support points of the second measure
		const std::vector<double> y;
		// power of the cost
		const double p;

		double diff(const double yi, const int k) const
		{
			// increase of the dual solution
			return std::pow(std::fabs(yi - x[k + 1]), p) - std::pow(std::fabs(yi - x[k]), p);
		}

		double virtual quantileDiff(const double csum, const int i, const int k) const
		{
			return quantileDiff(i, k);
		}
		
		double quantileDiff(const int i, const int k) const
		{
			return diff(y[i], k);
		}

	public:
		StepFunctionDual(const arma::vec& x, const std::vector<double>& y, 
		                 const std::vector<double>& breaks, const double p) : Dual(breaks), x(x), y(y), p(p)
		{
		}
};

#endif
