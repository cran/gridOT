#ifndef GRIDOT_LINEAR_INTERPOLATION_DUAL_H_INCLUDED
#define GRIDOT_LINEAR_INTERPOLATION_DUAL_H_INCLUDED

#include"step_function_dual.h"

// dual calculator that linearly interpolates between the breakpoints
class LinearInterpolationDual : public StepFunctionDual
{
	protected:
		double quantileDiff(const double csum, const int i, const int k) const
		{
			// factor of linear interpolation
			double t = (csum - breaks[i - 1]) / (breaks[i] - breaks[i - 1]);
	
			return diff((1 - t) * y[i - 1] + t * y[i], k);
		}
		
	public:
		using StepFunctionDual::StepFunctionDual;
};

#endif
