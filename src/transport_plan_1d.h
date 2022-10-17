#ifndef GRIDOT_TRANSPORT_PLAN_1D_H_INCLUDED
#define GRIDOT_TRANSPORT_PLAN_1D_H_INCLUDED

#include<vector>
#include"transport_struct.h"

template<typename Vector1, typename Vector2>
std::vector<Transport> transportPlan1d(const Vector1& mu, const int n, const Vector2& nu, const int m, const double threshold)
{
	std::vector<Transport> plan;
	plan.reserve(n + m - 1);

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
				return plan;
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
				return plan;
			}
		}

		t = std::min(r, s);

		plan.emplace_back(i, j, t);

		r -= t;
		s -= t;

	}
}

#endif
