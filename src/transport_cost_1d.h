#ifndef GRIDOT_TRANSPORT_COST_1D_H_INCLUDED
#define GRIDOT_TRANSPORT_COST_1D_H_INCLUDED

#include<cmath>

template<typename Vector1, typename Vector2, typename Vector3, typename Vector4>
double transportCost1d(const Vector1& x, const Vector2& mu, const int n, const Vector3& y, const Vector4& nu, const int m, 
                       const double p, const double threshold)
{
	int i = 0;
	int j = 0;

	double r = mu[i];
	double s = nu[j];
	double t = 0;

	double cost = 0;

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
				return cost;
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
				return cost;
			}
		}

		t = std::min(r, s);

		cost += t * std::pow(std::fabs(x[i] - y[j]), p);

		r -= t;
		s -= t;

	}

	return cost;
}

#endif
