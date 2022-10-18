#ifndef GRIDOT_DUAL_H_INCLUDED
#define GRIDOT_DUAL_H_INCLUDED

#include<vector>

// class that encapsulates information to calculate the dual solution given the new weights of the first measure
class Dual
{
	protected:
		// break points of the quantile function
		const std::vector<double> breaks;
		
		// calculates the increase of the dual solution 
		// @param csum current cumulative weight sum
		// @param i index of the break point
		// @param k index of the dual
		double virtual quantileDiff(const double csum, const int i, const int k) const = 0;
		
		double virtual quantileDiff(const int i, const int k) const = 0;

	public:
		Dual(const std::vector<double>& breaks) : breaks(breaks)
		{
		}
		
		virtual ~Dual() = default;
		
		// calculates the dual solution
		// @param weightStart start of the new weights of the first measure
		// @param weightEnd end of the new weights of the first measure
		// @param dualStart start of where the new dual solution is written to
		template<typename WeightIterator, typename DualIterator>
		void calculateDual(WeightIterator&& weightStart, WeightIterator&& weightEnd, DualIterator&& dualStart) const;
};

template<typename WeightIterator, typename DualIterator>
void Dual::calculateDual(WeightIterator&& weightStart, WeightIterator&& weightEnd, DualIterator&& dualStart) const
{
	DualIterator dual = dualStart;
	*dual = 0.0;
	++dual;
	
	WeightIterator weight = weightStart;
	weightEnd = std::prev(weightEnd); 
	if (weight == weightEnd)
	{
		return;
	}
	
	const int n = breaks.size();
	int k = 0;
	
	double csum = *weight;
	double br = breaks[0];
	double dsum = 0.0;
	
	// all that are <= the first break are mapped to 0
	while (csum <= br)
	{
		dsum += quantileDiff(0, k);
		*dual = dsum;
		
		++dual;
		++k;
		++weight;
		if (weight == weightEnd)
		{
			return;
		}
		
		csum += *weight;
	}
	
	for (int i = 1; i < n; ++i) 
	{
		br = breaks[i];
		
		// find first break that is greater
		while (csum <= br)
		{	
			dsum += quantileDiff(csum, i, k);
			*dual = dsum;
			
			++dual;
			++k;
			++weight;
			if (weight == weightEnd)
			{
				return;
			}
			
			csum += *weight;
		}
	}
	
	// all that are > last break are mapped to n - 1
	while (weight != weightEnd)
	{
		dsum += quantileDiff(n - 1, k);
		*dual = dsum;
		
		++dual;
		++k;
		++weight;
	}
}

#endif
