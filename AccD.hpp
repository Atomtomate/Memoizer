#ifndef ACC_D_H_
#define ACC_D_H_
#include <bitset>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>


/*!	This header provides a general templated for the accumulation of function values over arbitrary dimensions.
 *
 *	As of now, this is limited to the mean value, but generalization should be straight forward
 */	
namespace utility{
template<unsigned int D>
using VecD = std::array<RealT, D>;
using AccT= boost::accumulators::accumulator_set<RealT, boost::accumulators::stats<boost::accumulators::tag::mean> >;

/*!	this namespace is needed to avoid partial specialization of member functions inside AccD
 *	the static member expandLoop recursivly opens loops over all dimension
 */
namespace detail
{
	template<unsigned int D, unsigned int ND>
	struct Internal
	{
		static void expandLoop(RealT (*integrand)(VecD<D> x),const VecD<D> &min,const VecD<D> &incs, const std::array<unsigned long, D> &N,\
						VecD<D> &xVec, AccT &acc)
		{
			RealT xi=min[ND];
			for(unsigned int n=0; n<N;n++)
			{
				xVec[ND] = xi;
				Internal<D,ND-1>::expandLoop(integrand, min, incs, N, xVec, acc);
				xi+=incs[ND];
			}
		}
	};

	template<unsigned int D>
	struct Internal<D,0>
	{
		static void expandLoop(RealT (*integrand)(VecD<D> x),const VecD<D> &min,const VecD<D> &incs, const std::array<unsigned long, D> &N,\
							VecD<D> &xVec, AccT &acc)
		{
			RealT xi=min[0];
			for(unsigned int n=0; n<N;n++)
			{
				xVec[0] = xi;
				acc(integrand(xVec));
				xi+=incs[0];
			}
		}
	};
}

template<unsigned int D>
class AccD
{
	public:
		template<unsigned int Di, unsigned int NDi>
		using expL = detail::Internal<Di,NDi>;
		//TODO: 2 test types: 3 or 4 known sums, mean over >1000 VERY large OR mean over >MAXINT number of samples numbers (check for overflow)
				RealT expandLoop(RealT (*integrand)(VecD<D> x), const VecD<D> &min, const VecD<D> &max, const std::array<unsigned long, D> &N) const
		{
			AccT acc;		// accumulate mean (1/N sum_x1 ... sum_xD f(x1,...,xD))
			std::bitset<D> stage(0);												// this is needed to emulate nested loops
			VecD<D> incs; 												// increments in each dimension
			VecD<D> xVec 	= min;
			for(unsigned int d=0;d<D;d++) incs[d] = (max[d]-min[d])/N[d];

			detail::Internal<D,D>::expandLoop(integrand, min, incs, N, xVec, acc);	
			return boost::accumulators::mean(acc);
		}

	private:

};

		
}
#endif
