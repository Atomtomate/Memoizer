#ifndef KINT_H_
#define KINT_H_

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/framework/accumulator_base.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/framework/parameters/sample.hpp>
#include <boost/accumulators/framework/parameters/weight.hpp>
#include <boost/accumulators/framework/accumulators/external_accumulator.hpp>
#include <boost/accumulators/statistics/weighted_sum.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>

#include <iostream>

/*!	This class provides methods for integration over momenta or energies
 *	
 *	Problem:
 *		Parameters are: 
 *		- a: lattice spacing
 *		- \f$ D(\eps) \f$: density of states (DOS)
 *		- Momentum, energy dependent function. In most cases this will be the impurity Green's \f$ G(k,i\omega_n) \f$
 *		- Some constant offset. In most cases \f$ i \omega_n + \mu - \Simga(i \omega_n)\f$
 *		- function and limits for the integrand e.g.: 
 *		\f[
 *			\int\limits_{-\frac{\pi}{a}}^\frac{\pi}{a} d^D k \frac{1}{i \omega_n - \epsilon_k + mu - \Sigma (i \omega_n)} 
 *		\f]
 *		or
 *		\f[
 *			\int\limits_{-\infty}^\infty d^D \epsilon \frac{d(\epsilon)}{i \omega_n - \epsilon_k + mu - \Sigma (i \omega_n)} 
 *		\f]
 *
 *	Methods:
 *
 *	- summation (this constructs some general weighted D dimensional sum)
 *		- riemann sum
 *		- gaussian smearing
 *		.
 *
 *	- reformulation to ODE and solution via boost ODE int
 * 
 *	- Linear tetrahedral method for 1,2,3 D
 * 		- issues: can break symmetry, gamma point must be included
 *		- use MFEM?
 *		.
 *
 *	- TODO: let user decide wether to use kahan summation or not
 */	
namespace utility{
using KahanTag = boost::accumulators::tag::weighted_sum_kahan;
using KahanSumAccT = boost::accumulators::stats<KahanTag>;
template <typename T>
using AccT = boost::accumulators::accumulator_set<T, KahanSumAccT, T >;


namespace detail
{
	/*!	this recursively constructs the depth of the nested for loops
	 *	for the innermost loop the partially specialized struct Internal<D,0> is called
	 */
	template<unsigned int D, unsigned int ND, typename T, typename RetT >
	struct Internal
	{
		static void sumKPoints(RetT (*integrand)(std::array<T, D> x),const std::array<T, D> &min,const std::array<T, D> &incs,\
			 const std::array<unsigned long, D> &N,std::array<T, D> &xVec, AccT<RetT> &acc)
		{
			RetT xi=min[ND];
			for(unsigned int n=0; n<N[ND];n++)
			{
				xVec[ND] = xi;
				Internal<D,ND-1, T, RetT>::sumKPoints(integrand, min, incs, N, xVec, acc);
				acc(boost::accumulators::extract_result<KahanTag>(acc), boost::accumulators::weight = incs[ND]);

				xi+=incs[ND];
			}
		}
	};

	template<unsigned int D, typename T, typename RetT>
	struct Internal<D,0,T,RetT>
	{
		static void sumKPoints(RetT (*integrand)(std::array<T, D> x),const std::array<T, D> &min, \
			const std::array<T, D> &incs, const std::array<unsigned long, D> &N, std::array<T, D> &xVec, AccT<RetT> &acc)
		{
			RetT xi=min[0];
			for(unsigned int n=0; n<N[0];n++)
			{
				xVec[0] = xi;
				acc(integrand(xVec), boost::accumulators::weight = incs[0]);
				xi+=incs[0];
			}
		}
	};
}


//TODO: 2 test types: 3 or 4 standard intrgrals, mean over >1000 VERY large OR mean over >MAXINT number of samples numbers (check for overflow)
//TODO: let user specify weights, TODO: let user choose between kahan and normal sum

template<unsigned int D>
class KInt
{
	public:
		/*!	@brief	instantiates an accumulator for D nested for loops \f$ \sum\limits_{min_1}^{max_1} ... \sum\limits_{min_D}^{max_D} f \f$
		 *			with \f$ f(T x_1, ... ,T x_D) -> T \f$
		 *
		 *	@param	integrand	function pointer over which to sum
		 *	@param	min			array of start values (one element for each dimension)
		 *	@param	max			array of values for the upper limit (one element for each dimension)
		 *	@param	N			array with number of steps (one element for each dimension)
		 *
		 *	@return	accumulated value
		 */
		template<typename T, typename RetT>
		RetT sumKPoints(RetT (*integrand)(std::array<T, D> x), std::array<T, D> &min, \
			const std::array<T, D> &max, const std::array<unsigned long, D> &N) const
		{
			//acc<summand type, method, weight type>
			AccT<RetT> acc;
			std::array<T, D> incs; 						// increments in each dimension
			std::array<T, D> xVec 	= min;
			for(unsigned int d=0;d<D;d++) incs[d] = (max[d]-min[d])/N[d];

			detail::Internal<D,D-1,T,RetT>::sumKPoints(integrand, min, incs, N, xVec, acc);	
			return boost::accumulators::extract_result<KahanTag>(acc);
		}

	private:

};

		
}
#endif
