#include "AccD.hpp"

#include <iostream>
#include <cmath>
using namespace utility;
// compile with:clang -std=c++14 test_AccD.cpp -I[PATH_TO_BOOST]  -I[PATH_TO_SPRNG5]/include -L[PATH_TO_SPRNG5]/lib -lboost_timer -lboost_system -lstdc++ -lm


inline RealT f1(std::array<RealT,1> x)
{
	return x[0]; //std::sin(x[0]);
}

inline RealT f2(std::array<RealT,2>x)
{
	return std::sin(x[0])*x[1];
}

int main(void){
	//naive 1D test
	double s = 0.0;
	double si = 4.0/10000.0;
	for(double i=0.0;i<4.0;i+=si) s+= i;
	s = 4.0*s/10000.0;
	std::cout << "direct: " << s << std::endl;

	// 1D test
	KInt<1> k1;
	std::array<RealT,1> min1 = {0};
	std::array<RealT,1> max1 = {4};
	std::array<unsigned long,1> N1 = {10000};
	std::cout << "weighted kahan: " << k1.sumKPoints(f1, min1, max1, N1) << std::endl;;
	
	// 2D test
	KInt<2> k2;
	std::array<RealT,2> min2 = {0,0};
	std::array<RealT,2> max2 = {4,9};
	std::array<unsigned long,2> N2 = {10000,10000};
	std::cout << "2D weighted kahan: " <<k2.sumKPoints(f2, min2, max2, N2) << std::endl;;

}
