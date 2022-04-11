/*
 * Implementation of quadratic program generator from:
 * 
 * P. H. Calamai, L. N. Vicente and J. J. Judice, A new technique for generating quadratic programming test programming test problems
 * Mathematical Programming 61 (1993), 215-231
 * 
 * 
 * Author: Wendel Melo
 * 
 * Date: 21-May-2013
 * 
*/

#ifndef _IQD_QPGEN_HPP
#define _IQD_QPGEN_HPP


# include <cstdio>
# include <cmath>
#if IQD_HAVE_RANDOMA
	# include "sfmt.h"
#endif
# include "iquad.hpp"
# include "IQD_tools.hpp"


namespace iquad
{
	
	
    enum IQD_QPGEN_TYPES
    {
		IQD_QPGEN_CONCAVE 	= 491,
		IQD_QPGEN_BILINEAR	= 492,
		IQD_QPGEN_CONVEX	= 493,
		IQD_QPGEN_INDEFINITE	= 494
    };
    
    
    
    
    


    class IQD_QPGenerator
    {
		#if IQD_HAVE_RANDOMA
			IQD_Random random;
		#endif
			
			int seedToRandomNumbers;

		public:
		
		//input parameters
		int L;
		int condd; //to condition number of D (10**condm)
		int lowerBoundToGenerateD;
		double denv; //density of v (denv*v gives the number of nonzeros at v)
		double factor; //we multiply all coeficients in the problem by this factor...
		
		IQD_QPGenerator();
		
		IQD_QPGenerator(const int seed);
		
		int initialize(const int seed);
		
		int setSeed(const int *seed);
	
		int generateProblem(const int nConcaveSubProbs, const int nBilinearSubProbs, const int nConvexSubProbs, IQD_IQuadProb &prob);
    };



}


#endif
