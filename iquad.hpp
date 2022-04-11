/*iquad.hpp: header to use iquad library. 
 * Here, we have basic definitions to iquad (Indefinity Quadratic) solver
 * 
 * Problem addressed:
 * 
 * Min c'x + 0.5 x'Qx + f(x)
 * s. t.:
 *  g(x) <= 0
 * 
 *	A x + 0.5 x'QQ x {<=, =, >=} b
 *	
 *	l <= x <= u
 * 
 *	x_i is integer for i in I
 * 
 * 
 * functions f(x) and g(x) are convex. The unique nonconvex part is restricted 
 * to a quadratic term in objective function 0.5x'Qx. All constraints are convex.
 * 
 * 
 */


#ifndef _IQUAD_HPP
#define _IQUAD_HPP



#include "IQD_iquadBasic.hpp"

#include "IQD_branchAndBound.hpp"






#endif




