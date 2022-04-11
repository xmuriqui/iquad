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


#ifndef _IQD_IQUAD_BASIC_HPP
#define _IQD_IQUAD_BASIC_HPP

#include <cstdlib>
#include <cstring>

#include <iostream>
#include <list>


//#include "SPM_SparseMatrix.hpp"
#include "MIP_minlpProblem.hpp"
#include "OPT_solvers.hpp"

#include "IQD_config.hpp"


namespace iquad
{
	
	
	
	
	
    #define IQD_VERSION "0.2.2"
	#define IQD_IDEALIZER "Jon Lee (University of Michigan, USA)"
    #define IQD_AUTHOR "Wendel Melo"
    #define IQD_EMAIL "wendelmelo@ufu.br"
	#define IQD_AUTHOR_FILIATION "Federal University of Uberlandia, Brazil"
	#define IQD_COLLABORATORS "Marcia Fampa (Federal University of Rio de Janeiro, Brazil)"
    #define IQD_DATE_VERSION "19-Nov-2015"
    
    #define IQD_MULTITHREADING 1
    
    #define IQD_DEBUG_MODE WAXM_DEBUG_MODE
	
	//if IQD_DEBUG_MODE is off, automatically IQD_NODES_DEBUG_MODE will be off too. otherwise, we decide if turn on or off
	#if IQD_DEBUG_MODE
		#define IQD_NODES_DEBUG_MODE 0
	#else
		#define IQD_NODES_DEBUG_MODE 0
	#endif
	
	#define IQD_INFINITY_TO_DIRECTIONS_BOUNDS 1e10;//used when a direction is unbounded in the split proccess
    #define IQD_INFINITY MIP_INFINITY
    //#define IQD_ZERO_TOLERANCE 1.0e-6
	
	#define IQD_TEST_CONVEX_PROBLEM 0
	#define IQD_CONVEX_DIR_MOVING "_convex"
	
	#define IQD_USE_MURIQUI_ON_INT_ORIGINAL_PROBLEM 0
	
	#define IQD_APPLY_JON_TRANSFORMATION 0
	
	#define IQD_USE_CPLEX_ON_NONCONVEX_PROBLEM 0
	
	#define IQD_MURIQUI_GAP_MIN_TREATMENT 1
	
	
	#define IQD_FILTER_PROBLEMS 0
	
	
	
	//typedef spm::SPM_SparseElement<double> IQD_SparseElement;
	//typedef spm::SPM_SparseRow<double> IQD_SparseRow; 
    //typedef spm::SPM_SparseMatrix<double> IQD_SparseMatrix;
    
    typedef minlpproblem::MIP_SparseMatrix IQD_SparseMatrix;
    
    
    template <class myClass>
    inline myClass IQD_max(const myClass &a, const myClass &b ) 
    {
		return a >= b ? a : b;
    }
    
    template <class myClass>
    inline myClass IQD_min(const myClass &a, const myClass &b)
    {
		return a <= b ? a : b;
    }
    
    template <class myClass>
    inline myClass IQD_abs(const myClass &a)
    {
		return a >= 0 ? a : -a;
    }
    
    
    
    enum IQD_RETURNCODES
    {
		IQD_OPTIMAL_SOLUTION 		= 0,
		
		IQD_UNDEFINED				= -101,
		IQD_SPLITING_ERROR			= -102,
		IQD_INFEASIBLE_PROBLEM 		= -103,
		IQD_UNBOUNDED_PROBLEM 		= -104,
		
		IQD_MEMORY_ERROR 			= -105,
		IQD_BAD_DEFINITIONS 		= -106,
		IQD_MAX_ITERATIONS_STOP 	= -107,
		IQD_MAX_TIME_STOP 			= -108,
		IQD_GENERIC_SOLVER_ERROR 	= -109,
		IQD_UNDEFINED_ERROR			= -110,
		IQD_LIBRARY_NOT_AVAILABLE 	= -111,
		IQD_CALLBACK_FUNCTION_ERROR	= -112,
		
		IQD_SDP_SOLVING_ERROR 		= -113,
		IQD_QCP_SOLVER_ERROR 		= -114,
		IQD_NLP_SOLVER_ERROR 		= -115,
		IQD_SOLVER_ERROR			= -116,
		
		IQD_NO_BOUNDED_FEASIBLE_SET = -117,
		
		IQD_NAME_ERROR				= -118,
		IQD_VALUE_ERROR				= -119
    };
    
    enum IQD_VARTYPE
    {
		IQD_VT_CONTINUOUS =	201,
		IQD_VT_BINARY = 	202,
		IQD_VT_INTEGER = 	203,
		IQD_VT_CONTINTEGER = 	204 //a special type of variable treated as continuous in the subproblems, but branched like a integer variable...
    };
    
    /*enum IQD_CONSTRTYPE
    {
		IQD_CT_EQ = QCO_CT_EQ,	// =
		IQD_CT_GREQ = QCO_CT_GREQ,	// >=
		IQD_CT_LEEQ = QCO_CT_LEEQ	// <=
    }; */
	
    //ways to SPLIT the non convex quadratic objective function...
    enum IQD_SPLITQUAD
    {
		//IQD_SQ_SDP 		= 101,	//via semidefinite programming
		IQD_SQ_SCHUR 		= 102,	//via real schur decomposition
		IQD_SQ_DIAG_SDP		= 103,	//via diagonal matrix using semidefinite programming
		//IQD_SQ_DIAG_SDP2	= 203,	//via diagonal matrix using semidefinite programming
		IQD_SQ_DIAG_DOM		= 104,	//via Diagonally Dominant
		IQD_SQ_IDENTITY		= 105,  //via identity matrix
		//IQD_SQ_USER_MATRIX	= 106,  //via matrix from user
		//IQD_SQ_USER_DIAGONAL	= 107,	//via diagonal from user
		
		//IQD_SQ_2_BLOCK_SDP	= 108,  //n/2 2 x 2 sdp blocks...
		//IQD_SQ_MIN_EIG_SDP = 109, //semidefinite programming minimizing highest eigenvalue
		//IQD_SQ_MILP_BASIS = 110, //via MILP problem on schur decomposition
		IQD_SQ_MIX_SCHUR_DIAG_SDP = 111
    };
	
	
	//strategy to handle integer problems
	enum IQD_INTSTRATEGY
	{
		IQD_IS_ON_BB_ENUMERATION = 500,
		IQD_IS_ON_EACH_SUBPROBLEM = 501
	};
	
	
	
	
	
	
	
    /*enum IQD_PROBLEMTYPE
    {
		IQD_PT_QP	= 140,
		IQD_PT_MIQP	= 141,
		IQD_PT_QCP	= 142,
		IQD_PT_MIQCP	= 143,
		IQD_PT_NLP	= 144,
		IQD_PT_MINLP	= 145
    }; */
	
	
	
    enum IQD_SDP_SOLVERS
    {
		IQD_SS_CSDP = 150,
		IQD_SS_MOSEK= optsolvers::OPT_MOSEK
	};
	
	
	enum IQD_MILP_SOLVERS
	{
		IQD_IP_GLPK		= optsolvers::OPT_GLPK,
		IQD_IP_CPLEX	= optsolvers::OPT_CPLEX,
		IQD_IP_GUROBI	= optsolvers::OPT_GUROBI,
		IQD_IP_XPRESS	= optsolvers::OPT_XPRESS,
		IQD_IP_MOSEK	= optsolvers::OPT_MOSEK
	};
	
	
	enum IQD_QP_SOLVERS
	{
		IQD_QP_CPLEX	= optsolvers::OPT_CPLEX,
		IQD_QP_GUROBI	= optsolvers::OPT_GUROBI,
		IQD_QP_XPRESS	= optsolvers::OPT_XPRESS,
		IQD_QP_MOSEK	= optsolvers::OPT_MOSEK,
		IQD_QP_IPOPT	= optsolvers::OPT_IPOPT,
		IQD_QP_WORHP	= optsolvers::OPT_WORHP
	};
	
	
	enum IQD_QCP_SOLVERS
	{
		//IQD_QCP_CPLEX	= IQD_QP_CPLEX,
		//IQD_QCP_GUROBI	= IQD_QP_GUROBI,
		IQD_QCP_XPRESS	= IQD_QP_XPRESS,
		IQD_QCP_MOSEK	= IQD_QP_MOSEK,
		IQD_QCP_IPOPT	= IQD_QP_IPOPT,
		IQD_QCP_WORHP	= IQD_QP_WORHP
	};
	
	
	enum IQD_NLP_SOLVERS
	{
		IQD_NLP_MOSEK	= IQD_QCP_MOSEK,
		IQD_NLP_IPOPT	= IQD_QCP_IPOPT,
		IQD_NLP_WORHP	= IQD_QCP_WORHP
	};
	
	
	enum IQD_MINLP_SOLVERS
	{
		IQD_MS_MURIQUI	=	180,
		IQD_MS_BONMIN	=	181
    };
    
	
	//that enumeration is used when user wants check the relaxation problem being solved by iquad in the root node
	enum IDQ_RELAX_PROB_WRITING_MODE
	{
		IQD_RPWM_NO_WRITING 		= 800,
		IQD_RPWM_WRITING
	};
	
	
	int IQD_getStrEnumValue(const char *svalue, IQD_SDP_SOLVERS &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_QP_SOLVERS &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_QCP_SOLVERS &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_MILP_SOLVERS &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_NLP_SOLVERS &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_MINLP_SOLVERS &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_SPLITQUAD &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_INTSTRATEGY &value);
	
	
    
    
    /* That class represents the problem addressed:
     *
     * Min c'x + 0.5 x'Qx + f(x) + d
     * subject to:
	 * 
     * 	l_c <= a_i x + 0.5 x'Q_i x + g_i(x) <= u_c
     * 	
     * 	l_x <= x <= u_x
     * 
     * 	x_i is integer for i \in I
     * 
     * 
     * We observe that functions f(.) and g_i(.) are convex. 
     * The unique nonconvex term in the problem above is the quadratic term in 
     * objective function 0.5 x'Qx.
     *  
     */
	
// 	
    typedef minlpproblem::MIP_MINLPProb IQD_IQuadProb;
	typedef minlpproblem::MIP_NonLinearEval IQD_NonLinearEval;
	
	typedef optsolvers::OPT_GeneralSolverParams IQD_GeneralSolverParams;
	
	typedef minlpproblem::MIP_PROBLEMTYPE IQD_PROBLEMTYPE;
	//typedef enum MIP_VARTYPE IQD_VARTYPE;
	
	
	
	
	
	
	
	
	
    
}



#endif
