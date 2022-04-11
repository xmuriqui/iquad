
#ifndef _IQD_BRANCHANDBOUND_HPP
#define _IQD_BRANCHANDBOUND_HPP

#include <ctime>
#include <iostream>


#include "BBL_branchAndBound.hpp"
#include "iquad.hpp"


namespace iquad
{
	
    typedef branchAndBound::BBL_Array<double>	IQD_DoubleArray;
	//typedef BBL_Array<unsigned int>	IQD_UIntArray;
	typedef branchAndBound::BBL_ArraySize<unsigned int>	IQD_UIntArraySize;
    typedef branchAndBound::BBL_NodeBounds 	IQD_NodeBounds;
    //typedef BBL_Node 		IQD_Node;
    typedef branchAndBound::BBL_NodeList 	IQD_NodeList;
    typedef branchAndBound::BBL_NodeListManager IQD_NodeListManager;
	typedef branchAndBound::BBL_MTNodeListManager IQD_MTNodeListManager;
    
	typedef branchAndBound::BBL_Mutex IQD_Mutex;
	
	
	
    
    enum IQD_EXP_STRATEGY
    {
		IQD_ES_DEPTH			= branchAndBound::BBL_ES_DEPTH,
		IQD_ES_WIDTH			= branchAndBound::BBL_ES_WIDTH,
		IQD_ES_BEST_LIMIT		= branchAndBound::BBL_ES_BEST_LIMIT,
		IQD_ES_DEPTH_BEST_LIMIT = branchAndBound::BBL_EES_DEPTH_BEST_LIMIT
    };
    
    enum IQD_BRANCH_STRATEGY
    {
		//IQD_BS_LARGEST_GAP	= 501,
		//IQD_BS_LARGEST_DIF	= 502,
		
		
		IQD_BS_MAX_CONSTRS_VIOLATED	= 501,  //we choose y variable appearing in the max number of constraint violated...
		
		IQD_BS_MOST_VIOLATED_CONSTR = 502
	};
	
	
	
	
	enum IQD_AUX_BOUNDS_CALC_STRATEGY
	{
		IQD_ABCS_SUBPROBLEMS = 300,
		IQD_ABCS_ORIG_BOUNDS = 301
	};
	
    
    int IQD_getStrEnumValue(const char *svalue, IQD_EXP_STRATEGY &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_BRANCH_STRATEGY &value);
	
	int IQD_getStrEnumValue(const char *svalue, IQD_AUX_BOUNDS_CALC_STRATEGY &value);
	
	int IQD_getStrEnumValue(const char *svalue, IDQ_RELAX_PROB_WRITING_MODE &value);
	
	/*
	 * Class to solve the auxiliary problems to try imporve the bounds for auxiliary variable y when a feasible solution in the B&B procedure. Here, when a feasible is found, we solve:
	 * 
	 * 
	 *  \bar{l_y} = 	min y
				s.t.   x \in X 	(x is feasible)
					-( u_y - l_y ) y + u_y l_y   =  w       (from secant inequalities)
					c'x + 0.5 x' P x  +   0.5 -\lambda_i w_i   <= zu - \epsilon        (objective cut ) 
					l_y <= y <= u_y
      
      
      \bar{u_y} = 	max y
				s.t. :   x \in X     (x is feasible )
					  -( u_y - l_y ) y + u_y l_y   =  w       (from secant inequalities)
					  c'x + 0.5 x' P x  + 0.5 -\lambda_i w_i   <= zu       (objective cut )
					  l_y <= y <= u_y
	 * 
	 * 
	 * Note secant are equalities. Since w variables are not in the objective function, we need to fix w in the secant using equalities
	 * 
	 * 
	 */ 
	
	
	
	class IQD_RefProb2;
	class IQD_AuxBoundsUpdate;
	
	
	//New branch-and-bound implementatio using the BBL framework
	class IQD_BranchAndBound
    {
	public:
		
		/** input parameters to branch-and-bound procedure */
		
		bool in_build_split_without_vars_also;
		bool in_certain_of_nonconvex_Q_on_obj_f; //just set this flag to true if you are sure Q is nonconvex. It will avoid spending time to check if you are using diagonal splits. If this flag is false, iquad will check ayway...
		bool in_considerate_integrality_at_calculating_bounds;
		
		bool in_reorganize_lists;
		bool in_set_secants_as_equality_constraints; //if set as true, set the secant constraints as equalities instead of inequalities. Some people say inequalities constraints are better to some solvers, but equalities  can avoid some numerical problems in cases where we can have some lambda closed to zero, specially in mix way of split (IQD_SQ_MIX_SCHUR_DIAG_SDP). if set as false, secant constraints will be set as false...
		bool in_set_special_solver_parameters;
		bool in_priorize_obj_f_on_branch_choosing;
		bool in_recalculate_aux_bounds_on_zu_updt;
		bool in_use_real_dif_on_branch_choosing;
		bool in_use_heuristics_on_integer_probs;
		
		IQD_SDP_SOLVERS in_sdp_solver;
		//IQD_MILP_SOLVERS in_milp_solver;
		IQD_QP_SOLVERS in_qp_solver;
		IQD_QCP_SOLVERS in_qcp_solver;
		IQD_NLP_SOLVERS in_nlp_solver;
		//IQD_MINLP_SOLVERS in_minlp_solver;
		
		IQD_SPLITQUAD in_split_way;
		IQD_EXP_STRATEGY in_exp_strategy;
		IQD_BRANCH_STRATEGY in_branch_strategy;
		IQD_INTSTRATEGY in_integrality_strategy;
		
		IQD_AUX_BOUNDS_CALC_STRATEGY in_aux_bounds_calc_strategy;
		
		IDQ_RELAX_PROB_WRITING_MODE in_relax_prob_writing_mode;
		
		unsigned int in_lists_reorganization_frequency;
		
		unsigned int in_number_of_node_sublists;
		unsigned int in_number_of_threads; 
		unsigned int in_printing_frequency;
		unsigned int in_print_level;
		
		long unsigned int in_max_iters;
		
		double in_absolute_convergence_tol;
		double in_alpha;						//parameter to define branch point 1.0 (use only the solution of relaxation), 0.0 (use only the middle point between bound). (0.0  1.0) a balance between both.
		//double in_eps_added_sd_solver_matrices;
		double in_alpha_to_integer_branch_choosing; //alpha to balance the index of integer branching:  1.0 (use only the integer gap), 0.0 (use only a possible objective difference). (0.0 1.0) a balance between both. This parameter is effective on diagonal splits
		
		double in_absolute_eps_added_to_diag_P;
		double in_relative_eps_added_to_diag_P;
		
		double in_factor_to_prior_obj_f_on_branch;
		
		double in_factor_to_recalculate_aux_bounds_on_zu_updt;
		
		double in_abs_feas_tol;
		double in_rel_feas_tol;
		
		double in_infinity_value;
		double in_integer_tol;
		
		double in_lower_bound;
		double in_max_time;
		double in_max_cpu_time;
		double in_upper_bound;
		double in_relative_convergence_tol;
		
		double in_zero_tol_to_split; //tolerance to consider some value zero in several cases inside computations...
		
		double in_hours_to_save_file;
		
		
		
		
		/** output parameters */
		bool out_feasible_sol;
		unsigned int out_number_of_feas_sols;
		unsigned long int out_number_of_iterations;
		unsigned long int out_number_of_open_nodes;
		double *out_best_sol;
		double out_best_obj;
		IQD_RETURNCODES out_return_code;
		int out_return_subcode;
		double out_cpu_time_to_fisrt_feas_sol;
		double out_splitting_cpu_time;
		double out_cpu_time, out_clock_time;
		double out_lower_bound, out_upper_bound;
		double out_root_node_lower_bound, out_root_node_upper_bound;
		
		
		
		//special parameter to muriqui minimization gap problems
		
		#if IQD_MURIQUI_GAP_MIN_TREATMENT
			bool mrq_gap_min;
		#endif
		
		
		
		IQD_BranchAndBound();
		
		~IQD_BranchAndBound();
		
		
		void resetInputParameters(void);
		
		void resetOutput(void);
		
		void desallocateSol(void);
		
		
		int run(IQD_IQuadProb& prob, IQD_GeneralSolverParams* solverParams = NULL, IQD_GeneralSolverParams* sdpParams = NULL);
		
		
		int setIntegerParameter(const char *name, const long int value);
		
		int setDoubleParameter(const char *name, const double value);
		
		//that method should be used to set enumeration parameters
		int setStringParameter(const char *name, const char *value);
		
		
		
		
	protected:
		
		bool useWarmStart;
		
		double timeStart;
		clock_t clockStart;
		
		//IQD_RefProb *refProb;
		//IQD_RefProb2 *refProb2;
		
		
		int allocateSol(const unsigned int n);
		
		//vix has the value v_i' x
		double calculateBranchingPoint(const double vix, const double l, const double u, const double alpha);
		
		//void chooseIndexToBranch(const int maxVarsToBranch, const double* sol, const double* y, double* gaps, unsigned int& numVarsToBranch, unsigned int* indices);
		
		void makeSolverInitialSolution(const int nBranchVars, const unsigned int* branchVars, const int nBounds, const iquad::IQD_NodeBounds* bounds, const int nPrimal, const double* solDad, double* initialSol);
		
		void printSolverInfo(const unsigned int nthreads);
		
		friend class IQD_Callbacks;
		
	};
    
    
	
	class IQD_Callbacks : public branchAndBound::BBL_UserCallbacks
	{
		bool origVarsUptd;
		bool useWarmStart;
		
		unsigned int nthreads;
		double hoursToSaveFile;
		int indsecconstrs;
		int mtotal;
		
		int nI; //number of integer variables
		int *intVars; //indices of integer variables
		int *yintVars; //yintVars[i] has the index of y directly related to intVars[i]. Note that can happen when v = e_k (e.g. diagonal splits)
		
		bool *constrEval;
		
		//bool **constrEvals;
		
		char *QsplitType; //zero: original Q split  one: split on Q whithout integer vars
		int **columns;
		double **values;
		
		bool yBoundsUpdated;
		double lastZuRecalcYBounds;
		double *ly, *uy; //we keep a copy of bounds to y. We can change our copy during B&B procedure...
		
		IQD_SparseMatrix *Qir; //rows of Q related to integer variables...
		
		IQD_BranchAndBound *bb;
		IQD_RefProb2 *refProb2;
		optsolvers::OPT_Solver **solvers; //we need a array of pointers because solver has different sizes in subclasses. An array of OPT_Solver would not work...
		IQD_AuxBoundsUpdate *auxBoundsUpdate;
		IQD_GeneralSolverParams *solverParams;
		
		IQD_Mutex SEMAPH_ybounds;
		
	public:
		
		IQD_Callbacks(IQD_BranchAndBound *bb, IQD_RefProb2 *refProb2, IQD_GeneralSolverParams *params);
		
		virtual ~IQD_Callbacks();
		
		
		int allocateAuxBoundsUpdtObjects( unsigned int nthreads);
		
		int allocateMemory( const unsigned int nthreads, const unsigned int mtotal);
		
		
		int allocateSolverObjects( const unsigned int nthreads);
		
		
		int allocateThreadStructures( const unsigned int nthreads, const unsigned int mtotal);
		
		
		void myChooseIndexToBranch(const double* lambda, const int maxVarsToBranch, const double* sol, const double* y, const bool largestdif, double* gaps, unsigned int& numVarsToBranch, unsigned int* indices);
		
		
		
		void desallocateAuxBoundsUpdtObjects();
		
		void desallocateMemory();
		
		void desallocateSolverObjects();
		
		
		void initialize(IQD_BranchAndBound *bb, IQD_RefProb2 *refProb2, IQD_GeneralSolverParams *params);
		
		
		int setSolverPreRelax( const unsigned int thnumber, optsolvers::OPT_Solver *solver );
		
		int setSpecialSolverParameters( optsolvers::OPT_Solver *solver );
		
		
		int updateAuxVariablesBounds(const unsigned int thnumber, optsolvers::OPT_LPSolver* solver, const double* ly, const double* uy );
		
		
		
		
		virtual int beforeAll(const unsigned int numberOfThreads, double* lx, double* ux) override;
		
		
		virtual int beforeSolvingRelaxation( const unsigned int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, bool &pruneNode) override;
		
		
		virtual int solveSubProblem(const unsigned int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES &retCode, double &objValue, double &dualObjValue, double *sol, double *dualSol, bool &generalFeasibleSol, bool &pruneNode, double &nodeLowerBound, branchAndBound::BBL_BRANCH_STRATEGY &branchStrategy) override;
		
		
		virtual int chooseIndexToBranch(const int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, branchAndBound::BBL_RETURN_CODES retCode, double objValue, double dualObjValue, double *sol, double *dualSol, unsigned int &sizeIndices, unsigned int *indices, double *breakValues1, double *breakValues2, branchAndBound::BBL_Node* &nodes) override;
		
		
		virtual int endOfIteration(const int thnumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub, branchAndBound::BBL_Node& node, const double *nlx, const double *nux) override;
		
		
		//virtual int updatingBestSolution(const int threadNumber, double* sol, double &objValue, const double ub, const long unsigned int iter) override;
		
		
		virtual void newBestSolution( const int threadNumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter ) override;
		
	};
		
		
	
	
	
    
}



#endif
