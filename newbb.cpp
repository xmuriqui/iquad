/*
* New implementation of IQD_BranchAndBound, now using BBL frameworking :)
* 
* Good luck! 
*/


#include <math.h>

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <climits>

#include <iostream>
#include <iomanip>
#include <new>

#include "BBL_branchAndBound.hpp"
#include "BBL_tools.hpp"

#include "IQD_branchAndBound.hpp"
#include "IQD_tools.hpp"


#if IQD_HAVE_MURIQUI
    #include "muriqui.hpp"
    #include "MRQ_advanced.hpp"
#endif


#if IQD_SAVE_OUTPUT_FILE_BY_HOUR
    iquad::IQD_OutFile *IQD_outFileHour = NULL;
#endif



using namespace std;
using namespace branchAndBound;
using namespace minlpproblem;
using namespace optsolvers;
using namespace iquad;






int IQD_runMuriquiHeuristic( minlpproblem::MIP_MINLPProb &prob, const int milpSolver, const int nlpSolver, IQD_GeneralSolverParams *milpSolverParams, IQD_GeneralSolverParams *nlpSolverParams, const double maxTime, const double maxCPUTime, double &objf, double *solution )
#if IQD_HAVE_MURIQUI
{
    int r;
    muriqui::MRQ_HeuristicExecutor heur;
    
    
    
    heur.setSolvers( muriqui::MRQ_intToMILPSolver(milpSolver), muriqui::MRQ_intToNLPSolver(nlpSolver) );
    
    //prob.print();
    //IQD_getchar();
    
    if( maxTime < INFINITY )
        heur.setMaxTimes( maxTime );
    
    if( maxCPUTime < INFINITY )
        heur.setMaxCPUTimes( maxCPUTime );
    
    
    
    
    r =  heur.run(prob, milpSolverParams, nlpSolverParams, INFINITY, objf, solution);
    
    return r == muriqui::MRQ_HEURISTIC_SUCCESS || r == muriqui::MRQ_OPTIMAL_SOLUTION ? 0 : IQD_UNDEFINED_ERROR;
}
#else
{
    
        
    return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif




IQD_BranchAndBound::IQD_BranchAndBound()
{
    out_best_sol = NULL;
    
    resetInputParameters();
    resetOutput();
}


IQD_BranchAndBound::~IQD_BranchAndBound()
{
    desallocateSol();
}



void IQD_BranchAndBound::resetInputParameters(void)
{
    in_build_split_without_vars_also = false;
    in_certain_of_nonconvex_Q_on_obj_f = false;
    in_considerate_integrality_at_calculating_bounds = false;
    in_priorize_obj_f_on_branch_choosing = true;
    in_recalculate_aux_bounds_on_zu_updt = false;
    in_set_secants_as_equality_constraints = true;
    in_set_special_solver_parameters = false;
    in_use_real_dif_on_branch_choosing = true;
    in_use_heuristics_on_integer_probs = true;
    
    in_sdp_solver = IQD_SS_MOSEK;
    //in_milp_solver = IQD_IS_GUROBI;
    in_qp_solver = IQD_QP_MOSEK;
    in_qcp_solver = IQD_QCP_MOSEK;
    in_nlp_solver = IQD_NLP_MOSEK;
    //in_minlp_solver = IQD_MS_BONMIN;
    
    in_split_way = IQD_SQ_SCHUR;
    in_exp_strategy = IQD_ES_BEST_LIMIT;
    in_branch_strategy = IQD_BS_MOST_VIOLATED_CONSTR;
    in_integrality_strategy = IQD_IS_ON_BB_ENUMERATION;
    in_aux_bounds_calc_strategy = IQD_ABCS_SUBPROBLEMS;
    
    in_relax_prob_writing_mode = IQD_RPWM_NO_WRITING;
    
    in_print_level = 2;
    
    //in_eps_added_sd_solver_matrices = 1.0e-3;
    //in_eps_added_to_P_non_sdp = 0.0;
    in_absolute_eps_added_to_diag_P = 1.0e-6;
    in_relative_eps_added_to_diag_P = 0.0;
    in_alpha = 0.8;
    
    in_factor_to_prior_obj_f_on_branch = 0.01;
    in_factor_to_recalculate_aux_bounds_on_zu_updt = 0.1;
    
    in_abs_feas_tol = 1.0e-6;
    in_rel_feas_tol = 1.0e-6;
    
    in_infinity_value = IQD_INFINITY;
    in_integer_tol = 1.0e-4;
    
    in_lower_bound = -INFINITY;
    in_upper_bound = INFINITY;
    in_absolute_convergence_tol = 1.0e-4;
    in_relative_convergence_tol = 1.0e-3;
    
    in_zero_tol_to_split = 1.0e-6;
    
    in_alpha_to_integer_branch_choosing = 0.5;
    
    in_number_of_threads = 1;
    in_printing_frequency = 100;
    in_number_of_node_sublists = 100;
    
    in_reorganize_lists = true;
    in_lists_reorganization_frequency = 10000;
    
    in_max_iters = ULONG_MAX;
    
    in_max_cpu_time = INFINITY;
    in_max_time = INFINITY;
    
    in_hours_to_save_file = 0.5;
    
    #if IQD_MURIQUI_GAP_MIN_TREATMENT
        mrq_gap_min = false;
    #endif
}


void IQD_BranchAndBound::resetOutput(void)
{
    out_number_of_iterations = 0;
    out_return_code = IQD_UNDEFINED;
    out_return_subcode = IQD_UNDEFINED;
    
    desallocateSol();
    
    out_feasible_sol = false;
    out_best_sol = NULL;
    out_best_obj = IQD_INFINITY;
    out_lower_bound = out_root_node_lower_bound = -IQD_INFINITY;
    out_upper_bound = out_root_node_upper_bound = IQD_INFINITY;
    out_number_of_open_nodes = 0;
    out_cpu_time = out_clock_time = out_splitting_cpu_time = -1.0;
}


int IQD_BranchAndBound::allocateSol( const unsigned int n )
{
    desallocateSol();
    out_best_sol = (double *) malloc( n * sizeof(double) );
    
    if( !out_best_sol )
        return IQD_MEMORY_ERROR;
    
    for(unsigned int i = 0; i < n; i++)
        out_best_sol[i] = NAN;
    
    return 0;
}


void IQD_BranchAndBound::desallocateSol(void)
{
    IQD_secFree( out_best_sol );
}




int IQD_BranchAndBound::setIntegerParameter(const char *name, const long int value)
{
    int ret = 0;
    
    
    if( IQD_setAtt<bool>( IQD_STRATT(in_build_split_without_vars_also), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_certain_of_nonconvex_Q_on_obj_f), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_considerate_integrality_at_calculating_bounds), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_reorganize_lists), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_set_secants_as_equality_constraints), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_set_special_solver_parameters), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_priorize_obj_f_on_branch_choosing), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_recalculate_aux_bounds_on_zu_updt), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_use_real_dif_on_branch_choosing), name, value )  == 0  )
    {}
    else if( IQD_setAtt<bool>( IQD_STRATT(in_use_heuristics_on_integer_probs), name, value )  == 0  )
    {}
    
    else if( IQD_setAtt<unsigned int>( IQD_STRATT(in_lists_reorganization_frequency), name, value )  == 0  )
    {}
    else if( IQD_setAtt<unsigned int>( IQD_STRATT(in_number_of_node_sublists), name, value )  == 0  )
    {}
    else if( IQD_setAtt<unsigned int>( IQD_STRATT(in_number_of_threads), name, value )  == 0  )
    {}
    else if( IQD_setAtt<unsigned int>( IQD_STRATT(in_printing_frequency), name, value )  == 0  )
    {}
    else if( IQD_setAtt<unsigned int>( IQD_STRATT(in_print_level), name, value )  == 0  )
    {}
    else if( IQD_setAtt<long unsigned int>( IQD_STRATT(in_max_iters), name, value )  == 0  )
    {}
    else
        ret = IQD_NAME_ERROR;
    
    
    return ret;
}




int IQD_BranchAndBound::setDoubleParameter(const char *name, const double value)
{
    int ret = 0;
    
    
    if( IQD_setAtt<double>( IQD_STRATT(in_absolute_convergence_tol), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_alpha), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_alpha_to_integer_branch_choosing), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_absolute_eps_added_to_diag_P), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_relative_eps_added_to_diag_P), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_factor_to_prior_obj_f_on_branch), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_factor_to_recalculate_aux_bounds_on_zu_updt), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_abs_feas_tol), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_rel_feas_tol), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_infinity_value), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_integer_tol), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_lower_bound), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_max_time), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_max_cpu_time), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_upper_bound), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_relative_convergence_tol), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_zero_tol_to_split), name, value ) == 0 )
    {}
    else if( IQD_setAtt<double>( IQD_STRATT(in_hours_to_save_file), name, value ) == 0 )
    {}
    else
        ret = IQD_NAME_ERROR;
    
    
    return ret;
}



int IQD_BranchAndBound::setStringParameter(const char *name, const char *value)
{
    int r = 0, ret;
    
    
        
    if(  (r = IQD_setStrAtt( IQD_STRATT(in_sdp_solver), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_qp_solver), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_qcp_solver), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_nlp_solver), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_split_way), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_exp_strategy), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_branch_strategy), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_integrality_strategy), name, value ) ) >= 0  ||
    (r = IQD_setStrAtt( IQD_STRATT(in_aux_bounds_calc_strategy), name, value ) ) >= 0 ||
    (r = IQD_setStrAtt( IQD_STRATT(in_relax_prob_writing_mode), name, value ) )
    )
    {
        ret = r == 0 ? 0 : IQD_VALUE_ERROR;
    }
    else
        ret = IQD_NAME_ERROR;
    
    
    return ret;
}



double IQD_BranchAndBound::calculateBranchingPoint(const double vix, const double l, const double u, const double alpha)
{
    //printf("Calculating calculateBranchingPoint. vix: %f l: %f u: %f alpha: %0.2f\n", vix, l, u, alpha);
    return alpha*vix + (1.0 - alpha)*(l + u)/2.0;
}



void IQD_BranchAndBound::makeSolverInitialSolution(const int nBranchVars, const unsigned int* branchVars, const int nBounds, const iquad::IQD_NodeBounds* bounds, const int nPrimal, const double* solDad, double* initialSol)
{
    int i, j, k;
    
    IQD_copyArray(nPrimal, solDad, initialSol);
    
    for(k = 0; k < nBranchVars; k++)
    {
        for(i = 0; i < nBounds; i++)
        {
            if( bounds[i].ind == branchVars[k])
            {
                j = branchVars[k];
                
                if( initialSol[j] < bounds[i].l )
                    initialSol[j] = bounds[i].l;
                else if( initialSol[j] > bounds[i].u )
                    initialSol[j] = bounds[i].u;
                //if initial solution is already in the bounds, we let in this way
                break;
            }
            else if( bounds[i].ind > branchVars[k] )
            {
                break; //bounds are sorted by index
            }
        }
    }
}





void IQD_BranchAndBound::printSolverInfo(const unsigned int nthreads)
{
    char sdp[20], qp[20], miqcp[20], nlp[20];
    
    if(in_sdp_solver == IQD_SS_MOSEK)
        sprintf(sdp, "mosek");
    else if(in_sdp_solver == IQD_SS_CSDP)
        sprintf(sdp, "csdp");
    else
        sprintf(sdp, "invalid");
    
    
    if( in_qp_solver == IQD_QP_MOSEK )
        sprintf(qp, "mosek");
    else if( in_qp_solver == IQD_QP_XPRESS )
        sprintf(qp, "xpress");
    else if( in_qp_solver == IQD_QP_CPLEX )
        sprintf(qp, "cplex");
    else if( in_qp_solver == IQD_QP_GUROBI )
        sprintf(qp, "gurobi");
    else if( in_qp_solver == IQD_QP_IPOPT )
        sprintf(qp, "ipopt");
    else
        sprintf(qp, "invalid");
    
    
    if( in_qcp_solver == IQD_QCP_MOSEK )
        sprintf(miqcp, "mosek");
    else if( in_qcp_solver == IQD_QCP_XPRESS )
        sprintf(miqcp, "xpress");
    else if( in_qcp_solver == IQD_QCP_IPOPT )
        sprintf(miqcp, "ipopt");
    else
        sprintf(miqcp, "invalid");
    
    if(in_nlp_solver == IQD_NLP_MOSEK)
        sprintf(nlp, "mosek");
    else if(in_nlp_solver == IQD_NLP_IPOPT)
        sprintf(nlp, "ipopt");
    
    
    //printf("Threads:%d\nSDP Solver: %s\nMIQCP Solver: %s\nNLP Solver: %s\n", nthreads, sdp, miqcp, nlp);
    
    std::cout << "Threads: " << nthreads << "\nSDP Solver: " << sdp << "\nQP Solver: " << qp << "\nQCP Solver: " << miqcp << "\nNLP Solver: " << nlp << "\n";
}






int IQD_BranchAndBound::run(IQD_IQuadProb& prob, IQD_GeneralSolverParams* solverParams, IQD_GeneralSolverParams* sdpParams)
{
    timeStart = IQD_getTime();
    clockStart = clock();
    
    
    const IQD_PROBLEMTYPE probType = prob.getProblemType();
    
    unsigned int nthreads, ndual;
    int ret;
    
    double *lxy = NULL, *uxy = NULL;
    
    
    clock_t clockAux;
    
    BBL_BranchAndBound bbcore;
    
    IQD_Callbacks *callbakcs = NULL;
    //IQD_RefProb refProb;
    IQD_RefProb2 refProb;
    
    
    IQD_welcome();
    
    
    
    
    nthreads = in_number_of_threads == 0 ? BBL_getNumCores() : in_number_of_threads;
    
    
    
    if(in_print_level > 1)
    {
        std::cout << IQD_PREPRINT "Problem type: ";
        
        if( probType == MIP_PT_QP )
            std::cout << "QP";
        else if( probType == MIP_PT_MIQP )
            std::cout << "MIQP";
        else if( probType == MIP_PT_QCP )
            std::cout << "QCP";
        else if( probType == MIP_PT_MIQCP )
            std::cout << "MIQCP";
        else if( probType == MIP_PT_NLP )
            std::cout << "NLP";
        else if( probType == MIP_PT_MINLP )
            std::cout << "MINLP";
        
        std::cout << "\n";
    }
    
    
    
    if(in_print_level > 1)
        printSolverInfo( nthreads );
    
    
    if( prob.nlEval )
    {
        ret = prob.nlEval->initialize(nthreads, prob.n, prob.m, prob.J.getNumberOfElements(),  prob.lagH.getNumberOfElements());
        
        
        if(ret != 0)
        {
            if( in_print_level >= 1 )
                IQD_PRINTMEMERROR;
            
            out_return_code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
    }
    
    ret = allocateSol(prob.n);
    if( ret != 0 )
    {
        if( in_print_level >= 1 )
            IQD_PRINTMEMERROR;
        
        out_return_code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    if(in_print_level > 1)
        std::cout << IQD_PREPRINT "Spliting matrices\n";
    
    
    //refProb.zero_tol = in_zero_tol_to_split;
    //refProb.eps_diag_matrix_from_sdp = in_eps_added_sd_solver_matrices;
    //refProb.setWayOfSplitQ(in_split_way);
    //refProb.prob = &prob;
    
    refProb.zero_tol = in_zero_tol_to_split;
    refProb.abs_eps_diag_P = in_absolute_eps_added_to_diag_P;
    refProb.rel_eps_diag_P = in_relative_eps_added_to_diag_P;
    //refProb.eps_diag_matrix_from_sdp = in_eps_added_sd_solver_matrices; //do not add this..
    refProb.splitWay = in_split_way;
    refProb.max_cpu_time = in_max_cpu_time;
    refProb.max_time = in_max_time;
    refProb.oprob = &prob;
    
    clockAux = clock();
    
    //ret = refProb.splitQ( in_sdp_solver, 1986, in_qcp_solver, in_nlp_solver, 1986, sdpParams, NULL, solverParams );
    
    
    
    ret = refProb.buildRefProb( in_build_split_without_vars_also,  in_certain_of_nonconvex_Q_on_obj_f, in_sdp_solver, in_qp_solver, in_qcp_solver, in_nlp_solver, sdpParams, solverParams );
    
    
    
    out_splitting_cpu_time =  ((double) (clock() - clockAux))/CLOCKS_PER_SEC;
    
    if( ret != 0 )
    {
        if( in_print_level > 1 )
        {
            if( ret == IQD_NO_BOUNDED_FEASIBLE_SET )
                IQD_PRINTERRORMSG("Error on calculating bounds on directions! DIrection is unbounded! Try set lower and upper bounds to your variables!");
            else
                std::cerr << IQD_PREPRINT << "Error " << ret << " at splitting matrix Q" << IQD_GETFILELINE << endl ;
        }
        
        
        out_return_code = ret == IQD_NO_BOUNDED_FEASIBLE_SET ? IQD_NO_BOUNDED_FEASIBLE_SET : IQD_SPLITING_ERROR;
        out_return_subcode = ret;
        
        goto termination;
    }
    
    /*{
        
        cout << "#####################################################################################################" << endl;
        
        cout << "v':" << endl;
        refProb.v.printAllMatrix();
        
        cout << endl << "lambdaObj: ";
        for(int j = 0; j < refProb.dimv; j++)
            cout << refProb.lambdaObj[j] << " " ;
        cout << endl;
        
        cout << "Q: " << endl;
        prob.Q.printAllMatrix();
        
        cout << "P: " << endl;
        refProb.rprob.Q.printAllMatrix();
        
        
        for(int i = 0; i < refProb.rprob.m; i++ )
        {
            if( refProb.rprob.QC[i].getNumberOfElements() > 0 )
            {
                cout << endl << "lambdaConstr " << i << ": ";
                for(int j = 0; j < refProb.dimv; j++)
                    cout << refProb.lambdaConstr[i][j] << " " ;
                cout << endl;
                
                cout << "QC " << i << ": " ;
                if( i >= prob.m )
                    cout << "(-1.0)*" ;
                
                cout << endl;
                prob.QC[i < prob.m ? i : refProb.origConstr[i]].printAllMatrix();
                
                cout << "P: " << endl;
                refProb.rprob.QC[i].printAllMatrix();
            }
        }
        
        refProb.rprob.writeMIQCPModelInAMPLFile( "refProb.mod" );
        
        getchar();
    } */
    
    
    
    
    /*
    printf("Q: \n");
    prob.Q.printAllMatrix();
    
    if(refProb->v)
    {
        printf("v: \n");
        for(i = 0; i < prob.n; i++)
        {
            for(j = 0; j < prob.n; j++)
                printf("%0.16f, ", refProb->v[i][j]);
            
            printf("\n");
        }
    }
    
    printf("lambda:\n");
    for(i = 0; i < n; i++)
        printf("%0.16f ;", refProb->lambda[i]);
    printf("\n");
    
    printf("P: \n");
    refProb->P->printAllMatrix();
    
    if(refProb->R)
    {
        printf("R: \n");
        for(i = 0; i < prob.n; i++)
        {
            for(j = 0; j < prob.n; j++)
            {
                printf("%f ", refProb->R[ i*n + j ] );
            }
            printf(";\n");
        }
    }
    */ 
    /*
    printf("A: \n");
    refProb->prob->qc->A->printAllMatrix();
    
    printf("b:\n");
    for(i = 0; i < prob.mqc; i++)
        printf("%0.16f ,", refProb->prob->qc->b[i]);
    
    printf("\nc: \n");
    for(i = 0; i < prob.n; i++)
        printf("%0.16f ,", prob.c[i]);
    
    printf("-------------------------------------------------------------------\n");
    getchar();
    */
    
    
    ndual =  refProb.getNDual(); //2*refProb2.rprob.m + 2*refProb2.rprob.n + 4*refProb2.dimv;  //  
    
    
    
    //if(in_print_level > 2)
        //cout << IQD_PREPRINT << "Using solver " << solver->solverName << " to compute relaxations." << endl;
    
    
    if(in_print_level > 4)
    {
        //const unsigned int n = prob.n;
        
        std::cout << IQD_PREPRINT;
        
        for(int w = 0; w < refProb.dimv; w++)
            std::cout << "ly[" << w << "]: " << refProb.ly[w] << " uy[" << w << "]: " << refProb.uy[w] << "\t";
        
        std::cout << std::endl;
    }
    
    
    if( IQD_getTime()-timeStart > in_max_time )
    {
        out_return_code = IQD_MAX_TIME_STOP;
        goto termination;
    }
    
    if( ( (double) (clock() - clockStart))/CLOCKS_PER_SEC > in_max_cpu_time )
    {
        out_return_code = IQD_MAX_TIME_STOP;
        goto termination;
    }
    
    
    
    lxy = (double *) malloc( 2*(prob.n + refProb.dimv) * sizeof(double) );
    
    callbakcs = new (nothrow) IQD_Callbacks(this, &refProb, solverParams);
    
    if( !lxy || !callbakcs )
    {
        if(in_print_level > 0)
            std::cerr << IQD_PREPRINT << "Memory error" << IQD_GETFILELINE << endl;
        
        out_return_code = iquad::IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    #if IQD_SAVE_OUTPUT_FILE
        if( IQD_outFileHour )
            IQD_outFileHour->lastHoursToSaveFile = 0.0;
    #endif
    
    
    
    //note that the first n indices in lx and ux has the bounds for the original variables. The next indices are bounds for y, not to w... We do it because is under y space we perform the spatial branching...
    
    uxy = &lxy[ prob.n + refProb.dimv ];
    
    IQD_copyArray( prob.n, prob.lx, lxy );
    IQD_copyArray( prob.n, prob.ux, uxy );
    
    IQD_copyArray( refProb.dimv, refProb.ly, &lxy[prob.n] );
    IQD_copyArray( refProb.dimv, refProb.uy, &uxy[prob.n] );
    
    
    useWarmStart = (probType == MIP_PT_NLP || probType == MIP_PT_MINLP);
    
    bbcore.in_call_updating_best_solution_callback = false;
    bbcore.in_call_new_best_solution_callback = in_recalculate_aux_bounds_on_zu_updt;
    
    bbcore.in_use_dual_obj_to_bound_prunning = true;
    bbcore.in_consider_relax_infeas_if_solver_fail = true;
    
    bbcore.in_reorganize_lists = in_reorganize_lists;
    bbcore.in_store_parent_primal_solution_on_nodes = bbcore.in_store_parent_dual_solution_on_nodes = useWarmStart;
    
    bbcore.in_exp_strategy = BBL_int2extExpStrategy( in_exp_strategy );
    
    
    bbcore.in_lists_reorganization_frequency = in_lists_reorganization_frequency;
    bbcore.in_number_of_node_sublists = in_number_of_node_sublists;
    bbcore.in_printing_frequency = in_printing_frequency;
    bbcore.in_print_level = in_print_level;
    
    
    bbcore.in_number_of_threads = nthreads;
    
    
    bbcore.in_max_cpu_time = in_max_cpu_time -  ( (double) (clock() - clockStart))/CLOCKS_PER_SEC;
    bbcore.in_max_time = in_max_time - (IQD_getTime() - timeStart );
    bbcore.in_max_iterations = in_max_iters;
    
    
    bbcore.in_absolute_convergence_tol = in_absolute_convergence_tol;
    
    bbcore.in_lower_bound = in_lower_bound;
    bbcore.in_upper_bound = in_upper_bound;
    bbcore.in_infinity_value = in_infinity_value;
    
    bbcore.in_relative_convergence_tol = in_relative_convergence_tol;
    
    
    bbcore.run(prob.n + refProb.dimv, refProb.rprob.m + 2*refProb.dimv, ndual, lxy, uxy, *callbakcs );
    
    
    out_feasible_sol = bbcore.out_feasible_sol;
    
    if( out_feasible_sol )
    {
        IQD_copyArray( prob.n, bbcore.out_best_sol, out_best_sol );
        
        IQD_copyArray( prob.n, out_best_sol, prob.x );
        
        prob.objValue = out_best_obj = bbcore.out_best_obj;
    }
    
    out_number_of_iterations = bbcore.out_number_of_iterations;
    
    out_cpu_time_to_fisrt_feas_sol = bbcore.out_cpu_time_to_fisrt_sol;
    out_number_of_feas_sols = bbcore.out_number_of_feas_sols;
    
    out_lower_bound = bbcore.out_lower_bound;
    out_upper_bound = bbcore.out_upper_bound;
    out_root_node_lower_bound = bbcore.out_obj_opt_at_root_relax;
    
    out_number_of_open_nodes = bbcore.out_number_of_open_nodes;
    
    out_return_subcode = bbcore.out_return_subcode == BBL_UNDEFINED ? IQD_UNDEFINED : bbcore.out_return_subcode;
    
    switch( bbcore.out_return_code )
    {
        case BBL_OPTIMAL_SOLUTION:
            out_return_code = IQD_OPTIMAL_SOLUTION;
            break;
        case BBL_INFEASIBLE_PROBLEM:
        case BBL_NO_FEASIBLE_SOLUTION_FOUND:
            out_return_code = IQD_INFEASIBLE_PROBLEM;
            break;
        case BBL_MAX_TIME_STOP:
            out_return_code = IQD_MAX_TIME_STOP;
            break;
        case BBL_MAX_ITERATIONS_STOP:
            out_return_code = IQD_MAX_ITERATIONS_STOP;
            break;
        default:
            out_return_code = IQD_UNDEFINED_ERROR;
            break;
    }
    
    
    out_cpu_time = ( (double) (clock() - clockStart))/CLOCKS_PER_SEC;
    out_clock_time = IQD_getTime() - timeStart;
    
    
    
    
    
    
    
termination:
    
    #if IQD_SAVE_OUTPUT_FILE_BY_HOUR
        if( IQD_outFileHour )
        {
            IQD_outFileHour->writeToCompleteHours( 5, in_max_cpu_time, in_hours_to_save_file );
            IQD_outFileHour->flush();
        }
    #endif
    
    
    if(callbakcs)	delete callbakcs;
    
    if(lxy)			free(lxy);
    
    
    
    return out_return_code;
}






IQD_Callbacks::IQD_Callbacks(IQD_BranchAndBound *bb, IQD_RefProb2 *refProb2, IQD_GeneralSolverParams *params):BBL_UserCallbacks()
{
    initialize( bb, refProb2, params );
}


IQD_Callbacks::~IQD_Callbacks()
{
    desallocateMemory();
    desallocateSolverObjects();
    desallocateAuxBoundsUpdtObjects();
}


int IQD_Callbacks:: allocateAuxBoundsUpdtObjects( unsigned int nthreads)
{
    auxBoundsUpdate = new (nothrow) IQD_AuxBoundsUpdate[nthreads];
    
    return auxBoundsUpdate ? 0 : IQD_MEMORY_ERROR;
}



int IQD_Callbacks::allocateThreadStructures( const unsigned int nthreads, const unsigned int mtotal)
{
    const int dimv = refProb2->dimv;
    const unsigned int nprimal = refProb2->rprob.n ;
    const unsigned int size2 = IQD_max<unsigned int>( nprimal, dimv );
    
    unsigned int size =  IQD_max( nprimal + (bb->useWarmStart ? refProb2->getNDual() : 0) , 2u* dimv) ;  //we need nprimal + ndual if we do warm start, and we need 2*dimv to calculate gap and y on branching...
    
    
    
    size = IQD_max( size, mtotal );
    
    QsplitType = (char *) calloc( nthreads, sizeof(char) );
    values = (double **) calloc( nthreads, sizeof(double*) );
    columns = (int **) calloc( nthreads, sizeof(int*) );
    if( !QsplitType || !values || !columns )
    {
        if( bb->in_print_level > 0 )
            IQD_PRINTMEMERROR;
        return IQD_MEMORY_ERROR;
    }
    
    
    for( unsigned int i = 0; i < nthreads; i++ )
    {
        values[i] = (double *) malloc( size * sizeof(double) );
        columns[i] = (int *) malloc( size2 * sizeof(int) );
        
        if( !values[i] || !columns[i] )
        {
            if( bb->in_print_level > 0 )
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
    }
    
    
    /*values[0] = (double *) malloc(  nthreads*size * sizeof(double)  );
    columns[0] = (int *) malloc( nthreads*size2 * sizeof(int)  );
    
    if( !values[0] || !columns[0] )
    {
        if( bb->in_print_level > 0 )
            IQD_PRINTMEMERROR;
        return IQD_MEMORY_ERROR;
    }
    
    
    for(unsigned int i = 1; i < nthreads; i++)
        values[i] = &(values[0][ i*size ]);
    
    for(unsigned int i = 1; i < nthreads; i++)
        columns[i] = &(columns[0][ i*size2 ]); */
    
    
    return 0;
}


int IQD_Callbacks::allocateMemory( const unsigned int nthreads, const unsigned int mtotal)
{
    const int dimv = refProb2->dimv;
    const unsigned int nprimal = refProb2->rprob.n ;
    //const unsigned int size2 = IQD_max<unsigned int>( nprimal, dimv );
    const IQD_SparseMatrix &v = refProb2->v;
    
    
    int r;
    
    unsigned int size =  IQD_max( nprimal + (bb->useWarmStart ? refProb2->getNDual() : 0) , 2u* dimv) ;  //we need nprimal + ndual if we do warm start, and we need 2*dimv to calculate gap and y on branching...
    
    
    size = IQD_max(size, mtotal);
    
    
    
    /*values = (double **) calloc( nthreads, sizeof(double*) );
    columns = (int **) calloc( nthreads, sizeof(int*) );
    if( !values || !columns )
    {
        if( bb->in_print_level > 0 )
            IQD_PRINTMEMERROR;
        return IQD_MEMORY_ERROR;
    }
    
    
    
    values[0] = (double *) malloc(  nthreads*size * sizeof(double)  );
    columns[0] = (int *) malloc( nthreads*size2 * sizeof(int)  );
    
    if( !values[0] || !columns[0] )
    {
        if( bb->in_print_level > 0 )
            IQD_PRINTMEMERROR;
        return IQD_MEMORY_ERROR;
    }
    
    
    for(unsigned int i = 1; i < nthreads; i++)
        values[i] = &(values[0][ i*size ]);
    
    for(unsigned int i = 1; i < nthreads; i++)
        columns[i] = &(columns[0][ i*size2 ]); */
    
    
    r = allocateThreadStructures(nthreads, mtotal);
    
    if( r != 0 )
    {
        if( bb->in_print_level > 0 )
            IQD_PRINTERRORNUMBER(r);
        return IQD_MEMORY_ERROR;
    }
    
    
    
    if( mtotal > 0 )
    {
        const int m = refProb2->oprob->m;
        const int mdis = refProb2->mdis;
        const int *origConstr = refProb2->origConstr;
        //const double *lc = refProb2->oprob->lc;
        //const double *uc = refProb2->oprob->uc;
        double **lambdaConstr = refProb2->lambdaConstr;
        
        constrEval = (bool *) malloc( mtotal * sizeof(bool) );
        if( !constrEval )
        {
            if( bb->in_print_level > 0 )
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
        
        for(int i = 0; i < m; i++)
        {
            constrEval[i] = lambdaConstr[i] != NULL; //if lambdaConstr[i] is NULL, is because constraint i is convex and we will not need evaluate... unless constraint is double bounded, but we treat this case in the next for...
        }
            
        
        for(int i = 0; i < mdis; i++)
        {
            constrEval[ origConstr[m+i] ] = true;
        }
        
        
        for(unsigned int i = m; i < mtotal; i++)
            constrEval[i] = false;
        
        
        /*constrEvals = (bool **) calloc( nthreads, sizeof(bool *) );
        if( !constrEvals )
        {
            if( bb->in_print_level > 0 )
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
    
        
        constrEvals[0] = (bool *) malloc( nthreads*mtotal * sizeof(bool) );
        
        if(  !constrEvals[0] )
        {
            if( bb->in_print_level > 0 )
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
        
        for(unsigned int i = 1; i < nthreads; i++)
            constrEvals[i] = &( constrEvals[0][i*mtotal] );
        
        bool *auxb = constrEvals[0];
        
        for(unsigned int i = 0; i < nthreads*mtotal; i++)
            auxb[i] = true; //actually we just need to set true for nonconvex constraints...
        */
    }
    
    
    ly = (double *) malloc( 2*dimv * sizeof(double) );
    
    if( !ly )
    {
        if( bb->in_print_level > 0 )
            IQD_PRINTMEMERROR;
        return IQD_MEMORY_ERROR;
    }
    
    uy = &ly[dimv];
    
    
    IQD_copyArray( dimv, refProb2->ly, ly );
    IQD_copyArray( dimv, refProb2->uy, uy );
    
    
    nI = refProb2->oprob->getNumberOfIntegerVars();
    if( nI > 0 )
    {
        intVars = (int *) malloc( 2*nI * sizeof(int) );
        if( !intVars )
        {
            if( bb->in_print_level > 0 )
                IQD_PRINTMEMERROR;
                
            return IQD_MEMORY_ERROR;
        }
        
        yintVars= &intVars[nI];
        refProb2->oprob->getIntegerIndices(intVars);
        
        IQD_setAllarray(nI, yintVars, -1);
        
        
        int *cols = columns[0];
        IQD_setAllarray( nprimal, cols, -1 );
        
        for(int i = 0; i < dimv; i++)
        {
            //IQD_SparseRow &rv = refProb2->v[i];
            
            if( v.getNumberOfElementsAtRow(i) == 1 ) //if( rv.getNumberOfElements() == 1 )
            {
                #if IQD_DEBUG_MODE
                    assert( IQD_abs( v(i)[0]-1.0) <= 0.00001 ); //assert( IQD_abs(rv[0].getValue()-1.0) <= 0.00001 ); //checking if the value is 1.0. So the original and artifical varable are directly linked
                #endif
                
                #if IQD_DEBUG_MODE
                    assert( cols[ v[i][0] ] == -1 ); //assert( cols[ rv[0].getColumn() ] == -1 ); //if this test fail, is because we have two lines equals in v. Try fix it. otherwise, yintVars[i] should have store two or ore indices
                #endif
                
                cols[ v[i][0] ] = i; //cols[ rv[0].getColumn() ] = i; //variable in rv[0] is related to variable i
            }
        }
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( cols[ ind ] >= 0 )
                yintVars[i] = cols[ind];
        }
        
        
        /*std::cout << "v: " << std::endl;
        refProb2->v.printSparseMatrix();
        for(int i = 0; i < nI; i++)
        {
            std::cout << "intVars[" << i << "]: " << intVars[i] << " yintVars["<<i<<"]: " << yintVars[i] << std::endl;
        } */
    }
    
    
    return 0;
}


int IQD_Callbacks::allocateSolverObjects( const unsigned int nthreads)
{
    const IQD_PROBLEMTYPE probType = refProb2->rprob.getProblemType();
    int code;
    
    
    
    solvers = (OPT_Solver **) calloc( nthreads, sizeof(OPT_Solver*) );
    
    if(!solvers)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    
    if( MIP_isQuadOnlyProblemType( probType ) )
    {
        for( unsigned int i = 0; i < nthreads; i++ )
        {
            solvers[i] = OPT_newQPSolver( bb->in_qp_solver );
            
            if( !solvers[i] )
            {
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
        }
        
    }
    else if( MIP_isQuadConstrProblemType(probType) )
    {
        for( unsigned int i = 0; i < nthreads; i++ )
        {
            solvers[i] = OPT_newQCPSolver( bb->in_qcp_solver );
            
            if( !solvers[i] )
            {
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
        }
    }
    else //if(MIP_isNonlinearProblemType(probType) )
    {
        for( unsigned int i = 0; i < nthreads; i++ )
        {
            solvers[i] = OPT_newNLPSolver( bb->in_nlp_solver );
            
            if( !solvers[i] )
            {
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            
            ( (OPT_NLPSolver *) solvers[i] )->setThreadNumber(i);
        }
        
    }
    
    
    
    code = 0;
    
termination:
    
    if( code != 0 )
    {
        desallocateSolverObjects();
    }
    

    return code;
}


static inline double IQD_calDif( const double lambda, const double y, const double w )
{
    return lambda *( y*y + w );
}



void IQD_Callbacks:: myChooseIndexToBranch( const double *lambda, const int maxVarsToBranch, const double* sol, const double* y, const bool largestdif, double* gaps, unsigned int& numVarsToBranch, unsigned int* indices)
{
    const int n = refProb2->oprob->n;
    const int dimv = refProb2->dimv;
    //const double *lambda = rowIndex < 0 ? refProb2->lambdaObj : refProb2->lambdaConstr[rowIndex];
    const double zero_tol = refProb2->zero_tol;
    //int i;
    
    
    //by now, we choose the variable with the highest gap beetwen y and w.
    
    /*if( IQD_isDiagonalSplitQWay(bb->in_split_way) )
    {
        //we use the diference beetwen w_i and x_i²
        
        for(i = 0; i < dimv; i++)
        {
            if( lambda[i] >= -zero_tol )
            {
                gaps[i] = 0.0;
                continue;
            }
            
            gaps[i] = -sol[i]*sol[i] - sol[n+i];
            #if IQD_DEBUG_MODE
                if(gaps[i] < - 100*sqrt(zero_tol) )
                {
                    printf("Negative gap at IQD_BranchAndBound::chooseIndexToBranch. gaps[%d]: %e lambda: %f\n", i, gaps[i], lambda[i]);
                    getchar();
                }
            #endif
        }
    }
    else */
    {
        //we use the diference beetwen w_i and -(v_i'x)²
        for(int i = 0; i < dimv; i++)
        {
            if( lambda[i] >= -zero_tol )
            {
                gaps[i] = 0.0;
                continue;
            }
            
            gaps[i] = -y[i]*y[i] - sol[n+i];
            
            #if IQD_DEBUG_MODE
                if(gaps[i] < - 100*sqrt(zero_tol))
                {
                    printf("Negative gap at IQD_BranchAndBound::chooseIndexToBranch. gaps[%d]: %e lambda: %0.10f zero_tol: %0.10f dimv: %d\n", i, gaps[i], lambda[i], zero_tol, dimv);
                    
                    auto nzvi = refProb2->v.getNumberOfElementsAtRow(i);
                    int *cols = refProb2->v[i];
                    double *vals = refProb2->v(i);
                    
                    cout << "v: ";
                    for(decltype(nzvi) j = 0; j < nzvi; j++ )
                    {
                        cout << cols[j] << ": " << vals[j] << "  ";
                    }
                    cout << "\n";
                    
                    for(int j = 0; j < dimv; j++)
                        cout << "lambda[" << j << "]: " << lambda[j] << " sol[" << n+j << "]: " << sol[n+j] << " y[" << j << "]: " << y[j] << " gaps[" << j << "]: " << gaps[j] << endl;
                    
                    if( nthreads == 1 )
                        IQD_getchar();
                }
            #endif
        }
    }
    
    
    //for(int i = 0; i < dimv; i++)
        //cout << fixed << setprecision(8) <<  "lambda[" << i << "]: " << lambda[i] << " gaps["<< i << "]: " << gaps[i] << " \t";// endl;
    //std::cout << endl;
    
    
    if(largestdif)
    {
        for(int i = 0; i < dimv; i++)
            gaps[i] = -lambda[i] *gaps[i];
    } 
    
    
    indices[0] = 0; //we start from the first index
    
    for(int i = 1; i < dimv; i++)
    {
        if( gaps[i] > gaps[ indices[0] ] )
        {
            indices[0] = i;
        }
    }
    
    numVarsToBranch = 1;
    
    
    
    
    #if IQD_DEBUG_MODE
        for(unsigned int w = 0; w < numVarsToBranch; w++)
        {
            if( gaps[ indices[w] ] <= 0.0 )
            {
                
                for(int i = 0; i < dimv; i++)
                    cout << fixed << setprecision(8) <<  "lambda[" << i << "]: " << lambda[i] << " gaps["<< i << "]: " << gaps[i] << " \t\n";
                std::cout << "\n";
                
                
                std::cout << "constraint lambda: \n";
                for(int i = 0; i < refProb2->oprob->m; i++)
                {
                    std::cout << i << ":  " ;
                    
                    const  double *l = refProb2->lambdaConstr[i];
                    
                    if( l )
                    {
                        for(int j = 0; j < dimv; j++)
                            cout << l[j] << " ";
                    }
                    std::cout << "\n";
                }
                
                
                std::cout << "gap nao positivo em indices["<<w<<"]: " << indices[w] << " gap: " << gaps[ indices[w] ] << "\n";
                assert( false );
            }
        }
    #endif
    
    
    //return 0;
}




void IQD_Callbacks:: desallocateAuxBoundsUpdtObjects( )
{
    IQD_secDeleteArray(auxBoundsUpdate);
}



void IQD_Callbacks:: desallocateMemory( )
{
    if( values )
    {
        if( values[0] )
            free( values[0] );
        
        free(values);
        values = NULL;
    }
    
    if( columns )
    {
        if( columns[0] )
            free( columns[0] );
        
        free( columns );
        columns = NULL;
    }
    
    /*if( constrEvals )
    {
        if( constrEvals[0] )
            free( constrEvals[0] );
        
        free(constrEvals);
        constrEvals = NULL;
    } */
    
    IQD_secFree(constrEval);
    
    IQD_secFree(ly);
    
    IQD_secFree(intVars);
    
    IQD_secFree(QsplitType);
    
    IQD_secDelete(Qir);
}



void IQD_Callbacks:: desallocateSolverObjects( )
{
    if( solvers )
    {
        for(unsigned int i = 0; i < nthreads; i++)
        {
            if( solvers[i] )
                delete solvers[i];
        }
        
        free(solvers);	solvers = NULL;
        nthreads = 0;
    }
    
}



void IQD_Callbacks::initialize(IQD_BranchAndBound *bb, IQD_RefProb2 *refProb2, IQD_GeneralSolverParams *params)
{
    this->bb = bb;
    this->refProb2 = refProb2;
    this->solverParams = params;
    
    nthreads = 0;
    
    ly = uy = NULL;
    
    Qir = NULL;
    
    constrEval = NULL;
    //constrEvals = NULL;
    columns = NULL;
    values = NULL;
    solvers = NULL;
    auxBoundsUpdate = NULL;
    
    nI = 0;
    intVars = NULL;
    yintVars = NULL;
}



int IQD_Callbacks::setSolverPreRelax( const unsigned int thnumber, OPT_Solver *solver )
{
    //const bool setObj = true;
    int ret;
    //const double *oly = refProb2->ly, *ouy = refProb2->uy;
    
    
    //refProb2->rprob.writeMIQCPModelInAMPLFile("novorefprob2.mod");
    
    
    ret = OPT_setMINLPProblem( refProb2->rprob, solver, true, true, true, bb->in_integrality_strategy == IQD_IS_ON_EACH_SUBPROBLEM);
    
    if( ret != 0 )
    {
        if(bb->in_print_level > 0)
            IQD_PRINTERRORMSG("Solver error");
        
        return IQD_SOLVER_ERROR;
    }
    
    
    if(solverParams)
    {
        ret = solver->setParameters( *solverParams );
        
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                std::cerr << IQD_PREPRINT << "Warning: Error on set list of parameters!" << IQD_GETFILELINE << endl;
        }
    }
    
    
    ret = solver->setNumberOfThreads(1);
    #if IQD_DEBUG_MODE
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                std::cerr << IQD_PREPRINT << "Warning: Error on set list of parameters!" << IQD_GETFILELINE << endl;
        }
    #endif
    
    
    //( (OPT_LPSolver *) solver)->generateModelFile("iquad_prerelax.lp");
    
    
    
    ret = solver->getNumberOfConstraints( indsecconstrs );
    
    #if IQD_DEBUG_MODE
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                std::cerr << IQD_PREPRINT << "Warning: Error on set list of parameters!" << IQD_GETFILELINE << endl;
        }
    #endif
    
    
    ret = solver->addConstraints( refProb2->dimv ); //add space for secants...
    
    if( ret != 0 )
    {
        if(bb->in_print_level > 0)
            std::cerr << IQD_PREPRINT << "Solver error" << IQD_GETFILELINE << endl;
        
        return IQD_SOLVER_ERROR;
    }
    
    if(bb->in_set_special_solver_parameters)
        setSpecialSolverParameters(solver);
    
    
    //std::cout << "refProb2->dimv: " << refProb2->dimv << "\n";
    //( (OPT_LPSolver *) solver)->generateModelFile("preiquad.lp");
    //MRQ_getchar();
    
    return 0;
}



int IQD_Callbacks:: setSpecialSolverParameters( optsolvers::OPT_Solver *solver )
{
    int code = 0;
    const int scode = solver->getSolverCode();
    
    
    if( scode == OPT_MOSEK )
    {
        code = solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_REL_GAP", 1.0e-4 );
        
        code += solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_MU_RED", 1.0e-6 );
        
        code += solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_PFEAS", 1.0e-6);
        
        code += solver->setDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_DFEAS", 1.0e-6 );
        
        code += solver->setDoubleParameter( "MSK_DPAR_INTPNT_TOL_INFEAS", 1.0e-6 );
        
        code += solver->setIntegerParameter( "MSK_IPAR_CHECK_CONVEXITY", 0 );
        
        if( refProb2->rprob.hasNLConstraints() || refProb2->rprob.hasObjNLTerm() )
        {
            #if OPT_HAVE_MOSEK
                code += solver->setIntegerParameter( "MSK_IPAR_INTPNT_STARTING_POINT", MSK_STARTING_POINT_SATISFY_BOUNDS );
            #endif
        }
        
    }
    else if( scode == OPT_IPOPT )
    {
        code = solver->setStringParameter("expect_infeasible_problem", "yes");
        //code += solver->setDoubleParameter("tol", 1.0e-6);
        //code += solver->setDoubleParameter("acceptable_tol", 1.0e-4);
        //code += solver->setDoubleParameter("mu_target", 1.0e-8);
    }
    
    
    return code == 0 ? 0 : IQD_BAD_DEFINITIONS;
}






//note that the first n indices in lx and ux has the bounds for the original variables. The next indices are bounds for y, not to w...
int IQD_Callbacks::beforeAll(const unsigned int numberOfThreads, double *lx, double *ux)
{
    const int n = refProb2->oprob->n;
    const int dimv = refProb2->dimv;
    int ret;
    
    nthreads = numberOfThreads;
    
    lastZuRecalcYBounds = IQD_INFINITY;
    yBoundsUpdated = false;
    
    origVarsUptd = refProb2->oprob->getNumberOfIntegerVars() > 0 ;//&& !IQD_isDiagonalSplitQWay(bb->in_split_way); //do not put if before allocateMemory, because
    
    
    mtotal =  refProb2->rprob.getNumberOfConstraints() + dimv; //set mtotal before allocate memory...
    
    
    if( refProb2->rprob.nlEval )
    {
        ret = ( (IQD_RefProb2NonLinearEval *)  refProb2->rprob.nlEval)->allocate(numberOfThreads);
        
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
    }
    
    
    ret = allocateMemory( numberOfThreads, mtotal ) + allocateSolverObjects( numberOfThreads );
    
    if( ret != 0 )
    {
        if(bb->in_print_level > 0)
            IQD_PRINTMEMERROR;
        return IQD_MEMORY_ERROR;
    }
    
    
    for(unsigned int i = 0; i < numberOfThreads; i++)
    {
        ret += setSolverPreRelax( i, solvers[i] );
    }
    
    if( ret != 0 )
    {
        if(bb->in_print_level > 0)
            IQD_PRINTERRORNUMBER(ret);
        return IQD_SOLVER_ERROR;
    }
    
    
    if( bb->in_recalculate_aux_bounds_on_zu_updt )
    {
        ret = allocateAuxBoundsUpdtObjects( numberOfThreads );
        
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
        
        for(unsigned int i = 0; i < numberOfThreads; i++)
            ret += auxBoundsUpdate[i].allocateAuxArrays(n + dimv);
        
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
        
        for(unsigned int i = 0; i < numberOfThreads; i++)
            ret += auxBoundsUpdate[i].buildPreProblem( refProb2->rprob, dimv, bb->in_qp_solver, bb->in_qcp_solver, bb->in_nlp_solver );
        
        
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
    }
    
    
    
    
    
    hoursToSaveFile = bb->in_hours_to_save_file;
    
    #if IQD_DEBUG_MODE
        
        useWarmStart = bb->useWarmStart;
        
        if( useWarmStart )
            assert( MIP_isNonlinearProblemType( refProb2->rprob.getProblemType() ) );
    #endif
    
    
    
    if( bb->in_use_heuristics_on_integer_probs && refProb2->oprob->getNumberOfIntegerVars() > 0 )
    {
        bool feasSol;
        const double MAXCPUTIME = 120;
        
        double obj;
        double *sol = values[0];
        
        
        ret = IQD_runMuriquiHeuristic( refProb2->rprob, bb->in_qp_solver, bb->in_nlp_solver, NULL, solverParams, INFINITY, MAXCPUTIME, obj, sol );
        
        if( ret == 0 )
        {
            const double *sol = refProb2->rprob.x;
            
            
            ret = refProb2->oprob->isFeasibleToConstraints( 0, sol, true, constrEval,  bb->in_abs_feas_tol, bb->in_rel_feas_tol, feasSol, solvers[0]->constr );
            
            
            if( ret == 0 && feasSol )
            {
                double obj;
                
                ret = refProb2->oprob->objEval(0, true, sol, obj);
                
                if( ret == 0 )
                {
                    tryUpdateBestSolution(0, refProb2->rprob.x, obj, 1);
                    
                    
                    if( bb->in_print_level > 2 )
                        std::cout << IQD_PREPRINT << "Feasible soluiton found by Muriqui. Obj value: " << obj << std::endl;
                }
            }
        }
        
    }
    
    
    if( bb->in_build_split_without_vars_also && refProb2->lambdaObjNoInt)
    {
        //we store coefficients related to integer variables in objective function
        const int *xtype = refProb2->oprob->xtype;
        const IQD_SparseMatrix &Q = refProb2->oprob->Q;
        double *values = this->values[0];
        
        
        Qir = new (std::nothrow) IQD_SparseMatrix;
        
        if( !Qir )
        {
            if(bb->in_print_level > 0)
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
        
        
        ret = Qir->initialize( nI, n, false );
        if( ret != 0 )
        {
            if(bb->in_print_level > 0)
                IQD_PRINTMEMERROR;
            return IQD_MEMORY_ERROR;
        }
        
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            //IQD_setAllarray( n, values, 0.0 );
            
            Q.getFullRow(ind, values, false, false, 1.0);
            
            
            for(int j = ind + 1; j < n; j++)
            {
                if( MIP_isIntegerType( xtype[j] ) )
                    continue; //we are only interested in continuous vars. Otherwise, we will count two times terms between integer vars...
                
                //IQD_SparseRow &row = Q[j];
                const int* Qrcols = Q[j];
                const double* Qrvals = Q(j);
                unsigned int nel = Q.getNumberOfElementsAtRow(j); //row.getNumberOfElements();
                
                for(unsigned int k = 0; k < nel; k++)
                {
                    if( Qrcols[k] == ind ) //if( (int) row[k].getColumn() == ind )
                    {
                        values[j] = Qrvals[k]; //row[k].getValue();
                        break;
                    }
                }
                
            }
            
            
            
            values[ind] *= 0.5; //diagonal position should be multiply by 0.5
            
            ret = Qir->setRowStructureAndValues(i, values);
            
            if( ret != 0 )
            {
                if(bb->in_print_level > 0)
                    IQD_PRINTMEMERROR;
                return IQD_MEMORY_ERROR;
            }
            
        }
        
        
        //std::cout << "Qir: " << "Qir.symmetric: " << Qir->symmetric  << std::endl;
        //Qir->printSparseMatrix();
        //IQD_getchar();
    }
    
    
    //for(int i = 0; i < refProb2->dimv; i++)
        //cout << "oly[" << i << "]: " << refProb2->ly[i] << " ouy[" << i << "]: " << refProb2->uy[i] << endl;
    
    
    return 0;
}



int IQD_Callbacks:: updateAuxVariablesBounds(const unsigned int thnumber, OPT_LPSolver *solver, const double *ly, const double *uy )
{
    const int dimv = refProb2->dimv;
    const int n = refProb2->oprob->n;
    
    int r, code = 0;
    double lb, ub;
    
    for(int i = 0; i < dimv; i++)
    {
        if( uy[i] < 0.0 )
            ub = - uy[i] * uy[i];
        else if ( ly[i] > 0.0 )
            ub = - ly[i] * ly[i];
        else
            ub = 0.0;
        
        lb = IQD_min( -ly[i]*ly[i], -uy[i]*uy[i] );
        
        r = solver->setVariableBounds( i+n, lb, ub );
        
        if( r != 0 )
            code = IQD_SOLVER_ERROR;
    }
    
    
    
    return code;
}



int IQD_Callbacks:: beforeSolvingRelaxation( const unsigned int threadNumber, branchAndBound::BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, bool &pruneNode)
{
    const IQD_SparseMatrix &v = refProb2->v;
    pruneNode = false;
    
    //for(int i = 0; i < refProb2->rprob.n; i++)
        //std::cout << "nlx[" << i << "]: " << nlx[i] << " nux[" << i << "]: " << nux[i] << std::endl;
    
    //node.print();
    //IQD_getchar();
    
    
    if( yBoundsUpdated )
    {
        const int n = refProb2->oprob->n;
        const int dimv = refProb2->dimv;
        
        double *nly = &nlx[n];
        double *nuy = &nux[n];
        
        for( int i = 0; i < dimv; i++ )
        {
            if( ly[i] > nly[i] )
            {//we already have a better lower bound;
                //cout << "vou atualizar nly[" << i << "]! nly[" << i <<"]: " << nly[i] << " ly[" << i << "]: " << ly[i] << endl;
                
                nly[i] = ly[i];
                
                if( nly[i] > nuy[i] )
                {
                    pruneNode = true;
                    return 0;
                }
                
            }
        }
        
        
        for( int i = 0; i < dimv; i++ )
        {
            if( uy[i] < nuy[i] )
            {//we already have a better upper bound;
                //cout << "vou atualizar nuy[" << i << "]! nuy[" << i <<"]: " << nuy[i] << " uy[" << i << "]: " << uy[i] << endl;
                
                nuy[i] = uy[i];
                
                if( nly[i] > nuy[i] )
                {
                    pruneNode = true;
                    return 0;
                }
            }
        }
    }
    
    
    
    if( origVarsUptd && !IQD_isDiagonalSplitQWay(bb->in_split_way) )
    {
        const int n = refProb2->oprob->n;
        const int dimv = refProb2->dimv;
        
        double *nly = &nlx[n];
        double *nuy = &nux[n];
        
        double sumlb, sumub;
        
        //we check if new bounds for original variables can change the bounds for y
        
        for(int i = 0; i < dimv; i++)
        {
            //IQD_SparseRow &rowv = refProb2->v[i];
            const int* vrcols = v[i];
            const double* vrvals = v(i);
            const int nel = v.getNumberOfElementsAtRow(i); //rowv.getNumberOfElements();
            
            
            sumlb = sumub = 0.0;
            
            for(int j = 0; j < nel; j++)
            {
                const int ind  = vrcols[j]; //rowv[j].getColumn();
                const double v = vrvals[j]; //rowv[j].getValue();
                
                if( v < 0.0 )
                {
                    sumlb += v * nux[ind];
                    sumub += v * nlx[ind];
                }
                else //v >= 0.0
                {
                    sumlb += v * nlx[ind];
                    sumub += v * nux[ind];
                }
            }
            
            if( sumlb > nly[i]  )
            {
                //cout << "Atualizei nly["<< i << "]! Antigo: " << nly[i] << " Novo: " << sumlb << std::endl;
                nly[i] = sumlb;
            }
            
            if( sumub < nuy[i] )
            {
                //cout << "Atualizei nuy["<< i << "]! Antigo: " << nuy[i] << " Novo: " << sumub << std::endl;
                nuy[i] = sumub;
            }
        }
        
    }
    
    
    return 0;
}




int IQD_Callbacks::solveSubProblem( const unsigned int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES &retCode, double &objValue, double &dualObjValue, double *sol, double *dualSol, bool &generalFeasibleSol, bool &pruneNode, double &nodeLowerBound, BBL_BRANCH_STRATEGY &branchStrategy)
{
    const IQD_IQuadProb *prob = refProb2->oprob;
    
    const unsigned int n = prob->n;
    const int dimv = refProb2->dimv;
    const int indvconstrs = refProb2->indvconstrs;
    
    const unsigned int nprimal = refProb2->rprob.getNumberOfVars();
    const unsigned int ndual = refProb2->getNDual();
    
    const double maxcputime = bb->in_max_cpu_time;
    const double *ly = &nlx[n], *uy = &nux[n];
    
    
    //bool *constrEval = constrEvals[thnumber];
    double *vals = values[thnumber] ;
    int *cols = columns[thnumber];
    
    const IQD_SparseMatrix &v = refProb2->v;
    
    
    OPT_QPSolver* solver = (OPT_QPSolver *) solvers[thnumber];
    
    int ret = 0;
    
    
    
    
    //generalFeasibleSol and pruneNode are already initialized as false by BBL_BranchAndBound... retCode is already initialized as BBL_UNDEFINED
    
    #if IQD_MURIQUI_GAP_MIN_TREATMENT
        if( bb->mrq_gap_min )
        {
            //we are addressing the gap minimization problem from igma
            const double zerotol = bb->in_integer_tol;
            const double onetol = 1 - bb->in_integer_tol;
            //int r, nneg = 0 ;
            
            for(int i = 0; i < nI; i++)
            {
                //double objCoef;
                const int yind = yintVars[i]; //we have a diagonal split
                
                
                if( uy[yind] < onetol )
                {
                    
                    if(ly[yind] > zerotol)
                    {
                        //we can prune this node because we have a binary variable ha cannot be zero nor 1
                        pruneNode = true;
                        IQD_PRINTERRORMSG("PODEI");
                        IQD_getchar();
                        return 0;
                    }
                    
                    //vals[i] = 1.0;
                }
                /*else if( ly[yind] > zerotol )
                {
                    vals[i] = -1.0;
                    nneg++;
                }
                else
                {
                    vals[i] = 0.0;
                }*/
                
            }
            
            
            /*r = solver->setLinearObjCoefs(nI, intVars, vals);
            r += solver->setObjConstant(nneg);
            #if IQD_DEBUG_MODE
                if( r != 0 )
                {
                    IQD_PRINTERRORNUMBER(r);
                    return IQD_SOLVER_ERROR;
                }
            #endif */
            
            
            
        }
    #endif
    
    
    
    
    
    
    if( refProb2->lambdaObjNoInt )
    {
        const bool intfix = IQD_isIntVarsFixed(nI, intVars, nlx, nux);
        const double *lambdaObjNoInt = refProb2->lambdaObjNoInt;
        char qstype = QsplitType[ thnumber ];
        int r;
        
        
        if( intfix )
        {
            int nzLambda = 0;
            const double zero_tol = refProb2->zero_tol;
            const IQD_SparseMatrix *PnoInt = refProb2->PnoInt;
            
            
            
            if( prob->hasLinCoefObj() )
                IQD_copyArray(n, prob->c, vals);
            else
                IQD_setAllarray(n, vals, 0.0);
            
            
            for( int i = 0; i < nI; i++ )
            {
                const int ind = intVars[i];
                const double f = nlx[ind];
                
                if( f != 0.0 )
                {
                    /*IQD_SparseRow &row = (*Qir)[i];
                    
                    unsigned int nel = row.getNumberOfElements();
                    
                    for(unsigned int j = 0; j < nel; j++)
                        vals[ row[j].getColumn() ] += f*row[j].getValue(); */
                    
                    Qir->accumulateRowInArray(i, vals, f);
                }
            }
            
            if( prob->objFactor != 1.0 )
                IQD_multiplyAllArray( n, vals, prob->objFactor );
            
            
            #pragma ivdep
            #pragma GCC ivdep
            for(int i = 0; i < (int)n; i++)
                cols[i] = i;
            
            r = solver->setObjLinearPart(n, vals);
            #if IQD_DEBUG_MODE
                if(r != 0)
                {
                    if( bb->in_print_level > 0 )
                        IQD_PRINTERRORNUMBER(r);
                    return IQD_SOLVER_ERROR;
                }
            #endif
            
            
            
            if(qstype != 1)
            {
                
                for( int j = 0; j < dimv; j++ )
                {
                    if( lambdaObjNoInt[j] < -zero_tol )
                    {
                        cols[nzLambda] = n+j;
                        vals[nzLambda] = -0.5 * lambdaObjNoInt[j];
                        nzLambda++;
                    }
                }
                
                
                r = solver->setObjLinearCoefs( nzLambda, cols, vals );
                #if IQD_DEBUG_MODE
                    if(r != 0)
                    {
                        if( bb->in_print_level > 0 )
                            IQD_PRINTERRORNUMBER(r);
                        return IQD_SOLVER_ERROR;
                    }
                #endif
                
                r = solver->setObjQuadMatrix(*PnoInt);
                if(r != 0)
                {
                    if( bb->in_print_level > 0 )
                        IQD_PRINTERRORNUMBER(r);
                    return IQD_SOLVER_ERROR;
                }
                
                
                qstype = 1;
                
                
                std::cout << "Alternei para P sem vars inteiras" << std::endl;
                IQD_getchar();
            }
            
        }
        else
        {
            if( qstype != 0 )
            {
                const double zero_tol = refProb2->zero_tol;
                const double *lambdaObj = refProb2->lambdaObj;
                const MIP_SparseMatrix &P = refProb2->rprob.Q;
                int nzLambda = 0;
                double *pv;
                
                
                if( prob->hasLinCoefObj() )
                {
                    
                    if( prob->objFactor == 1.0)
                    {
                        pv = prob->c;
                    }
                    else
                    {
                        IQD_copyArray( n, prob->c, vals );
                        IQD_multiplyAllArray(n, vals, prob->getObjFactor() );
                        
                        pv = vals;
                    }
                }
                else
                {
                    IQD_setAllarray(n, vals, 0.0);
                    pv = vals;
                }
                
                
                r = solver->setObjLinearPart(n, pv);
                #if IQD_DEBUG_MODE
                    if( r != 0 )
                    {
                        if( bb->in_print_level > 0 )
                            IQD_PRINTERRORNUMBER(r);
                        return IQD_SOLVER_ERROR;
                    }
                #endif
                
                
                for(int j = 0; j < dimv; j++)
                {
                    if( lambdaObj[j] < -zero_tol )
                    {
                        cols[nzLambda] = n+j;
                        vals[nzLambda] = -0.5*lambdaObj[j];
                        
                        nzLambda++;
                    }
                }
                
                
                r = solver->setObjLinearCoefs( nzLambda, cols, vals );
                #if IQD_DEBUG_MODE
                    if( r != 0 )
                    {
                        if( bb->in_print_level > 0 )
                            IQD_PRINTERRORNUMBER(r);
                        return IQD_SOLVER_ERROR;
                    }
                #endif
                
                
                r = solver->setObjQuadMatrix(P);
                if( r != 0 )
                {
                    if( bb->in_print_level > 0 )
                        IQD_PRINTERRORNUMBER(r);
                    return IQD_SOLVER_ERROR;
                }
                
                qstype = 0;
                
                std::cout << "Alternei para P com vars inteiras" << std::endl;
                IQD_getchar();
            }
        }
        
    }
    
    
    if( origVarsUptd )
    { //so, we are branching on original variables space also:
        
        ret += solver->setnVariablesBounds( n, nlx, nux );
    }
    
    
    
    //updating auxiliary constraints...
    for(int i = 0, k = 0; i < dimv; i++)
    {
        //IQD_SparseRow &vrow = v[i];
        
        const unsigned int nel = v.getNumberOfElementsAtRow(i); //vrow.getNumberOfElements();
        
        if( nel == 1 )
        {
            const double d = v(i)[0]; //vrow[0].getValue();
            const int ind = v[i][0]; //vrow[0].getColumn();
            
            double lb = ly[i];
            double ub = uy[i];
            
            if( d != 1.0 )
            {
                lb = lb/d;
                ub = ub/d;
                
                if( d < 0 )
                    IQD_swap(lb, ub);
            }
            
            
            #if IQD_DEBUG_MODE
                ret +=
            #endif
            
            solver->setVariableBounds( ind, lb, ub );
            
        }
        else if( nel >= 1 )
        {
            #if IQD_DEBUG_MODE
                ret +=
            #endif
            
            solver->setConstraintBounds( indvconstrs + k, ly[i], uy[i] );
            
            k++;
        }
    }
    
    #if IQD_DEBUG_MODE
        if( ret != 0 )
        {
            if( bb->in_print_level > 0 )
                std::cerr << IQD_PREPRINT << "Error at setting subproblem on node" << IQD_GETFILELINE << endl ;
            
            return IQD_SOLVER_ERROR;
        }
    #endif
    
    ret = updateAuxVariablesBounds(thnumber, solver, ly, uy);
    
    if( ret != 0 )
    {
        if( bb->in_print_level > 0 )
            std::cerr << IQD_PREPRINT << "Error at setting subproblem on node" << IQD_GETFILELINE << endl ;
        
        return IQD_SOLVER_ERROR;
    }
    
    
    ret = IQD_setSecantConstraints(n, dimv, v, indsecconstrs, solver, ly, uy, bb->in_set_secants_as_equality_constraints, cols, vals );
    
    if( ret != 0 )
    {
        if( bb->in_print_level > 0 )
            std::cerr << IQD_PREPRINT << "Error at setting subproblem on node" << IQD_GETFILELINE << endl ;
        
        return IQD_SOLVER_ERROR;
    }
    
    
    if( useWarmStart )
    {
        if( node.buildInitialSolution(nprimal, vals) == 0  &&  node.getParentDualSolution(nprimal, ndual, &vals[nprimal] ) == 0)
        {
            
            ( (OPT_NLPSolver*) solver)->setInitialSolution( vals, &vals[nprimal], &vals[mtotal] );
        }
        
    }
    
    if( ub < bb->in_infinity_value )
        solver->setObjCutUpperBound( ub );//solver->setUpperObjCut( ub );
    
    
    #if IQD_DEBUG_MODE
        if(iter == 1u)
        {
            if( bb->in_relax_prob_writing_mode != IQD_RPWM_NO_WRITING)
            {
                const char name[] = "_refprobiquad.lp";
                
                solver->generateModelFile(name);
            }
        }
    #endif
    
    
    //solver->generateModelFile( "iquadnoraiz.lp" );
    //IQD_getchar();
    
    if( maxcputime < INFINITY )
    {
        const double tremin = maxcputime - ( double(clock() - bb->clockStart) )/CLOCKS_PER_SEC;
        
        solver->setMaxTime(tremin);
    }
    
    
    ret = solver->solve(false);
    
    
    
    //cout << "solver->return code: " << solver->retCode << " obj: " << solver->objValue << " origSolverRetCode: " << solver->origSolverRetCode << endl;
    //for(int i = 0; i < n + dimv; i++)
    //for(int i = 0; i < n ; i++)
        //std::cout << "x["<<i<<"]: " << solver->sol[i] << "   \t";// std::endl;
    //std::cout << std::endl;
    
    //cout << "gerei iquad.lp" << endl;
    //solver->generateModelFile("iquad.lp");
    
    
    if( ret == OPT_OPTIMAL_SOLUTION )
    {
        dualObjValue = solver->objValue;
        
        ret = prob->objEval(thnumber, true, solver->sol, objValue, 1.0); 
        
        if( ret != 0 )
            objValue = NAN;
        
        IQD_copyArray( nprimal, solver->sol, sol );
        if( useWarmStart )
        {
            solver->getDualSolution( dualSol,  &dualSol[mtotal], false );
            
            //IQD_copyArray( mtotal, solver->dualSolC, dualSol );
            //IQD_copyArray( 2*nprimal, solver->dualSolV, &dualSol[mtotal] );
        }
        
        
        
        if( nI > 0 && bb->in_integrality_strategy != IQD_IS_ON_EACH_SUBPROBLEM)
        {
            generalFeasibleSol = IQD_isIntegerSol( nI, intVars, solver->sol, bb->in_integer_tol );
        }
        else
            generalFeasibleSol = true; 
        
        
        if( generalFeasibleSol )
        {
            //we will store in solve the constraints value in solver object. Maybe it is not ideal, but in this way we take advantage the space already allocated and is easier to pass it to method to choose variable to do branching...
            
            if( refProb2->nonConvexConstrs )
            {
                ret = prob->isFeasibleToConstraints( thnumber, sol, true, constrEval, bb->in_abs_feas_tol, bb->in_rel_feas_tol, generalFeasibleSol, solver->constr );
                
                if( ret != 0 )
                    generalFeasibleSol = false;
            }
        }
        else
        {
            #if IQD_DEBUG_MODE
                IQD_setAllarray( prob->m, solver->constr, (double) NAN );
            #endif
        }
        
        
        if( iter == 1u )
        {
            if( generalFeasibleSol )//if( generalFeasibleSol && iter == 1u )
                bb->out_root_node_upper_bound = objValue;
        }
        
        
        //we just store in solver object to pass the information to chooseIndexToBranch method.
        solver->feasSol = generalFeasibleSol;
        
        
        retCode = BBL_OPTIMAL_SOLUTION;
        
        
        #if IQD_DEBUG_MODE
// 			//we use bb->in_zero_tol_to_split as relative tolerance
            if(ret == 0 && objValue + IQD_abs(solver->objValue * bb->in_zero_tol_to_split * 100	) + bb->in_zero_tol_to_split*1000 < solver->objValue)
            {
                std::cerr << IQD_PREPRINT << "Relaxation obj function is greater than real obj function at IQD_BranchAndBound::bbloop relaxation! Relaxation:" << setprecision(10) << solver->objValue << " real value: " << objValue << IQD_GETFILELINE << std::endl;
                
                double robj;
                
                refProb2->rprob.objEval(thnumber, true, solver->sol, robj );
                
                std::cout << "robj: " << robj << std::endl;
                std::cout << "dimv: " << dimv << std::endl;
                
                for(int j = 0; j < dimv; j++)
                {
                    double lb = NAN, ub = NAN;
                    
                    solver->getVariableBounds(n+j, lb, ub);
                    
                    
                    std::cerr << "x[" << j << "]: " << solver->sol[j] << " w[" << j << "]: " << solver->sol[n+j] << " wlb[" << j << "]: " << lb << " wub[" << j << "]: " << ub << std::endl;
                }
                
                if( objValue - solver->objValue >= 1.0 )
                    assert(false);
            }
        #endif
        
        
    }
    else if( ret == OPT_INFEASIBLE_PROBLEM )
    {
        retCode = BBL_INFEASIBLE_PROBLEM;
        //std::cout <<  IQD_PREPRINT <<  "Infeasible problem\n";
    }
    
    //we let other return codes being treated like BBL_UNDEFINED
    
    
    
    
    
    
    //cout << "retCode: " << retCode << " objValue: " << objValue << " dualObj: " << dualObjValue << " solver.objValue: " << solver->objValue << " generalFeasibleSol: " << generalFeasibleSol  << endl;
    
    
    //IQD_getchar();
    
    return 0;
}








int IQD_Callbacks::chooseIndexToBranch(const int thnumber, BBL_Node &node, const long unsigned int iter, const double lb, const double ub, double *nlx, double *nux, BBL_RETURN_CODES retCode, double objValue, double dualObjValue, double *sol, double *dualSol, unsigned int &sizeIndices, unsigned int *indices, double *breakValues1, double *breakValues2, BBL_Node* &nodes)
{
    const bool feasSol = solvers[thnumber]->feasSol;
    const double *constrs = solvers[thnumber]->constr;
    const double ZEROTOLTOMAXDIF = 1.0e-5;
    const char qstype = QsplitType[ thnumber ];
    const double *lambdaObj = qstype == 0 ?  refProb2->lambdaObj : refProb2->lambdaObjNoInt ;
    double **lambdaConstr = refProb2->lambdaConstr;
    const IQD_SparseMatrix &v = refProb2->v;
    
    bool branchOnObjLambda, intBranching;
    unsigned int ind;
    double maxIGap = 0.0;
    
    int *count = columns[thnumber];
    double *vals = values[thnumber] ;
    double *gaps = vals, *y = &vals[refProb2->dimv];
    
    const double factorToPrioObj = bb->in_factor_to_prior_obj_f_on_branch;
    double denFactor = 1.0e-6;
    
    
    
    
    if( IQD_abs(objValue) > denFactor  )
    {
        denFactor = 0.0; //we do not need sum nothing to denominator...
    }
    
    
    //cout << "Escolhendo indice para branching" << endl;
    
    
    refProb2->calculatey(sol, y);
    
    
    /*{
        const int n = refProb2->oprob->n;
        const int m = refProb2->oprob->m;
        const double *ly = &nlx[n], *uy = &nux[n];
        
        
        cout << setprecision(10);
        
        cout << "feasSol: " << feasSol << " feasTol: " << bb->in_feas_tol << endl;
        
        for(int i = 0; i < m; i++)
            cout << "constrs[" << i << "]: " << constrs[i] << endl;
        
        
        const double *w = &sol[ n ];
        //cout << "w: " << endl;
        for( int i = 0; i < refProb2->dimv; i++ )
            cout << "ly[" << i << "]: " << ly[i] << " uy[" << i << "]: " << uy[i]  <<  " \ty[" << i << "]: " << y[i] << " \tw[" << i << "]: " << w[i] << endl;
    } */
    
    //std::cout << "feasSol: " << feasSol << std::endl;
    
    //if( feasSol )
        //cout << "Solucao viavel!" << endl;
    
    intBranching = false;
    if( !feasSol && bb->in_integrality_strategy != IQD_IS_ON_EACH_SUBPROBLEM )
    {
        //cheking if there is a integer variable in a non integer value
        const double itol = bb->in_integer_tol;
        
        for(int i = 0; i < nI; i++)
        {
            //const double v = sol[ intVars[i] ];
            
            gaps[i] = IQD_intGap( sol[ intVars[i] ] ); // IQD_abs(round(v) - v);
            
            if( gaps[i] > itol )
            {
                intBranching = true; //if we do a integer branching, we need the gaps calculated. So, do not BREAK here
                
                if( gaps[i] > maxIGap )
                    maxIGap = gaps[i];
            }
            
            //std::cout << "gaps[" << intVars[i] << "]: " << gaps[i] << std::endl;
        }
        
    }
    
    
    if( intBranching )
    {
        //we will choose a integer variable to perform the branch. To do it, we perform a balance between the integer value and a possible objective diference between a strtict auxiliary variable related 
        const int n = refProb2->oprob->n;
        const double itol = bb->in_integer_tol;
        const double ialpha = bb->in_alpha_to_integer_branch_choosing;
        
        int ind2;
        double maxObjDif = 0.0;
        double v, maxv = 0.0, objDif;
        
        
        //printf("oi 1\n");
        
        
        for( int i = 0; i < nI; i++ )
        {
            const int yind = yintVars[i];
            
            if( yind >= 0 && gaps[i] > itol)
            {
                objDif = IQD_calDif( lambdaObj[yind], y[yind], sol[n+yind] );
                //lambda[yind]*( y[yind]*y[yind] + sol[n+yind] ) ;
                
                if( objDif > maxObjDif )
                    maxObjDif = objDif;
                
                //std::cout << "index: " << intVars[i] << " gap: " << gaps[i] << " objDif: " << objDif << std::endl;
            }
        }
        
        //std::cout << "maxobjdif: " << maxObjDif << std::endl;
        
        
        #if IQD_DEBUG_MODE
            assert( maxIGap > 0.0 );
        #endif
        
        //calculating now the balance. Note the maximum gap is 0.5. So, we multiply by 2, in this way, gap weight is 1.0, the same weight of objective dif. Note, we normalized the obj dif dividing by maxObjDif. In this way, the weight of obj dif is 1 also.
        const double igapfactor = ialpha/maxIGap;
        ind = UINT_MAX;
        for( int i = 0; i < nI; i++ )
        {
            //std::cout << "i: " << i << " intVars["<<i<<"]: " << intVars[i] << " gaps["<<i<<"]: " << gaps[i] << " itol: " << itol <<  std::endl;
            
            if( gaps[i] > itol )
            {
                const int yind = yintVars[i];
                
                v = igapfactor *gaps[i];//ialpha* 2*gaps[i];
                
                if( yind >= 0 && maxObjDif > ZEROTOLTOMAXDIF )
                {
                    v += (1.0 - ialpha ) * IQD_calDif( lambdaObj[yind], y[yind], sol[n+yind] )/maxObjDif; //we have to consider the obj difference
                }
                
                //std::cout << "v: " << v << std::endl;
                
                if( v > maxv )
                {
                    maxv = v;
                    ind = yintVars[i] >= 0 ? yintVars[i] + n: intVars[i] ; //ind has the index of variable to integer branching. If the variable is related to auxiliary variable y, we prefer branch on y
                    ind2 = intVars[i];
                }
            }
        }
        
        #if IQD_DEBUG_MODE
            assert( ind < UINT_MAX );
        #endif
        
        //std::cout << "ind: " << ind << " ind2: " << ind2 << std::endl;
        
        indices[0] = ind;
        breakValues1[0] = floor( sol[ind2] );
        breakValues2[0] = breakValues1[0] + 1.0;
        
        sizeIndices = 1;
    }
    else
    {
        branchOnObjLambda = feasSol || ( bb->in_priorize_obj_f_on_branch_choosing  &&  (objValue - dualObjValue)/(IQD_abs(objValue) + denFactor) > factorToPrioObj )  ;
        
        
        
        if( branchOnObjLambda )
        {
            //std::cout << "branching OnObjLambda " << std::endl;
            myChooseIndexToBranch( lambdaObj, 1, sol, y, bb->in_use_real_dif_on_branch_choosing, gaps,  sizeIndices, &ind );
        }
        else
        {
            const int mo = refProb2->oprob->m;
            const double *lc = refProb2->oprob->lc;
            const double *uc = refProb2->oprob->uc;
            
            const double absFeasTol = bb->in_abs_feas_tol;
            const double relFeasTol = bb->in_rel_feas_tol;
            
            
            if( bb->in_branch_strategy == IQD_BS_MAX_CONSTRS_VIOLATED )
            {
                const int dimv = refProb2->dimv;
                const int *origConstr = refProb2->origConstr;
                double **lambdaConstr = refProb2->lambdaConstr;
                const double zero_tol = refProb2->zero_tol;
                
                int maxinds;
                
                double *lambda;
                
                
                
                for(int j = 0; j < dimv; j++)
                    count[j] = 0;
                
                for(int i = 0; i < mo; i++)
                {
                    if( !constrEval[i] )
                        continue;
                    
                    const double ltol = absFeasTol + IQD_abs( lc[i] * relFeasTol );
                    
                    
                    bool viol = false;
                    
                    if( constrs[i] < lc[i] - ltol )
                    {
                        viol = true;
                        if( uc[i] < MIP_INFINITY )
                            lambda = lambdaConstr[ origConstr[i] ];
                        else
                            lambda = lambdaConstr[ i ];
                    }
                    else
                    {
                        const double utol = absFeasTol + IQD_abs( uc[i] * relFeasTol );
                        
                        if( constrs[i] > uc[i] + utol )
                        {
                            viol = true;
                            lambda = lambdaConstr[i];
                        }
                    }
                    
                    
                    if(viol)
                    {
                        
                        #if IQD_DEBUG_MODE
                            assert(lambda != NULL);
                        #endif
                        
                        for( int j = 0; j < dimv; j++ )
                        {
                            if( lambda[j] < -zero_tol )
                                count[j]++;
                        }
                    }
                }
                
                maxinds = 0;
                for(int j = 0; j < dimv; j++)
                {
                    if( count[j] > maxinds )
                    {
                        maxinds = count[j];
                        ind = j;
                    }
                }
                
                sizeIndices = 1;
            }
            else
            {
                int indMostViol;
                double viol, mostViol = 0.0;
                const int *origConstr = refProb2->origConstr;
                
                //bb->in_feas_tol;
                
                #if IQD_DEBUG_MODE
                    indMostViol = INT_MAX;
                #endif
                
                
                
                for(int i = 0; i < mo; i++ )
                {
                    if( !constrEval[i] )
                        continue;
                    
                    
                    const bool lambdac = lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY ? lambdaConstr[ origConstr[i] ] != NULL: lambdaConstr[i] != NULL ; //if constraint is double bounded, part referent to lc is in the index origConstr[i]. if lambdaConstr[k] is not NULL, is because constraint k is nonconvex...
                    
                    
                    
                    
                    const double ltol = absFeasTol + IQD_abs( lc[i] * relFeasTol );
                    
                    
                    //std::cout << "constrs["<<i<<"]: " << constrs[i] << " lc["<<i<<"]: " << lc[i] << " uc["<<i<<"]: " << uc[i] << " ltol: " << ltol << " lambdac: " << lambdac << std::endl ;
                    
                    
                    if( lambdac && constrs[i] < lc[i] - ltol )
                    {
                        //we divide by |lc| to try have a normalization
                        viol = (lc[i] - constrs[i])/( IQD_abs(lc[i]) + 1.0);
                        
                        if( viol > mostViol )
                        {
                            mostViol = viol;
                            
                            
                            //if uc is not infinity, that is a new constraint in the reformulated problem...
                            
                            indMostViol = uc[i] < MIP_INFINITY ? refProb2->origConstr[i] : i;
                        }
                    }
                    else
                    {
                        const double utol = absFeasTol + IQD_abs( uc[i] * relFeasTol );
                        
                        
                        if( lambdaConstr[i] != NULL && constrs[i] > uc[i] + utol )
                        {
                            //we divide by |uc| to try have a normalization
                            viol = (constrs[i] - uc[i])/( IQD_abs(uc[i]) + 1.0);
                            
                            if( viol > mostViol )
                            {
                                mostViol = viol;
                                indMostViol = i;
                            }
                        }
                    }
                }
                
                //cout << "Ramificando sobre restricao: " << indMostViol << endl;
                
                #if IQD_DEBUG_MODE
                    
                    if( indMostViol == INT_MAX )
                    {
                        OPT_QPSolver* solver = (OPT_QPSolver *) solvers[thnumber];
                        
                        solver->generateModelFile( "iquad.lp" );
                        
                        for(int i = 0; i < refProb2->oprob->n; i++)
                            printf("x[%d]: %f\n", i, sol[i]);
                    }
                    
                    assert( indMostViol != INT_MAX ); //no violation detected...
                    assert( refProb2->lambdaConstr[indMostViol] != NULL );
                #endif
                
                
                
                
                myChooseIndexToBranch( refProb2->lambdaConstr[indMostViol], 1, sol, y, bb->in_use_real_dif_on_branch_choosing,  gaps, sizeIndices, &ind);
            }
            
        }
        
        
        indices[0] = ind + refProb2->oprob->n; //ind has the index in y. So, we have sum the n first indices from the original variables...
        
        if( v.getNumberOfElementsAtRow(ind) == 1 && MIP_isIntegerTypeToBranch(refProb2->rprob.xtype[ v[ind][0] ]) ) //		if( refProb2->v[ ind ].nElements == 1 && MIP_isIntegerTypeToBranch( refProb2->rprob.xtype[  refProb2->v[ ind ][0].col ] )  ) //if( IQD_isDiagonalSplitQWay(bb->in_split_way) && MIP_isIntegerTypeToBranch( refProb2->rprob.xtype[ indices[0] ] )  )
        {
            //so, we can use a integer programming branching
            const int ind2 = v[ind][0]; //refProb2->v[ ind ][0].col;
            
            
            #if IQD_DEBUG_MODE
                //assert( refProb2->v[ indices[0] ].nElements == 1 );
                const double vvalue = v(ind)[0]; //refProb2->v[ind][0].v;
                assert(vvalue >= 0.99999 && vvalue <= 1.00001);
            #endif
            
            breakValues1[0] = floor( sol[ind2] ); //y_i and x_i are linked and has the same value...
            breakValues2[0] = breakValues1[0] + 1.0;
        }
        else
        {
            //we have to do a spatial branching...
            
            breakValues1[0] = breakValues2[0] = bb->calculateBranchingPoint( y[ind], nlx[ indices[0] ], nux[ indices[0] ], bb->in_alpha );
        }
        
    }
    
    //cout <<  "sizeIndices: " << sizeIndices << " indices[0]: " << indices[0] << " breakValues1[0]: " << breakValues1[0] << " breakValues2[0]: " << breakValues2[0] << endl ;
    
    
    return 0;
}





int IQD_Callbacks::endOfIteration( const int thnumber, const long unsigned int iter, const double cpuTime, const double wallTime, const double lb, const double ub, BBL_Node &node, const double *nlx, const double *nux)
{
    #if IQD_DEBUG_MODE
        
        static double lastlb = -BBL_INFINITY;
        
        if( iter == 1 )
            lastlb = -BBL_INFINITY;
        
        
        if( lb < lastlb )
        {
            OPT_LPSolver* solver = (OPT_LPSolver *) solvers[thnumber];
            
            const int n = refProb2->oprob->n;
            const int dimv = refProb2->dimv;
            
            cout << "Bizarro! Limite inferior diminuiu. iter: " << iter << " antigo lb: " << lastlb << " atual lb: " << lb << " ub: " << ub << endl;
            
            cout << "solver.rcode: " << solver->retCode << " solver.objValue: " << solver->objValue << " solver.dualObjValue: " << solver->dualObjValue << endl;
            
            
            cout << "node:  lb: " << node.getLowerBound() << endl;
            node.print();
            
            for(int i = 0; i < dimv; i++)
                cout << "nlx[" << n+i << "]: " << nlx[n+i] << " nux[" << n+i << "]: " << nux[n+i] << " \t" ;
            
            
            //getchar();
            assert( false );
        }
        
        
        lastlb = lb;
        
    #endif
    
    
    
    #if IQD_SAVE_OUTPUT_FILE_BY_HOUR
        
        if( thnumber == 0 && IQD_outFileHour )
        {
            double cpu_time = ((double) (clock() - bb->clockStart)) /CLOCKS_PER_SEC;
            
            if( cpu_time >= hoursToSaveFile*60*60  ) //only thread 0 can save information...
            {
                IQD_outFileHour->writeDouble(hoursToSaveFile); //hours spent
                IQD_outFileHour->writeDouble( lb ); //lower bound
                IQD_outFileHour->writeDouble( ub ); //upper bound
                IQD_outFileHour->writeInt( (int) iter );
                IQD_outFileHour->writeDouble( IQD_getTime()- bb->timeStart ); //watch time
                
                IQD_outFileHour->flush();
                
                IQD_outFileHour->lastHoursToSaveFile = hoursToSaveFile;
                
                hoursToSaveFile += bb->in_hours_to_save_file;
            }
        }
    
    #endif
    
    
    //printOpenNodesList();
    //IQD_getchar();
    
    return 0;
}




void IQD_Callbacks:: newBestSolution( const int threadNumber, const double *newSol, const double oldBestObj, const double newBestObj, const long unsigned int iter )
{
    
    
    if( bb->in_recalculate_aux_bounds_on_zu_updt )
    {
        
        if(  (lastZuRecalcYBounds - newBestObj)/IQD_abs(newBestObj) >= bb->in_factor_to_recalculate_aux_bounds_on_zu_updt  )
        {
            const int n = refProb2->oprob->n;
            const int indvconstrs = refProb2->indvconstrs;
            const int dimv = refProb2->dimv;
            IQD_AuxBoundsUpdate &auxBoundsUpdt = auxBoundsUpdate[threadNumber];
            
            int ret;
            double *newly = values[threadNumber];
            double *newuy = &newly[dimv];
            
            
            
            
            //cout << "Entrei em IQD_Callbacks:: newBestSolution. threadNumber: " << threadNumber << " oldBestObj: " << oldBestObj << " newBestObj: " << newBestObj << " iter: " << iter << endl;
            
            
            ret = auxBoundsUpdt.calculateNewBounds(n, dimv, indvconstrs, refProb2->v, newBestObj, ly, uy, newly, newuy);
            
            if( ret != 0 )
            {
                if( bb->in_print_level > 0 )
                    IQD_PRINTERRORMSG("Failure to calculate new bounds to y");
            }
            
            
            /*for( int i = 0; i < dimv; i++ )
            {
                cout << "i: " << i << " ly: " << ly[i] << " uy: " << uy[i] << " newly: " << newly[i] << " newuy: " << newuy[i] << endl ;
            }*/
            
            
            SEMAPH_ybounds.lock(nthreads);
                
                if( lastZuRecalcYBounds > newBestObj ) //we can have several threads updating bounds in the same time...
                    lastZuRecalcYBounds = newBestObj;
            
                for( int i = 0; i < dimv; i++ )
                {
                    if( newly[i] > ly[i] )
                    {
                        ly[i] = newly[i];
                        yBoundsUpdated = true;
                    }
                    
                    if( newuy[i] < uy[i] )
                    {
                        uy[i] = newuy[i];
                        yBoundsUpdated = true;
                    }
                }
                
            SEMAPH_ybounds.unlock(nthreads);
            
            
            
            
            //( (OPT_LPSolver *) auxBoundsUpdt.solver )->generateModelFile( "iquadauxboundupdt.lp" );
            
            //getchar();
        }
    }
}

















