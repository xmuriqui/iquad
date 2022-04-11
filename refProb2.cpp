/**
* Author: Wendel Melo
* 
*/


#include <math.h>
#include <cstdlib>
#include <cassert>
#include <ctime>

#include <new>


#include "OPT_tools.hpp"

#include "iquad.hpp"
#include "IQD_tools.hpp"


#define IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME 100
#define IQD_N_TO_PRINT_MATRIX 20



using namespace minlpproblem;
using namespace optsolvers;
using namespace iquad;
using namespace std;




IQD_RefProb2NonLinearEval::IQD_RefProb2NonLinearEval(IQD_RefProb2 *refProb)
{
    initialize(refProb);
}


IQD_RefProb2NonLinearEval::~IQD_RefProb2NonLinearEval()
{
    desallocate();
}


int IQD_RefProb2NonLinearEval::allocate( const int nthreads)
{
    const int mo = refProb->oprob->m;
    
    
    //////////////////////////////////Begin of thread structures
    
    this->nthreads = nthreads;
    
    constrEvals = (bool **) calloc( nthreads, sizeof(bool*) );
    lambdas = (double**) calloc( nthreads, sizeof(double*) );
    
    if( !constrEvals || !lambdas )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        return IQD_MEMORY_ERROR;
    }
    
    
    for(int i = 0; i < nthreads; i++)
    {
        constrEvals[i] = (bool *) malloc( mo * sizeof(bool) );
        
        lambdas[i] = (double *) malloc( mo * sizeof(double) );
        
        if( !constrEvals[i] || !lambdas[i] )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            return IQD_MEMORY_ERROR;
        }
    }
    
    //////////////////////////////////End of thread structures
    
    
    
    {
        //getinf index of constraints that havd signal changed...
        const bool *chgSignal = refProb->chgSignal;
        const int mdis = refProb->mdis;
        
        
        nChgConstr = 0;
        for(int i = 0; i < mo + mdis; i++)
        {
            if( chgSignal[i] )
                nChgConstr++;
        }
        
        if( nChgConstr > 0 )
        {
            int mynChgConstr;
            
            
            indChgConstr = (int*) malloc( nChgConstr * sizeof(int) );
            if( !indChgConstr )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTMEMERROR;
                #endif
                return IQD_MEMORY_ERROR;
            }
            
            
            mynChgConstr = 0;
            for(int i = 0; i < mo + mdis; i++)
            {
                if( chgSignal[i] )
                {
                    indChgConstr[ mynChgConstr ] = i;
                    mynChgConstr++;
                }
            }
            
            #if IQD_DEBUG_MODE
                assert( mynChgConstr == nChgConstr );
            #endif
        }
    }
    
    return 0;
}



void IQD_RefProb2NonLinearEval::desallocate( )
{
    IQD_secFree(indChgConstr);
    
    
    //////////////////////////////////Begin of thread structures
    if( constrEvals )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if( constrEvals[i] )
                free(constrEvals[i]);
        }
        
        free(constrEvals);
        constrEvals = NULL;
    }
    
    if( lambdas )
    {
        for(int i = 0; i < nthreads; i++)
        {
            if( lambdas[i] )
                free(lambdas[i]);
        }
        
        free(lambdas);
        lambdas = NULL;
    }
    //////////////////////////////////End of thread structures
}



void IQD_RefProb2NonLinearEval::initialize(IQD_RefProb2 *refProb)
{
    nChgConstr = 0;
    indChgConstr = NULL;
    constrEvals = NULL;
    
    this->refProb = refProb;
}




int IQD_RefProb2NonLinearEval::eval_nl_obj_part( const int thnumber, const int n, const bool newx, const double* x, double& value)
{
    return refProb->oprob->nlEval->eval_nl_obj_part(thnumber, refProb->oprob->n, newx, x, value );
}


int IQD_RefProb2NonLinearEval::eval_nl_constrs_part( const int thnumber, const int n, const int m, const bool newx, const bool* constrEval, const double* x, double* values)
{
    const int mo = refProb->oprob->m;
    const int mdis = refProb->mdis;
    
    const int *origConstr = refProb->origConstr;
    
    bool *cEval = constrEvals[thnumber];
    
    int ret;
    
    
    IQD_copyArray( mo, constrEval, cEval);
    
    for(int i = mo; i < mo + mdis; i++)
    {
        if( constrEval[i] )
            cEval[ origConstr[i] ] = true;
    }
    
    
    ret = refProb->oprob->nlEval->eval_nl_constrs_part( thnumber, refProb->oprob->n, mo, newx, cEval, x, values);
    
    
    for(int i = mo; i < mo + mdis; i++)
        values[i] = values[ origConstr[i] ];
    
    
    for(int i = 0; i < nChgConstr; i++)
        values[ indChgConstr[i] ] *= -1.0;
    
    
    return ret;
}


int IQD_RefProb2NonLinearEval::eval_grad_nl_obj_part( const int thnumber, const int n, const bool newx, const double* x, double* values)
{
    return refProb->oprob->nlEval->eval_grad_nl_obj_part(thnumber, refProb->oprob->n, newx, x, values);
}


int IQD_RefProb2NonLinearEval::eval_grad_nl_constrs_part( const int thnumber, const int n, const int m, const int nz, const bool newx, const bool* constrEval, const double* x, MIP_SparseMatrix& jacobian)
{
    const int mo = refProb->oprob->m;
    const int mdis = refProb->mdis;
    
    const int *origConstr = refProb->origConstr;
    
    bool *cEval = constrEvals[thnumber];
    IQD_IQuadProb *oprob = refProb->oprob;
    
    int ret;
    
    
    IQD_copyArray( mo, constrEval, cEval);
    
    for(int i = mo; i < mo + mdis; i++)
    {
        if( constrEval[i] )
            cEval[ origConstr[i] ] = true;
    }
    
    
    ret = oprob->nlEval->eval_grad_nl_constrs_part( thnumber, oprob->n, mo, oprob->J.getNumberOfElements(), newx, cEval, x, jacobian );
    
    
    for(int i = mo; i < mo + mdis; i++)
        jacobian.setRowValues(i, jacobian(origConstr[i]) ); //jacobian[i].copyValuesInOrderFrom( jacobian[origConstr[i]] );
    
    
    for(int i = 0; i < nChgConstr; i++)
        jacobian.multiplyAllElementsAtRow(indChgConstr[i], -1.0); //jacobian[ indChgConstr[i] ].multiplyAllElements(-1.0);
    
    
    return ret;
}


int IQD_RefProb2NonLinearEval::eval_hessian_nl_lagran_part( const int thnumber, const int n, const int m, const int nz, const bool newx, const double* x, const double objFactor, const double* lambda, MIP_SparseMatrix& hessian)
{
    const int mo = refProb->oprob->m;
    const int mdis = refProb->mdis;
    
    double *mylambda = lambdas[thnumber];
    const int *origConstr = refProb->origConstr;
    const bool *chgSignal = refProb->chgSignal;
    
    
    
    IQD_copyArray( mo, lambda, mylambda );
    
    for( int i = 0; i < nChgConstr; i++ )
    {
        if( indChgConstr[i] < mo )
            mylambda[ indChgConstr[i] ] *= -1.0;
        else
            break; //index are in order...
    }
    
    
    //now, we run indices from dismembered constraints...
    for(int i = mo; i < mo + mdis; i++)
    {
        const int ind = origConstr[i];
        
        if( chgSignal[i] )
            mylambda[ ind ] -= lambda[i];
        else
            mylambda[ ind ] += lambda[i];
    }
    
    
    return refProb->oprob->nlEval->eval_hessian_nl_lagran_part( thnumber, refProb->oprob->n, mo, nz, newx, x, objFactor, mylambda, hessian );
    
}







IQD_RefProb2::IQD_RefProb2()
{
    initialize();
}


IQD_RefProb2::~IQD_RefProb2()
{
    desallocate();
}


void IQD_RefProb2::initialize()
{
    nonConvexConstrs = false;
    dimv = 0;
    mdis = 0;
    chgSignal = NULL;
    origConstr = NULL;
    lambdaObj = NULL;
    lambdaObjNoInt = NULL;
    lambdaConstr = NULL;
    PnoInt = NULL;
    //v = NULL;
    ly = NULL;
    uy = NULL;
    dimLambdaConstr = 0;
    eps_diag_matrix_from_sdp = 0.0;
    abs_eps_diag_P = 0.0;
    rel_eps_diag_P = 0.0;
    zero_tol = 0.0;
    max_cpu_time = INFINITY;
    max_time = INFINITY;
}



int IQD_RefProb2::allocateOrigConstrArrays(const int size)
{
    origConstr = (int *) malloc( size* sizeof(int) );
    chgSignal = (bool *) calloc( size, sizeof(bool) );
    
    if(!origConstr || !chgSignal)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        return IQD_MEMORY_ERROR;
    }
    
    return 0;
}



int IQD_RefProb2::allocateVAndLambda( const int nconstr, const int dimv, const int cols )
{
    int r;
    
    lambdaObj = (double *) malloc( dimv * sizeof(double) );
    
    ly = (double *) malloc( 2*dimv * sizeof(double) );
    
    if( !lambdaObj || !ly)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        return IQD_MEMORY_ERROR;
    }
    
    uy = &ly[dimv];
    
    
    
    if( nconstr > 0)
    {
        lambdaConstr = (double **) calloc( nconstr, sizeof(double *) );
        if( !lambdaConstr )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            return IQD_MEMORY_ERROR;
        }
        
        /*lambdaConstr[0] = (double *) malloc( nconstr * dimv * sizeof(double) );
        if( !lambdaConstr[0] )
            return IQD_MEMORY_ERROR;
        
        for(int i = 1; i < nconstr; i++)
            lambdaConstr[i] = &lambdaConstr[0][i*dimv]; */
        
        dimLambdaConstr = nconstr;
    }
    
    
    
    r = v.addNewRows( dimv );
    
    if( r != 0 )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        return IQD_MEMORY_ERROR;
    }
    
    v.setNumberOfColumns(cols);
    
    
    /*v = (double **) malloc( dimv * sizeof(double *) );
    
    if(!v)
        return IQD_MEMORY_ERROR;
    
    //v[0] = (double *) malloc( dim * dim * sizeof(double) );
    //we allocate memory in auxv instead of v[0] because we can reorder lines in v, for example in MILP basis problem...
    auxv = (double*) malloc( dimv*cols*sizeof(double) );
    
    if(!auxv)
        return IQD_MEMORY_ERROR;
    
    for(int i = 0; i < dimv; i++)
        v[i] = &(auxv[cols*i]) ; */
    
    
    
    this->dimv = dimv;
    
    return 0;
}





void IQD_RefProb2::desallocate()
{
    //IQD_secFree( auxv );
    //IQD_secFree( v );
    
    MIP_NonLinearEval *eval = rprob.getNonLinearEvaluationObject();
    
    if( eval )
    {
        delete eval;
        rprob.setNonLinearEvaluationObject(NULL); //it is not necessary but...
    }
    
    rprob.deallocateMatrices();
    
    
    v.desallocateMemory();
    
    IQD_secFree( chgSignal );
    IQD_secFree( origConstr );
    IQD_secFree( lambdaObj );
    IQD_secFree( lambdaObjNoInt );
    
    if( lambdaConstr )
    {
        for(int i = 0; i < dimLambdaConstr; i++)
        {
            if(lambdaConstr[i])	
                free(lambdaConstr[i]);
        }
        
        
        free(lambdaConstr);
        lambdaConstr = NULL;
    }
    
    IQD_secDelete( PnoInt );
    
    IQD_secFree(ly);
}




//do not use to split by shcur...
int IQD_RefProb2::splitMatrix(const int splitWay, IQD_SparseMatrix &Q, const double factor, const int nvecs, double **vecs, double *vecsWeight, int sdpSolver, IQD_GeneralSolverParams* sdpSolverParams, double *lambda, double *R, const bool initializeRwithZeros )
{
    const int n = Q.getNumberOfRows();
    int ret, code;
    
    IQD_SDPLambda sdpLambda;
    
    
    
    if( R && initializeRwithZeros )
    {
        const int nn = n*n;
        
        for( int i = 0; i < nn; i++ )
            R[i] = 0.0;
    }
    
    
    if( splitWay == IQD_SQ_MIX_SCHUR_DIAG_SDP )
    {
        
        ret = sdpLambda.calculateLambda( sdpSolver, Q, factor, nvecs, vecs, lambda, vecsWeight, sdpSolverParams, zero_tol);//, eps_diag_matrix_from_sdp );
        
        if( ret != 0 )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTERRORNUMBER(ret);
            #endif
            code = ret;
            goto termination;
        }
        
        
        if( R )
            calculateRbyLambdaAndv( n, nvecs, vecs, lambda, R, false );
        
    }
    else
    {
        //strict diagonal splits...
        
        //we assume the first n columns in vecs are the eye matrix..
            
        #if IQD_DEBUG_MODE
            for(int i = 0; i < n; i++)
                assert( vecs[i][i] == 1.0 );
        #endif
        
        for( int i = n; i < nvecs; i++ )
            lambda[i] = 0.0;
        
        
        if( splitWay == IQD_SQ_DIAG_SDP )
        {
            ret = sdpLambda.calculateLambda( sdpSolver, Q, factor, n, vecs, lambda, NULL, sdpSolverParams, zero_tol);//, eps_diag_matrix_from_sdp );
            
        }
        else if( splitWay == IQD_SQ_IDENTITY )
        {
            ret = calculateLambdaByIdentity( Q, factor, lambda );
            
        }
        else if( splitWay == IQD_SQ_DIAG_DOM )
        {
            ret = calculateLambdaByDiagDominant(Q, factor, 0.0, lambda);
        }
        else if( splitWay == IQD_SQ_SCHUR )
        {
            cerr << "iquad: invalid split way" << IQD_GETFILELINE << " Check your parameter" << endl;
            
            assert(false); //do not calculate schur split by this function!
            ret = IQD_BAD_DEFINITIONS;
        }
        else
        {
            cerr << "iquad: invalid split way" << IQD_GETFILELINE << " Check your parameter" << endl;
            
            assert(false);
            ret = IQD_BAD_DEFINITIONS;
        }
        
        if( ret != 0 )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTERRORNUMBER(ret);
            #endif
            code = ret;
            goto termination;
        }
        
        
        if(R)
        {
            for(int i = 0; i < n; i++)
            {
                if( lambda[i] < -zero_tol )
                R[ i*n + i ] = lambda[i];
            }
        }
        
    }
    
    
    /*if( eps_diag_P > 0.0 )
    {
        for(int i = 0; i < n; i++)
        {
            const int diag = i*n + i;
            if( IQD_abs( R[diag] ) > zero_tol )
                R[diag] -= eps_diag_P; //we subtrate from R to add on p;
        }
    } */
    
    
    
    code  = 0;
    
termination:
    
    return code;
}



int IQD_RefProb2::calculateBoundsOnDirections( const double max_cpu_time, const double max_time, const int qpSolver, const int qcpSolver, const int nlpSolver, IQD_GeneralSolverParams* solverParams, const int ndirs, double** dirs, bool* lbconc, bool* ubconv, double* lbdir, double* ubdir )
{
    const int n = oprob->getNumberOfVars();
    const int m = oprob->getNumberOfConstraints();
    
    MIP_SparseMatrix &A = oprob->A;
    MIP_SparseMatrix *QC = oprob->QC;
    
    bool nlConstr = false;
    int r, code = 0, maxnzs = 0;
    double lb, ub;
    int *rows = NULL, *cols = NULL;
    double *values = NULL;
    OPT_LPSolver *solver;
    
    clock_t clockStart;// = clock();
    double timeStart;// = IQD_getTime();
    
    
    if( max_cpu_time < INFINITY )
        clockStart = clock();
    
    if( max_time < INFINITY )
        timeStart = IQD_getTime();
    
    
    r = oprob->getProblemType( ) ;
    
    if( MIP_isNonlinearProblemType(r) )
    {
        solver = OPT_newNLPSolver(nlpSolver);
        maxnzs = n;
    }
    else if( MIP_isQuadConstrProblemType(r) )
        solver = OPT_newQCPSolver(qcpSolver);
    else
    {
        maxnzs = -1;
        solver = OPT_newQPSolver(qpSolver);
    }
    
    
    if( !solver )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    
    if( maxnzs >= 0 )
    {
        
        for(int i = 0; i < m; i++)
        {
            if( lbconc[i] || ubconv[i] )
            {
                int nzs = QC[i].getNumberOfElements();
                
                if( nzs > maxnzs )
                    maxnzs = nzs;
            }
        }
        
        rows = (int*) malloc(maxnzs* sizeof(int));
        cols = (int*) malloc(maxnzs* sizeof(int));
        values = (double*) malloc( maxnzs* sizeof(double) );
        
        if( !rows || !cols || !values )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
    }
    
    
    
    r = solver->initSolverEnv( n, m );
    if( r != 0 )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTERRORNUMBER(r);
        #endif
        code = IQD_SOLVER_ERROR;
        goto termination;
    }
    
    
    r = solver->addVariables(n) + solver->addConstraints(m);
    if( r != 0 )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTERRORNUMBER(r);
        #endif
        code = IQD_SOLVER_ERROR;
        goto termination;
    }
    
    
    
    for( int i = 0; i < n; i++ )
    {
        oprob->getVariableBounds(i, lb, ub);
        
        r = solver->setVariableBounds(i, lb, ub);
        #if IQD_DEBUG_MODE
            if( r != 0 )
            {
                IQD_PRINTERRORNUMBER(r);
                code = IQD_SOLVER_ERROR;
                goto termination;
            }
        #endif
    }
    
    
    for( int i = 0; i < m; i++ )
    {
        
        if( lbconc[i] || ubconv[i] )
        {
            r = solver->resetConstraintLinearPart(i, i, A); // solver->resetLinearConstraintPart( i, A[i]);
            
            if( r != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(r);
                #endif
                
                code = IQD_SOLVER_ERROR;
                goto termination;
            }
        
            oprob->getConstraintBounds(i, lb, ub);
            
            int nzs = QC[i].getNumberOfElements();
            
            if( nzs > 0 )
            {
                //QC[i].getStructureAndValues(rows, cols, values);
                
                
                //r = ( (OPT_QCPSolver *) solver)->setQuadConstraintMatrix( i, nzs, rows, cols, values );
                r = ( (OPT_QCPSolver *) solver)->setConstraintQuadMatrix(i, QC[i]);
                
                if( r != 0 )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTERRORNUMBER(r);
                    #endif
                    code = IQD_SOLVER_ERROR;
                    goto termination;
                }
                
            }
            
            if( oprob->hasConstraintNLTerm(i) )
            {
                OPT_NLPSolver *nlpsolver = (OPT_NLPSolver *) solver;
                
                nlpsolver->setNonLinearEvalObject( oprob->nlEval );
                
                nlpsolver->setConstrNLFlag(i, true );
                
                
                oprob->getJacobianRowStructure(i, &nzs, cols);
                
                r = nlpsolver-> setJacobianRowStructure( i, nzs, cols);
                
                if( r != 0 )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTERRORNUMBER(r);
                    #endif
                    code = IQD_SOLVER_ERROR;
                    goto termination;
                }
                
                
                nlConstr = true;
            }
            
            r = solver->setConstraintBounds( i, lbconc[i] ? lb : -OPT_INFINITY, ubconv[i] ? ub : OPT_INFINITY );
        }
        
    }
    
    
    
    if( nlConstr )
    {
        OPT_NLPSolver *nlpsolver = (OPT_NLPSolver *) solver;
        
        int nzs;
        
        for( int i = 0; i < n; i++ )
        {
            oprob->getLagrangianHessianRowStructure(i, &nzs, cols);
            
            r = nlpsolver->setLagrangianHessianRowStructure( i, nzs, cols);
            
            if( r != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(r);
                #endif
                code = IQD_SOLVER_ERROR;
                goto termination;
            }
        }
        
    }
    
    
    if(solverParams)
    {
        r = solver->setParameters( *solverParams );
        if( r != 0 )
        {
            IQD_PRINTERRORMSG( "Error at setting solver parameter");
        }
    }
    
    
    
    for(int i = 0; i< ndirs; i++)
    {
        r = solver->setObjLinearPart(n, dirs[i]);
        #if IQD_DEBUG_MODE
            if( r != 0 )
            {
                IQD_PRINTERRORNUMBER(r);
                code = IQD_SOLVER_ERROR;
                goto termination;
            }
        #endif
        
        
        solver->setObjSense( OPT_MINIMIZE );
        
        //solver->generateModelFile( "iquadbounds.lp" );
        //IQD_getchar();
        
        r = solver->solve();
        
        if( r == OPT_OPTIMAL_SOLUTION )
        {
            lbdir[i] = solver->objValue;
        }
        else
        {
            IQD_PRINTERRORMSG("Error at calculating bounds to spatial branch-and-bound");
            
            IQD_PRINTERRORNUMBER(r);
            
            if( r == OPT_INFEASIBLE_PROBLEM )
            {
                code = IQD_INFEASIBLE_PROBLEM;
                goto termination;
            }
            else if( r == OPT_UNBOUNDED_PROBLEM )
            {
                //code = IQD_UNBOUNDED_PROBLEM;
                
                code = IQD_NO_BOUNDED_FEASIBLE_SET;
                lbdir[i] = -IQD_INFINITY_TO_DIRECTIONS_BOUNDS;
                //goto termination;
            }
            else
            {
                code = IQD_SOLVER_ERROR;
                //goto termination;
            }
        }
        
        
        solver->setObjSense( OPT_MAXIMIZE );
        
        r = solver->solve();
        
        if( r == OPT_OPTIMAL_SOLUTION )
        {
            ubdir[i] = solver->objValue;
        }
        else
        {
            IQD_PRINTERRORMSG("Error at calculating bounds to spatial branch-and-bound");
            
            IQD_PRINTERRORNUMBER(r);
            
            //mosek declares a problem like dual infeasible if primal is infeasible or unbounded...
            if( r == OPT_UNBOUNDED_PROBLEM || r == OPT_INFEASIBLE_PROBLEM )
            {
                
                //code = IQD_UNBOUNDED_PROBLEM;
                code = IQD_NO_BOUNDED_FEASIBLE_SET;
                ubdir[i] = IQD_INFINITY_TO_DIRECTIONS_BOUNDS;
                //goto termination;
            }
            else
            {
                code = IQD_SOLVER_ERROR;
                //goto termination;
            }
        }
        
        
        if( max_cpu_time < INFINITY )
        {
            if( (double(clock()-clockStart) )/CLOCKS_PER_SEC >= max_cpu_time  )
            {
                code = IQD_MAX_TIME_STOP;
                goto termination;
            }
        }
        
        if( max_time < INFINITY )
        {
            if( IQD_getTime() - timeStart >= max_time )
            {
                code = IQD_MAX_TIME_STOP;
                goto termination;
            }
        }
        
    }
    
    
    
    
termination:
    
    if( solver )	delete solver;
    if( rows )		free(rows);
    if( cols )		free(cols);
    if( values )	free(values);
    
    return code;
}



int IQD_RefProb2::calculateLambdaByIdentity( IQD_SparseMatrix& Q, const double factor, double *lambda )
{
    const int n = Q.getNumberOfRows();
    int code;
    double lam;
    
    
    code = IQD_getEigenValues(Q, lambda);
    
    
    //lambda[0] has the lowest eigenvalue (negative)... But, if factor is negative, it will be the greatest
    if( factor > 0.0 )
    {
        lam = lambda[0]*factor;
    }
    else
    {
        //eigenvalues are ordered and we have a negative factor. So, the greatest now is the lowest
        lam = lambda[n-1]*factor;
    }
    
    //giving a tolerance... we use zero tolerance as relative tolerance too
    lam += IQD_min( -zero_tol, zero_tol*lam ); //P = Q - R, so we let the minus signal...
    
    for(int i = 0; i < n; i++)
        lambda[i] = lam;
    
    
    return code;
}



int IQD_RefProb2::calculateLambdaByDiagDominant(IQD_SparseMatrix& Q, const double factor, const double epsilon, double *lambda)
{
    const int n = Q.getNumberOfRows();
    const double afactor = IQD_abs(factor);
    
    int code;
    double *diag = NULL;
    
    
    diag = (double *) calloc( n, sizeof(double) );
    if(!diag)
    {
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    for(int i = 0; i < n; i++)
        lambda[i] = epsilon;
    
    
    for(int i = 0; i < n; i++ )
    {
        //IQD_SparseRow &myrow = Q[i];
        
        const int *rcols = Q[i];
        const double *rvals = Q(i);
        const int nel = Q.getNumberOfElementsAtRow(i); //myrow.getNumberOfElements();
        
        for(int j = 0; j < nel; j++)
        {
            const int col = rcols[j]; //myrow[j].getColumn();
            const double v = rvals[j]; //myrow[j].getValue();
            
            if( i == col )
            {
                diag[i] = v;
            }
            else
            {
                const double av = IQD_abs(v);
                lambda[i] += av;
                lambda[col] += av;
            }
        }
    }
    
    for(int i = 0; i < n; i++)
    {
        const double di = diag[i]*factor;
        const double li =  lambda[i]*afactor;
        
        if( di >= li )
            lambda[i] = 0.0;
        else
            lambda[i] = di - li; 
    }
    
    
    
    code = 0;
    
termination:
    
    if(diag)	free(diag);
    
    return code;
}



void IQD_RefProb2::calculateRbyLambdaAndv(const int n, const int dimv, double **v, const double *lambda, double *R, const bool initializeRwithZeros)
{
    int i, j, k;
    
    
    if( initializeRwithZeros )
    {
        
        IQD_setAllarray(n*n, R, 0.0);
    }
    
    
    for(k = 0; k < dimv; k++)
    {
        if( lambda[k] < -zero_tol ) 
        {
            const double *vk = v[k];
            
            for(i = 0; i < n; i++)
            {
                double *ri = &R[i*n];
                
                const double lvki = lambda[k] * vk[i];
                
                for(j = 0; j <= i; j++) //strict lower triangle
                {
                    //aux = lvki*vk[j];
                    ri[j] += lvki *vk[j]; //aux;
                    //R[j*n + i] += aux;
                }
                
                //ri[i] += lvki *vk[i]; //diagonal position
            }
        }
    }
    
    
    //setting upper triangle
    for( i = 1; i < n; i++ )
    {
        const double *ri = &R[i*n];
        
        for(j = 0; j < i; j++)
        {
            R[j*n + i] = ri[j];
        }
    }
    
    
}




void IQD_RefProb2::calculatey(const double *x, double *y)
{
    const int dimv = this->dimv;
    
    for(int i = 0; i < dimv; i++)
    {
        //y[i] = 0.0;
        
        //if( lambda[i] >= -zero_tol )
            //continue;
        
        //for(j = 0; j < n; j++)
            //y[i] += v[i][j]* x[j];
        
        y[i] = v.rowEvaluation(i, x); //v[i].evalTimesxt( x );
    }
}


static inline void IQD_correctBoundDirToInt( const double intTol, double &lbdir, double &ubdir )
{
    const double rolbdir = round(lbdir), roubdir = round(ubdir);
    
    if( IQD_abs(lbdir - rolbdir) <= intTol )
        lbdir = rolbdir;
    if( IQD_abs(ubdir - roubdir) <= intTol )
        ubdir = roubdir;
}


int IQD_RefProb2::buildRefProb(bool splitWithouIntVars, const bool certainOfNonconvexQ, const int sdpSolver, const int qpSolver, const int qcpSolver, const int nlpSolver,  IQD_GeneralSolverParams* sdpSolverParams, IQD_GeneralSolverParams* solverParams)
{
    const int n = oprob->getNumberOfVars();
    const int m = oprob->getNumberOfConstraints();
    
    const double *lc = oprob->lc;
    const double *uc = oprob->uc;
    
    const bool storeSchurVecs = !IQD_isDiagonalSplitQWay(splitWay);
    bool useEyeVector;
    int newm = 0, cnewm = 0;
    
    IQD_SparseMatrix &Q = oprob->Q;
    IQD_SparseMatrix &A = oprob->A;
    IQD_SparseMatrix *Qnoint = NULL;
    
    
    int splitWay = this->splitWay;
    int nlamobj = 0;
    int aux, ndirs, nnewdirs, k, ret, code;
    int nI, nauxconstrs = 0;
    //double lb, ub;
    
    bool *lbconcave = NULL, *ubconvex = NULL;
    int *intVars = NULL;
    double *M = NULL;
    double *lambda = NULL, **schvecs = NULL;
    double *eye = NULL, **dirs = NULL;
    double *lbdir = NULL, *ubdir, *vecsWeigth = NULL;
    //double *lamObj = NULL, **lamConstr;
    
    int *cols =  NULL; 
    double *values = NULL;
    
    clock_t clockStart = clock();
    double timeStart = IQD_getTime();
    
    
    
    
    if( splitWithouIntVars )
    {
        nI = oprob->getNumberOfIntegerVars();
        
        if( nI > 0 )
        {
            bool intOnQ = false;
            
            intVars = (int*) malloc(nI * sizeof(int));
            if( !intVars )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTMEMERROR;
                #endif
                
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            oprob->getIntegerIndices(intVars);
            
            
            //checking if some intever var appears in quadratic part of objective function. If not, we do not need split Q without those vars...
            //first, we check if some row about integer vars is not empty... remeber Q only store lower triangle
            for(int i = 0; i < nI; i++)
            {
                if(Q.getNumberOfElementsAtRow(intVars[i]	) > 0) //if( Q[ intVars[i] ].getNumberOfElements() > 0 )
                {
                    intOnQ = true;
                    break;
                }
            }
            
            //now, we have to check about upper triangle
            if( !intOnQ )
            {
                for( int i = 0; i < nI; i++)
                {
                    const int ind = intVars[i];
                    
                    //we just need chek rows below ind (Q only store lower triangle...)
                    for(int j = ind+1; j < n; j++)
                    {
                        //MIP_SparseRow &row = Q[j];
                        const int* rcols = Q[j];
                        const unsigned int nel = Q.getNumberOfElementsAtRow(j); // row.getNumberOfElements();
                        
                        for(unsigned int k = 0; k < nel; k++)
                        {
                            if(rcols[k] == ind) //if( (int) row[k].getColumn() == ind )
                            {
                                intOnQ = true;
                                //we want stop all fors...
                                j = n;
                                i = n;
                                break;
                            }
                        }
                    }
                }
            }
            
            if( !intOnQ )
                splitWithouIntVars = false;
            
        }
        else
            splitWithouIntVars = false;
    }
    
    
    lbconcave = (bool *) malloc( 2*m * sizeof(bool) );
    M = (double *) malloc( n*n * sizeof(double) );
    lambda = (double *) malloc(n* sizeof(double));
    
    if( !lbconcave || !M || !lambda )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    ubconvex = &lbconcave[m];
    
    
    
    
    
    
    nonConvexConstrs = false;
    
    for(int i = 0; i < m; i++ )
    {
        IQD_SparseMatrix &QCi = oprob->QC[i];
        
        
        lbconcave[i] = true; 
        ubconvex[i] = true;
        
        if( QCi.getNumberOfElements() > 0 )
        {
            #if IQD_DEBUG_MODE
                if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
                    std::cout << IQD_PREPRINT "Checking eigenvalues of matrix Q in constraint " << i << " ..." << std::endl;
            #endif
            
            //oprob->getConstraintBounds(i, lb, ub);
            
            
            ret = IQD_getEigenValues(QCi, lambda);
            
            
            if( uc[i] < MIP_INFINITY )
            {
                for(int j = 0; j < n; j++)
                {
                    if( lambda[j] < -zero_tol )
                    {
                        ubconvex[i] = false;
                        break;
                    }
                }
            }
            
            if( lc[i] > -MIP_INFINITY )
            {
                //I think eigenvalues are in ascendig order (but, I am not sure). So, we start from the end to try save computation...
                for(int j = n-1; j >= 0; j--)
                {
                    if( lambda[j] > zero_tol )
                    {
                        lbconcave[i] = false;
                        break;
                    }
                }
            }
            
            
            /*for(int j = 0; j < n; j++)
            {
                if( lambda[j] < -zero_tol && ub < MIP_INFINITY )
                    ubconv[i] = false;
                else if( lambda[j] > zero_tol && lb > -MIP_INFINITY )
                    lbconc[i] = false;
            } */
            
            
            
            if(lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY)
                newm++;
            
            
            if( (!lbconcave[i] && lc[i] > -MIP_INFINITY) || ( !ubconvex[i] && uc[i] < MIP_INFINITY) ) //if( (!lbconc[i] || !ubconv[i]) && (lb > -MIP_INFINITY || ub < MIP_INFINITY) )
                nonConvexConstrs = true;
            
            
            
            if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
            {
                if( max_cpu_time < INFINITY )
                {
                    if( (double(clock()-clockStart) )/CLOCKS_PER_SEC > max_cpu_time  )
                    {
                        code = IQD_MAX_TIME_STOP;
                        goto termination;
                    }
                }
                
                if( max_time < INFINITY )
                {
                    if( IQD_getTime() - timeStart > max_time )
                    {
                        code = IQD_MAX_TIME_STOP;
                        goto termination;
                    }
                }
            }
            
        }
        
        //we assume nonlinear part is convex...
    }
    
    
    //for(int i = 0; i < m; i++ )
        //std::cout << "lbconc["<<i<<"]: " << lbconc[i] << " ubconv["<<i<<"]: " << ubconv[i] << "\n" ;
    
    //std::cout << "nonConvexConstrs: " << nonConvexConstrs << "\n";
    
    
    useEyeVector = nonConvexConstrs || splitWithouIntVars || splitWay != IQD_SQ_SCHUR;
    
    dirs = (double **) calloc( (useEyeVector && storeSchurVecs? 2*n : n), sizeof(double*) );
    lbdir = (double *) malloc( (useEyeVector && storeSchurVecs? 4 : 2)*n* sizeof(double) );
    if(!dirs || !lbdir)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    ubdir = &lbdir[ (useEyeVector && storeSchurVecs? 2*n : n) ];
    
    
    
    if( useEyeVector )
    {
        eye = (double *) calloc( n*n, sizeof(double));
        
        if(!eye)
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        for(int i = 0; i < n; i++)
            eye[ i*n + i ] = 1.0;
        
        
        //setting the first n vector to do decomposition with eye
        for(int i = 0; i < n; i++)
            dirs[i] = &eye[i*n];
    }
    
    
    
    
    
    
    
    if( Q.getNumberOfElements() > 0 && !certainOfNonconvexQ )
    {
        const int start = useEyeVector ? n : 0;
        
        schvecs = (double **) calloc( n, sizeof(double*));
        if( !schvecs )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        
        schvecs[0] = (double*) malloc(n*n* sizeof(double));
        if( !schvecs[0] )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        for(int i = 1; i < n; i++)
            schvecs[i] = &schvecs[0][i*n];
        
        
        
        Q.copyMatrixTo(M, true, true, oprob->objFactor);
        
        #if IQD_DEBUG_MODE
            if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
                std::cout << IQD_PREPRINT << "Calculating schur decomposition on Q...\n";
        #endif
        
        /*std::cout << "calculando auto-valores..." << std::endl;
        IQD_getEigenValues(n, M, lambda);
        std::cout << "agora sim..." << std::endl; */
        
        ret = IQD_schurDecomposition(M, n, M, schvecs[0], lambda );
        
        if( ret != 0 )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTERRORNUMBER(ret);
            #endif
            code = ret;
            goto termination;
        }
        
        #if IQD_DEBUG_MODE
            if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
                std::cout << IQD_PREPRINT << "Schur decomposition on Q calulated \n";
        #endif
        
        
        /*cout << endl << "lambda obj schur: " ;
        for(int i = 0; i < n; i++)
            cout << lambda[i] << " ";
        cout << endl;
        
        cout << "schur vectors: " << endl;
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
                cout << schvecs[j][i] << " " ;
            cout << endl;
        }
        cout << endl; */
        
        
        for(int i = 0; i < n; i++)
        {
            if( lambda[i] < -zero_tol )
            {
                if( storeSchurVecs )
                    dirs[start + nlamobj] = schvecs[i];
                
                nlamobj++;
            }
        }
        
    }
    
    /*cout << endl << "nlobj: " << nlobj << endl;
    
    cout << "directions: " << endl ;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n + nlobj; j++)
            cout << dirs[j][i] << " ";
        cout << endl;
    }
    cout << endl; */
    
    
    //ndirs = !nonConvexConstrs && !splitWithouIntVars && splitWay == IQD_SQ_SCHUR ?  nlamobj :  n + nlamobj;
    
    ndirs = (useEyeVector ? n : 0) + (storeSchurVecs ? nlamobj : 0); 
    
    /*std::cout << "ndirs: " << ndirs << std::endl;
    for(int i = 0; i < ndirs; i++)
    {
        for(int j = 0; j < n; j++)
            std::cout << dirs[i][j] << " ";
        std::cout << std::endl;
    }
    IQD_getchar(); */
    
    
    #if IQD_DEBUG_MODE
        assert( ndirs > 0 || nlamobj == 0 );
    #endif
    
    rprob.deallocateMatrices();
    
    ret = rprob.addVariables(n);
    ret += rprob.addConstraints(m + newm);
    
    if( ret != 0 )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    
    if( oprob->hasLinCoefObj() )
    {
        const double *c = oprob->c;
        const double of = oprob->objFactor;
        
        for( int i = 0; i < n; i++)
            rprob.setObjLinearCoefficient(i, of*c[i]);
    }
    
    rprob.setObjConstant(oprob->d * oprob->objFactor);
    
    
    if( oprob->hasObjNLTerm() || oprob->hasNLConstraints() )
    {
        #if IQD_DEBUG_MODE
            assert( oprob->nlEval != NULL );
        #endif
        
        
        
        if( oprob->hasObjNLTerm() )
            rprob.setObjNonLinearTermFlag( true );
        
        
        if( oprob->getNumberOfJacobianNonZeros() > 0 )
        {
            int r, nz;
            for( int i = 0; i < m; i++ )
            {
                r = oprob->getNumberOfJacobianNonZerosAtRow(i, nz);
                
                #if IQD_DEBUG_MODE
                    if( r != 0 )
                    {
                        IQD_PRINTERRORNUMBER(r);
                        code = IQD_SPLITING_ERROR;
                        goto termination;
                    }
                #endif
                
                if( nz > 0 )
                {
                    //r = rprob.setJacobianRowStructure(i, oprob->J[i]);
                    r = rprob.setJacobianRowStructure(i, oprob->J.getNumberOfElementsAtRow(i), oprob->J[i]);
                    if( r != 0 )
                    {
                        #if IQD_DEBUG_MODE
                            IQD_PRINTERRORNUMBER(r);
                        #endif
                        code = r == minlpproblem::MIP_MEMORY_ERROR ? IQD_MEMORY_ERROR : IQD_SPLITING_ERROR;
                        goto termination;
                    }
                }
            }
        }
        
        
        
        if( oprob->getNumberOfLagrangianHessianNonZeros() > 0 )
        {
            int r, nz;
            
            for( int i = 0; i < n; i++ )
            {
                r = oprob->getNumberOfLagrangianHessianNonZerosAtRow(i, nz);
                
                #if IQD_DEBUG_MODE
                    if( r != 0 )
                    {
                        IQD_PRINTERRORNUMBER(r);
                        code = IQD_SPLITING_ERROR;
                        goto termination;
                    }
                #endif
                
                if( nz > 0 )
                {
                    //r = rprob.setLagrangianHessianRowStructure( i, oprob->lagH[i] );
                    r = rprob.setLagrangianHessianRowStructure( i, oprob->lagH.getNumberOfElementsAtRow(i), oprob->lagH[i] );
                    if( r != 0 )
                    {
                        #if IQD_DEBUG_MODE
                            IQD_PRINTERRORNUMBER(r);
                        #endif
                        code = r == minlpproblem::MIP_MEMORY_ERROR ? IQD_MEMORY_ERROR : IQD_SPLITING_ERROR;
                        goto termination;
                    }
                    
                }
            }
        }
        
        
        
        IQD_RefProb2NonLinearEval *nlEval = new (std::nothrow) IQD_RefProb2NonLinearEval(this);
        
        if( !nlEval )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        
        rprob.setNonLinearEvaluationObject( nlEval );
    }
    
    
    
    
    
    for( int i = 0; i < n; i++ )
        rprob.setVariableLowerBound(i, oprob->lx[i]);
    
    for( int i = 0; i < n; i++ )
        rprob.setVariableUpperBound(i, oprob->ux[i]);
    
    for( int i = 0; i < n; i++ )
    {
        MIP_VARTYPE vt = oprob->xtype[i] == MIP_VT_INTEGER ? MIP_VT_INTEGER : MIP_VT_CONTINUOUS;
        
        rprob.setVariableType(i, vt);
    }
    
    
    
    
    ret = allocateVAndLambda(m + newm, ndirs, n);
    ret+= allocateOrigConstrArrays(m + newm);
    if( ret != 0 )
    {
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    IQD_setAllarray( m+newm, origConstr, -1 );
    
    /*if( nlamobj > 0 )
    {
        k = nonConvexConstrs || splitWay == IQD_SQ_MIX_SCHUR_DIAG_SDP ? n : 0;
        for(int i = 0; i < n; i++)
        {
            if( lambda[i] < -zero_tol )
            {
                lambdaObj[k] = lambda[i];
                k++;
            }
        }
        
        #if IQD_DEBUG_MODE
            assert( k - n == nlamobj );
        #endif
    } */
    
    #if IQD_DEBUG_MODE
        if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
            std::cout << IQD_PREPRINT << "Calculating bounds on directions..." << std::endl;
    #endif
    
    ret = calculateBoundsOnDirections( max_cpu_time - (double(clock()-clockStart) )/CLOCKS_PER_SEC , max_time - (IQD_getTime() - timeStart), qpSolver, qcpSolver, nlpSolver, solverParams, ndirs, dirs, lbconcave, ubconvex, lbdir, ubdir );
    
    if( ret != 0 && ret != IQD_NO_BOUNDED_FEASIBLE_SET )
    {
        code = ret;
        goto termination;
    }
    
    
    if( oprob->getNumberOfIntegerVars() > 0 )
    {
        //we try to correct some numerical problem on integer bounds to directions
        const double intTol = 1.0e-4;
        const int *xtype = oprob->xtype;
        int neyedirs = 0;
        
        
        if( useEyeVector )
        {
            for(int i = 0; i < n; i++)
            {
                if( MIP_isIntegerType( xtype[i] ) )
                    IQD_correctBoundDirToInt(intTol, lbdir[i], ubdir[i]);
            }
            
            neyedirs = n;
        }
        
        
        for(int i = neyedirs; i < ndirs; i++)
        {
            const double *diri = dirs[i];
            int index, nnz = 0;
            
            
            for(int j = 0; j < n; j++)
            {
                if( IQD_abs( diri[j] ) >= zero_tol )
                {
                    nnz++;
                    if(nnz > 1)
                        break;
                    
                    index = j;
                }
            }
            
            
            if( nnz == 1 && MIP_isIntegerType( xtype[index] ) )
                IQD_correctBoundDirToInt(intTol, lbdir[i], ubdir[i]);
            
        }
        
        
    }
    
    
    #if IQD_DEBUG_MODE
        if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
            std::cout << IQD_PREPRINT << "Bounds calculated" << std::endl;
    #endif
    
    
    
    if( nonConvexConstrs || splitWithouIntVars || splitWay == IQD_SQ_MIX_SCHUR_DIAG_SDP )
    {
        vecsWeigth = (double *) malloc( ndirs * sizeof(double) );
        
        if( !vecsWeigth )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        #pragma ivdep
        #pragma GCC ivdep
        for(int i = 0; i < ndirs; i++)
            vecsWeigth[i] = (ubdir[i]*ubdir[i] + lbdir[i]*lbdir[i])/4.0 - ubdir[i]*lbdir[i]/2.0;
    }
    
    
    
    //setting quadratic part in objective reformulated problem
    if( nlamobj == 0 && !certainOfNonconvexQ )
    {
        IQD_setAllarray( ndirs, lambdaObj, 0.0 );
        
        ret = rprob.setObjQuadCoefsMatrix( Q );
        
        if( ret != 0 )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTERRORNUMBER(ret);
            #endif
            
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        
        if( oprob->objFactor != 1.0 )
            rprob.Q.multiplyAllElements( oprob->objFactor );
        
        
    }
    else
    {
        
        if( splitWay == IQD_SQ_SCHUR )
        {
            //IQD_setAllarray( n, lamObj, 0.0 );
            //calculateRbyLambdaAndv( n, ndirs, dirs, lamObj, M, true );
            
            
            int start = useEyeVector ? n : 0;
            for(int i = 0; i < n; i++)
            {
                if( lambda[i] < -zero_tol )
                {
                    lambdaObj[start] = lambda[i];
                    start++;
                }
            }
            
            #if IQD_DEBUG_MODE
                if( nonConvexConstrs )
                    assert( start - n == nlamobj );
            #endif
            
            
            
            //calculating P
            IQD_setAllarray(n*n, M, 0.0);
            
            for(int k = 0; k < n; k++)
            {
                if( lambda[k] > zero_tol ) //we are calculating P
                {
                    const double *vk = schvecs[k];
                    
                    for( int i = 0; i < n; i++ )
                    {
                        double *mi = &M[i*n];
                        
                        const double lvki = lambda[k] * vk[i];
                        
                        for(int j = 0; j <= i; j++ )
                            mi[j] += lvki * vk[j];
                    }
                }
            }
            
            
            //adding a small tolerance on diagonal of P.
            if( rel_eps_diag_P > 0.0 || abs_eps_diag_P )
                IQD_addEpsOnNonZeroDiag(n, M, zero_tol,  abs_eps_diag_P, rel_eps_diag_P);
            
            
            ret = rprob.setObjQuadCoefsMatrix(M, zero_tol, n, n);
            
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
        }
        else
        {
            ret = splitMatrix( splitWay, Q, oprob->objFactor, ndirs, dirs, vecsWeigth, sdpSolver, sdpSolverParams, lambdaObj, M);
            
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                
                code = ret;
                goto termination;
            }
            
            //M has the matrix R
            
            for(int i = 0; i < n; i++ )
            {
                double * Mrow = &M[i*n];
                #pragma ivdep
                #pragma GCC ivdep
                for(int j = 0; j <= i; j++)
                    Mrow[j] = -Mrow[j];
            }
            
            //M has the matrix -R
            //P = Q - R
            
            Q.accumulatedSumToMatrix(M, 1.0, false); //Q.accumulatedSumToMatrix(1.0, M, false);
            
            if( rel_eps_diag_P > 0.0 || abs_eps_diag_P > 0.0 )
                IQD_addEpsOnNonZeroDiag(n, M, zero_tol, abs_eps_diag_P, rel_eps_diag_P);
            
            
            ret = rprob.setObjQuadCoefsMatrix(M, zero_tol, n, n);
            
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            
            /*//M has the matrix R
            //P = Q - R
            //-P = -Q + R
                    
            Q.accumulatedSumToMatrix(-1.0, M);
            
            Q.printAllMatrix();
            
            
            
            //adding a smal tolerance on diagonals... Note, M has -P
            if( rel_eps_diag_P > 0.0 )
                IQD_addRelativeEpsOnNonZeroDiag(n, M, zero_tol, rel_eps_diag_P);
            
            ret = rprob.setObjQuadCoefsMatrix(M, zero_tol, n, n);
            
            if( ret != 0 )
            {
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            //Note we set -P matrix. So, we have to change the signals...
            
            rprob.Q.multiplyAllElements(-1.0); */
            
        }
        
    
    }
    
    
    if( max_cpu_time < INFINITY )
    {
        if( (double(clock()-clockStart) )/CLOCKS_PER_SEC > max_cpu_time  )
        {
            code = IQD_MAX_TIME_STOP;
            goto termination;
        }
    }
    
    if( max_time < INFINITY )
    {
        if( IQD_getTime() - timeStart > max_time )
        {
            code = IQD_MAX_TIME_STOP;
            goto termination;
        }
    }
    
    
    
    if( !IQD_isDiagonalSplitQWay(splitWay) )
        splitWay = IQD_SQ_MIX_SCHUR_DIAG_SDP;
    
    
    
    //spliting without integer vars on objective
    if( splitWithouIntVars && nlamobj > 0)
    {
        const int *xtype = oprob->xtype;
        
        
        lambdaObjNoInt = (double *) malloc( ndirs * sizeof(double) );
        
        Qnoint = new (nothrow) IQD_SparseMatrix;
        
        if( !lambdaObjNoInt || !Qnoint ) 
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTMEMERROR;
            #endif
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        ret = Qnoint->initialize(n,n, true);
        
        if( ret != 0 )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTERRORNUMBER(ret);
            #endif
            code = IQD_MEMORY_ERROR;
            goto termination;
        }
        
        
        for(int i = 0; i < n; i++)
        {
            if( !minlpproblem::MIP_isIntegerType( xtype[i] )  )
            {
                ret = Qnoint->setRowStructure(i, i, Q, true); //ret = Qnoint->setRowStructureAndValues( i, Q[i] );
                if( ret != 0 )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTERRORNUMBER(ret);
                    #endif
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
                
                //now, we will zerate postions related to integer vars. Note as Q store lower triangle, integer positions only appear in rows after row intVars[0] (columns indices must be slower than row indices)
                
                if( i > intVars[0] )
                {
                    //IQD_SparseRow &row = (*Qnoint)[i];
                    
                    const int *rcols = Qnoint->operator[](i);
                    double *rvals = Qnoint->operator()(i);
                    const unsigned int nel = Qnoint->getNumberOfElementsAtRow(i); //row.getNumberOfElements();
                    
                    for(unsigned int j = 0; j < nel; j++)
                    {
                        if( minlpproblem::MIP_isIntegerType(xtype[rcols[j]]) ) //( minlpproblem::MIP_isIntegerType( xtype[ row[j].getColumn() ] ) )
                        {
                            rvals[j] = 0.0; //row[j].setValue( 0.0 ); //removing integer variables...
                        }
                    }
                }
            }
        }
        
        
        /* std::cout << "Q: " << std::endl;
        Q.printSparseMatrix();
        
        std::cout << "Qnoint: " << std::endl;
        Qnoint->printAllMatrix(); */
        
        
        
        
        if( Qnoint->getNumberOfElements() > 0 )
        {
            
            if( oprob->objFactor != 1.0 )
                Qnoint->multiplyAllElements( oprob->objFactor );
            
            
            ret = splitMatrix( splitWay, *Qnoint, 1.0, ndirs, dirs, vecsWeigth, sdpSolver, sdpSolverParams, lambdaObjNoInt, M);
            
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    std::cerr << IQD_PREPRINT << "Error " << ret << " on splitMatrix" << IQD_GETFILELINE << endl;
                #endif
                
                code = ret;
                goto termination;
            }
            
            std::cout << "Rnoint: " << std::endl;
            IQD_printMatrix( n, n, M );
            
            
            //M has the matrix R
            //P = Q - R
            //-P = -Q + R
            
            Qnoint->accumulatedSumToMatrix(M, -1.0, false);
            
            //Now, M has -P
            
            //adding a small tolerance on diagonal of P. Note, M has -P
            if( rel_eps_diag_P > 0.0 || abs_eps_diag_P > 0.0 )
                IQD_addEpsOnNonZeroDiag(n, M, zero_tol, -abs_eps_diag_P, rel_eps_diag_P);
            
            
            delete Qnoint;
            PnoInt = new (std::nothrow) IQD_SparseMatrix;
            
            if( !PnoInt )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTMEMERROR;
                #endif
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            PnoInt->setSymmetricFlag(true);
            
            ret = PnoInt->setStructureAndValues(M, zero_tol, n, n);
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            PnoInt->multiplyAllElements(-1.0);
            
            std::cout << "PnoInt: " << std::endl;
            PnoInt->printAllMatrix();
            //PnoInt-> printSparseMatrix();
        }
        else
        {
            IQD_setAllarray(ndirs, lambdaObjNoInt, 0.0);
        }
        
        
        //IQD_getchar();
    }
    
    
    
    //taking care of constraints...
    
    
    for( int i = 0; i < m; i++ )
    {
        IQD_SparseMatrix &QC = oprob->QC[i];
        bool nlConstr;
        
        #if IQD_DEBUG_MODE
            std::cout << IQD_PREPRINT "Checking QC in constraint "<< i << "\n";
        #endif
        
        ret = oprob->getConstraintsNonLinearTermFlag(i, nlConstr);
        //ret += oprob->getConstraintBounds(i, lc[i], uc[i]);
        
        #if IQD_DEBUG_MODE
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                
                code = IQD_UNDEFINED_ERROR;
                goto termination;
            }
        #endif
        
        
        if( nlConstr )
        {
            rprob.setConstraintNonLinearTermFlag( i, true );
        }
        
        
        if( QC.getNumberOfElements() == 0 )
        {
            ret = rprob.setConstraintLinearPart(i, A.getNumberOfElementsAtRow(i),  A[i], A(i)); //ret = rprob.setConstraintLinearPart(i, A[i]);
            ret += rprob.setConstraintLowerBound(i, lc[i]);
            ret += rprob.setConstraintUpperBound(i, uc[i]);
            
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            continue;
        }
        
        
        if(lc[i] <= -MIP_INFINITY && uc[i] >= MIP_INFINITY)
            continue; //free constraint
        
        
        
        if( uc[i] < MIP_INFINITY )
        {
            if( ubconvex[i] )
            {
                ret = rprob.setConstraintQuadCoefsMatrix(i, QC);
            }
            else
            {
                #if IQD_DEBUG_MODE
                    if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
                        std::cout << IQD_PREPRINT << "Split quadratic constraint " << i << " on upper bound" << std::endl;
                #endif
                
                lambdaConstr[i] = (double *) malloc( ndirs * sizeof(double) );
                
                if( !lambdaConstr[i] )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTMEMERROR;
                    #endif
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
                
                //cout << "Aloquei lambdaConstr em ub" << endl;
                
                #if IQD_DEBUG_MODE
                    assert( ndirs == n + nlamobj );
                #endif
                
                ret = splitMatrix( splitWay, QC, 1.0, ndirs, dirs, vecsWeigth, sdpSolver, sdpSolverParams, lambdaConstr[i], M);
                
                if( ret != 0 )
                {
                    #if IQD_DEBUG_MODE
                        std::cerr << IQD_PREPRINT << "Error " << ret << " on splitMatrix" << IQD_GETFILELINE << std::endl;
                        
                        if( n < IQD_N_TO_PRINT_MATRIX)
                            QC.printSparseMatrix();
                    #endif
                    
                    code = ret;
                    goto termination;
                }
                
                //M has the matrix R
                //P = Q - R
                //-P = -Q + R
                
                QC.accumulatedSumToMatrix(M, -1.0, false);
                
                //Now, M has -P
                
                //adding a small tolerance on diagonal of P. Note, M has -P
                if( rel_eps_diag_P > 0.0 || abs_eps_diag_P > 0.0 )
                    IQD_addEpsOnNonZeroDiag(n, M, zero_tol, -abs_eps_diag_P, rel_eps_diag_P);
                
                
                
                ret = rprob.setConstraintQuadCoefsMatrix(i, M, zero_tol, n, n);
                
                
                //Note we set -P matrix. So, we have to change the signals...
                rprob.QC[i].multiplyAllElements(-1.0);
            }
            
            
            ret += rprob.setConstraintLinearPart(i, A.getNumberOfElementsAtRow(i), A[i], A(i)); //ret += rprob.setConstraintLinearPart(i, A[i]);
            ret += rprob.setConstraintUpperBound(i, uc[i]);
            
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    std::cerr << IQD_PREPRINT << "Error on minlpproblem" << IQD_GETFILELINE << endl;
                #endif
                
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
        }
        
        
        
        if( lc[i] > -MIP_INFINITY )
        {
            int cindex;
            
            
            if( uc[i] >= MIP_INFINITY )
            {//that is only a >= constraint...
                cindex = i; 
            }
            else
            {//that is a double bounded constraint... we had to add a new constraint
                
                cindex = m + cnewm;
                
                cnewm++;
                
                origConstr[i] = cindex;
                origConstr[cindex] = i;
                
                
                if( nlConstr )
                {
                    int r;
                    
                    rprob.setConstraintNonLinearTermFlag( cindex, nlConstr );
                    
                    r = rprob.setJacobianRowStructure( cindex, oprob->J.getNumberOfElementsAtRow(i), oprob->J[i] ); //r = rprob.setJacobianRowStructure( cindex, oprob->J[i] );
                    
                    if( r != 0 )
                    {
                        #if IQD_DEBUG_MODE
                            IQD_PRINTERRORNUMBER(ret);
                        #endif
                        
                        code = IQD_MEMORY_ERROR;
                        goto termination;
                    }
                }
                
            }
            
            
            ret = rprob.setConstraintLinearPart( cindex, A.getNumberOfElementsAtRow(i), A[i], A(i) ); //ret = rprob.setConstraintLinearPart( cindex, A[i] );
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            
            if( lbconcave[i] )
            {
                
                ret = rprob.setConstraintQuadCoefsMatrix( cindex, QC);
                if( ret != 0 )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTERRORNUMBER(ret);
                    #endif
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
                
                //we do not need multiply constraint by -1 because is convex...
                
                rprob.setConstraintLowerBound( cindex, lc[i]);
            }
            else
            {
                
                #if IQD_DEBUG_MODE
                    if( n >= IQD_N_TO_PRINT_MSG_ON_SPLIT_AND_CHECK_TIME )
                        std::cout << IQD_PREPRINT << "Split quadratic constraint " << i << " on lower bound" << std::endl;
                #endif
                
                
                //note we have to multiply this constraint by -1.0 in order to have <= constraint...
                
                lambdaConstr[cindex] = (double *) malloc( ndirs * sizeof(double) );
                
                if( !lambdaConstr[cindex] )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTMEMERROR;
                    #endif
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
                
                
                
                ret = splitMatrix( splitWay, QC, -1.0, ndirs, dirs, vecsWeigth, sdpSolver, sdpSolverParams, lambdaConstr[cindex], M);
                
                if( ret != 0 )
                {
                    #if IQD_DEBUG_MODE
                        std::cerr << IQD_PREPRINT << "Error " << ret << " on splitMatrix" << IQD_GETFILELINE << endl;
                    #endif
                    
                    code = ret;
                    goto termination;
                }
                
                
                if( n <= 20 )
                {
                    for(int w = 0; w < n + nlamobj; w++)
                        cout << "lambdaConstr[" << w << "]: " << lambdaConstr[cindex][w] << endl;
                }
                
                
                //cout << "R: " << endl;
                //IQD_printMatrix(n, n, M);
                
                
                
                //M has the matrix R
                //P = Q - R
                //-P = -Q + R
                
                //note we are summing -(-Q) + R, so the factor belowe is 1 instead of -1
                QC.accumulatedSumToMatrix(M, 1.0, false); //QC.accumulatedSumToMatrix(1.0, M, false);
                
                #if IQD_DEBUG_MODE
                if( n <= IQD_N_TO_PRINT_MATRIX )
                {
                    std::cout << "QC[" << i << "]: " << std::endl;
                    //QC.printAllMatrix();
                    QC.printSparseMatrix();
                }
                #endif
                
                
                //Now, M has -P
                
                //adding a small tolerance on diagonal of P. Note, M has -P
                if( rel_eps_diag_P > 0.0 || abs_eps_diag_P > 0.0 )
                    IQD_addEpsOnNonZeroDiag(n, M, zero_tol, -abs_eps_diag_P, rel_eps_diag_P);
                
                
                ret = rprob.setConstraintQuadCoefsMatrix( cindex, M, zero_tol, n, n);
                
                if( ret != 0 )
                {
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
                
                //Note we set -P matrix. So, we have to change the signals...
                rprob.QC[cindex].multiplyAllElements(-1.0);
                
                rprob.A.multiplyAllElementsAtRow(cindex, -1.0); //rprob.A[cindex].multiplyAllElements(-1.0);
                
                
                ret = rprob.setConstraintUpperBound( cindex, -lc[i]);
                
                if( ret != 0 )
                {
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
                
                chgSignal[cindex] = true;
                
                
                if( n <= 20 )
                {
                    cout << "P: " << endl;
                    rprob.QC[cindex].printAllMatrix(); IQD_printMatrix(n,n, M);
                }
                
            }
            
        } //end of if( lb > -MIP_INFINITY )
        
        
        if( max_cpu_time < INFINITY )
        {
            if( (double(clock()-clockStart) )/CLOCKS_PER_SEC > max_cpu_time  )
            {
                code = IQD_MAX_TIME_STOP;
                goto termination;
            }
        }
        
        if( max_time < INFINITY )
        {
            if( IQD_getTime() - timeStart > max_time )
            {
                code = IQD_MAX_TIME_STOP;
                goto termination;
            }
        }
        
    }
    
    mdis = cnewm;
    
    #if IQD_DEBUG_MODE
        assert( newm == cnewm );
    #endif
    
    
    
    
    //now, we check by vecs used to decompose matrices...
    
    nnewdirs = 0;
    
    for(int i = 0; i < ndirs; i++)
    {
        
        bool lused = false;
        
        if( lambdaObj[i] >= -zero_tol )
        {
            
            if( lambdaObjNoInt )
            {
                if( lambdaObjNoInt[i] < -zero_tol )
                    lused = true;
            }
            
            if( !lused )
            {
                for( int j = 0; j < m + newm; j++ )
                {
                    if( lambdaConstr[j] )
                    {
                        if( lambdaConstr[j][i] < -zero_tol )
                        {
                            lused = true;
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            lused = true;
        }
        
        
        if(lused)
        {
            //IQD_copyArray( n, dirs[i], this->v[nnewdirs] );
            
            v.setRowStructureAndValues( nnewdirs,  dirs[i], n, zero_tol);
            
            
            lambdaObj[nnewdirs] = lambdaObj[i];
            if( lambdaObjNoInt )
                lambdaObjNoInt[nnewdirs] = lambdaObjNoInt[i];
            
            for(int j = 0; j < m + newm; j++)
            {
                if( lambdaConstr[j] )
                    lambdaConstr[j][nnewdirs] = lambdaConstr[j][i]; //that only works because nnewdirs <= i
            }
            
            ly[nnewdirs] = lbdir[i];
            uy[nnewdirs] = ubdir[i];
            
            nnewdirs++;
        }
        
    }
    
    /*std::cout << "v: " << std::endl;
    v.printSparseMatrix();
    IQD_getchar(); */
    
    
    #if IQD_DEBUG_MODE
        assert(nnewdirs <= ndirs);
    #endif
    
    if( nnewdirs != ndirs )
    {
        for(int i = 0; i < m + newm; i++)
        {
            if( lambdaConstr[i] )
            {
                double *auxp = (double *) realloc( lambdaConstr[i], nnewdirs * sizeof(double) );
                
                if( auxp )
                    lambdaConstr[i] = auxp;
            }
        }
    }
    
    dimv = nnewdirs;
    
    
    
    
    
    
    
    //we should add nnewdirs auxiliary constarints. But, directions having only one nonzero value can be used to set variable bounds instead of auxiliary constraints, specially the eye columns. So, we calculate the real number of auxiliary constraints we addConstraints
    
    
    
    for(int i = 0; i < nnewdirs; i++)
    {
        if( v.getNumberOfElementsAtRow(i) > 1 ) //( v[i].getNumberOfElements() > 1 )
            nauxconstrs++;
    }
    
    
    
    
    cols = (int *) malloc( (n + nnewdirs)* sizeof(int) );
    values = (double *) malloc( (n + nnewdirs)* sizeof(double) );
    
    ret = rprob.addVariables( nnewdirs ); //w variables
    ret+= rprob.addConstraints( nauxconstrs ); //we could add the auyxiliary constarints after... we only do it here to try make things more organized... Anyway, secant constraint will be added directly in solver object in B&B procedure
    
    if( !cols || !values || ret != 0 )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    //auxiliary constraints...
    indvconstrs = m + newm;
    k = 0;
    for( int i = 0; i < nnewdirs; i++ )
    {
        //IQD_SparseRow &vrow = v[i];
        
        const int nel = v.getNumberOfElementsAtRow(i); // vrow.getNumberOfElements();
        
        if( nel == 1 )
        {
            const double d =  v(i)[0]; // vrow[0].getValue();
            const int ind = v[i][0]; //vrow[0].getColumn();
            
            double lb = ly[i];
            double ub = uy[i];
            
            
            if( d != 1.0 )
            {
                #if IQD_DEBUG_MODE
                    assert( d != 0.0 );
                #endif
                
                lb = lb/d;
                ub = ub/d;
                
                if( d < 0 )
                    IQD_swap(lb, ub);
            }
            
            
            ret = rprob.setVariableLowerBound(ind, lb);
            ret += rprob.setVariableUpperBound(ind, ub);
            
            #if IQD_DEBUG_MODE
                assert(ret == 0);
            #endif
            
        }
        else if( nel > 1 )
        {
            //ret = rprob.setConstraintLinearPart( i + aux, v[i], n, zero_tol );
            
            //ret = rprob.setConstraintLinearPart( indvconstrs + k, vrow );
            ret = rprob.setConstraintLinearPart( indvconstrs + k, nel, v[i], v(i) );
            
            if( ret != 0 )
            {
                #if IQD_MEMORY_ERROR
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                code = OPT_MEMORY_ERROR;
                goto termination;
            }
            
            
            //we could let those constraints free because the bounds will be updated by Branch-And-Bound
            rprob.setConstraintLowerBound( indvconstrs + k, ly[i] );
            rprob.setConstraintUpperBound( indvconstrs + k, uy[i] );
            
            k++;
        }
    }
    
    
    #if IQD_DEBUG_MODE
        assert( k == nauxconstrs );
    #endif
    
    
    for( int j = 0; j < nnewdirs; j++ )
    {
        if( lambdaObj[j] < -zero_tol )
        {
            ret = rprob.setObjLinearCoefficient( n+j, -0.5*lambdaObj[j] );
            
            #if IQD_DEBUG_MODE
                assert(ret == 0);
            #endif
        }
    }
    
    if( nonConvexConstrs )
    {
        for( int i = 0; i < m; i++ )
        {
            
            std::cout << "lbconv["<<i<<"]: " << lbconcave[i] << " ubconv["<<i<<"]: " << ubconvex[i] << std::endl;
            
            
            
            
            if( !lbconcave[i] || !ubconvex[i] )
            {
                //if we have a double bounded constraint and ub part is convex, we only build the relaxation for the new constraint added from disjunction (in the next for after this)...
                if( ubconvex[i] &&  ( lc[i] > -MIP_INFINITY && uc[i] < MIP_INFINITY)   )
                    continue;
                
                
                int nzs;
                double *lambdaci = lambdaConstr[i];
                
                
                rprob.getConstraintLinearPart(i, &nzs, cols, values);
                
                
                for(int j = 0; j < nnewdirs; j++)
                {
                    if(lambdaci[j] < -zero_tol)
                    {
                        cols[nzs] = n + j;
                        values[nzs] = -0.5*lambdaci[j];
                        nzs++;
                    }
                }
                
                ret = rprob.setConstraintLinearPart(i, nzs, cols, values);
                
                if( ret != 0 )
                {
                    #if IQD_MEMORY_ERROR
                        IQD_PRINTERRORNUMBER(ret);
                    #endif
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
                
            }
        }
        
        
        //now, we have to add lambda to new constraints generated by quadratic double bound constraints...
        
        aux = m + newm;
        for( int i = m; i < aux; i++)
        {
            //those new constraints are >= part of double bound constraint. This part can be convex (concave).
            
            if( !lbconcave[ origConstr[i] ]  )
            {
                int nzs;
                double *lambdaci = lambdaConstr[i];
                
                ret = rprob.getConstraintLinearPart(i, &nzs, cols, values);
                
                #if IQD_DEBUG_MODE
                    assert(ret == 0);
                #endif
                
                for(int j = 0; j < nnewdirs; j++)
                {
                    if( lambdaci[j] < -zero_tol )
                    {
                        cols[nzs] = n + j;
                        values[nzs] = -0.5*lambdaci[j];
                        nzs++;
                    }
                }
                
                ret = rprob.setConstraintLinearPart( i, nzs, cols, values );
                
                if( ret != 0 )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTERRORNUMBER(ret);
                    #endif
                    code = IQD_MEMORY_ERROR;
                    goto termination;
                }
            }
        }
    
    }
    
    /*cout << "dirs: \n";
    IQD_printMatrix(dimv, n, &dirs[0][0]);
    
    cout << "lambdaObj: \n";
    for(int i = 0; i < dimv; i++)
        cout << lambdaObj[i] << " ";
    
    cout << "\n";
    
    for( int j = 0; j < m + newm; j++ )
    {
        if( lambdaConstr[j] )
        {
            cout << j << ": ";
            
            for(int w = 0; w < dimv; w++)
                cout << lambdaConstr[j][w] << " ";
            
            cout << "\n";
        }
    }*/
    
    
    //rprob.print();
    
    //rprob.writeMIQCPModelInAMPLFile("rprob.mod");
    //IQD_getchar();
    
    
    code = 0;
    
termination:
    
    
    
    if(lbconcave)		free(lbconcave);
    if(M)			free(M);
    if(lambda)		free(lambda);
    if( schvecs )
    {
        if(schvecs[0])	free(schvecs[0]);
        free(schvecs);
    }
    
    if(eye)			free(eye);
    if(dirs)		free(dirs);
    if(lbdir)		free(lbdir);
    if(vecsWeigth)	free(vecsWeigth);
    
    if(intVars)		free(intVars);
    
    //if(lamObj)		free(lamObj);
    //if(lamConstr)	free(lamConstr);
    
    if(cols)		free(cols);
    if(values)		free(values);
    
    
    //rprob.writeMIQCPModelInAMPLFile( "rprob.mod" );
    
    
    return code;
}










unsigned int IQD_RefProb2::getNDual()
{
    return 2*(oprob->n + oprob->m + 2*dimv);
}













