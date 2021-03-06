

#include <cstdlib>
#include <climits>
#include <iostream>

#include <thread>

#include "OPT_solvers.hpp"
#include "OPT_tools.hpp"


#define OPT_CHECKMOSEKLICENSEERROR(r) r >= MSK_RES_ERR_LICENSE && r <= MSK_RES_ERR_LICENSE_NO_SERVER_LINE




using namespace optsolvers;
using namespace std;



OPT_Mosek::OPT_Mosek():OPT_NLPSolver()
{
    initialize();
}



OPT_Mosek::~OPT_Mosek()
{
    //desallocateMemory();
    deallocateSolverEnv();
}



int OPT_Mosek::allocateConstrStructures(const int m)
{
    


    if( m > maux )
    {
        /*bool *auxb = (bool *) realloc( nlConstr, m * sizeof(bool) );

        if(!auxb)
            return OPT_MEMORY_ERROR;

        nlConstr = auxb;
        
        for(int i = maux; i < m; i++)
            nlConstr[i] = false; */
        
        int r = OPT_realloc(nlConstr, m);
        OPT_IFERRORRETURN(r, r);
        
        if(m > maux)
            OPT_setAllArray(m-maux, &nlConstr[maux], false);
    }

    return OPT_NLPSolver::allocateConstrStructures(m);
}



void OPT_Mosek::updatehasNLConstr()
{
    int m = 0;
    
    getNumberOfNLConstraints(m);
    
    hasNLConstrs = false;
    
    for(int i = 0; i < m; i++)
    {
        if( nlConstr[i] )
        {
            hasNLConstrs = true;
            break;
        }
    }
}




void OPT_Mosek::deallocateSolverEnv()
{
#if OPT_HAVE_MOSEK
    if(task)
    {
        MSK_deletetask(&task);
        task = NULL; //maybe we dot need do it, but...
    }
    
    if(env)
    {
        MSK_deleteenv(&env);
        env = NULL; //maybe we dot need do it, but...
    }
#endif
    
    desallocateMemory();
    
    OPT_secFree(nlConstr);
}

#if 0
void OPT_Mosek::getDualSolution( double *dualConstrs, double *dualVarBounds, const bool correctSignal )
{
    
    if( dualConstrs )
    {
        int m = 0;
        MSKrescodee r = MSK_getnumcon(task, &m);
        
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
        
        OPT_copyArray(m, dualSolC, dualConstrs);
        
        // am not sure if we should correct this signal. Maybe we should correct signal of other solvers like ipopt and knitro because mosek convention is the same of linear programming solvers...
        if( correctSignal )
        {
            for(int i = 0; i < m; i++)
                dualConstrs[i] = -dualSolC[i];
        }
    }
    
    
    if( dualVarBounds )
    {
        int n;
        MSKrescodee r = MSK_getnumvar(task, &n);
    
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
        
        OPT_copyArray(2*n, dualSolV, dualVarBounds);
    }
    
}
#endif


int OPT_Mosek::getNumberOfIterations(long unsigned int &niter)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    int value;
    
    r = MSK_getnaintinf(task, "MSK_IINF_INTPNT_ITER", &value);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    niter = value;
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



OPT_LISTSOLVERS OPT_Mosek::getSolverCode()
{
    return OPT_MOSEK;
}


int OPT_Mosek::getVariableType( const int index, OPT_VARTYPE &varType )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    MSKvariabletypee vt;
    
    r = MSK_getvartype(task, index, &vt);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    varType = vt == MSK_VAR_TYPE_INT ? OPT_VT_INTEGER : OPT_VT_CONTINUOUS;
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


void OPT_Mosek::initialize()
{
    OPT_QCPSolver::initialize();
    
    hasNLConstrs = nlObj = false;
    nlConstr = NULL;
    
    nlObjFactor = 1.0;
    
    #if OPT_HAVE_MOSEK
        origSolverRetCode = MSK_SOL_STA_UNKNOWN;
        env = NULL;
        task = NULL;
    #endif
}


int OPT_Mosek::initSolverEnv(const int maxConstrs, const int maxVars, const int maxQuadNz)
#if OPT_HAVE_MOSEK
{
    int code;
    MSKrescodee r;
    
    //cout << "Entrei em initSolverEnv" << endl;
    
    nlObjFactor = 1.0;
    
    /* Create the mosek environment. */
    r = MSK_makeenv(&env, NULL);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    #if MSK_VERSION_MAJOR < 8 
    {
        /* Initialize the environment. */
        r = MSK_initenv(env);
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    }
    #endif
    
    
    r = MSK_maketask(env, maxConstrs, maxVars, &task);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    if( maxQuadNz > 0 )
    {
        r = MSK_putmaxnumqnz( task, maxQuadNz ) ;
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    }
    
    
    
    r = MSK_putintparam(task, MSK_IPAR_CHECK_CONVEXITY, MSK_CHECK_CONVEXITY_NONE); //turn off convexity checker
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }
    
    
    
    #if OPT_DEBUG_MODE
        //we need do it to generate a lp file
        r = MSK_putintparam(task, MSK_IPAR_WRITE_GENERIC_NAMES, MSK_ON);
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
        }
    #endif
    
    
    code = 0;
    
    
    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek:: setMaxCPUTime(const double time)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putdouparam(task, MSK_DPAR_OPTIMIZER_MAX_TIME,  time);
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setNumberOfThreads(const int nthreads)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    r = MSK_putintparam(task, MSK_IPAR_NUM_THREADS, nthreads);
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    /*r = MSK_putintparam(task, MSK_IPAR_INTPNT_MULTI_THREAD,   nthreads == 1 ? MSK_OFF : MSK_ON );
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
    }*/
    
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setOutputLevel( const int level )
#if OPT_HAVE_MOSEK
{
    const int r = MSK_putintparam(task, MSK_IPAR_LOG, level);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setRelativeDualTol( const double tol )
#if OPT_HAVE_MOSEK
{
    const MSKrescodee r = MSK_putdouparam(task, MSK_DPAR_INTPNT_NL_TOL_DFEAS, tol);  //default: 1.0e-8
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setRelativeOptimalityTol( const double tol )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    r = MSK_putdouparam(task, MSK_DPAR_INTPNT_NL_TOL_REL_GAP, tol); //default: 1.0e-6
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    r = MSK_putdouparam(task, MSK_DPAR_INTPNT_NL_TOL_MU_RED, tol * 1e-4 ); //default: 1.0e-12
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setRelativePrimalTol( const double tol )
#if OPT_HAVE_MOSEK
{
    const MSKrescodee r = MSK_putdouparam(task, MSK_DPAR_INTPNT_NL_TOL_PFEAS, tol); //default: 1.0e-8
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setObjCutUpperBound(const double objUBound)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    r = MSK_putdouparam(task, MSK_DPAR_UPPER_OBJ_CUT, objUBound);
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setObjCutLowerBound(const double objLBound)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    r = MSK_putdouparam(task, MSK_DPAR_LOWER_OBJ_CUT, objLBound);
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setDoubleParameter(const char *param, const double value)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    r = MSK_putnadouparam(task, param, value);
    if ( r != MSK_RES_OK )
    {
        printDblParamErrorMsg(r, param, value );
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setIntegerParameter(const char *param, const int value )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    r = MSK_putnaintparam(task, param, value);
    if ( r != MSK_RES_OK )
    {
        printIntParamErrorMsg(r, param, value );
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setStringParameter(const char *param, const char *value)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    r = MSK_putnastrparam(task, param, value);
    if ( r != MSK_RES_OK )
    {
        printStrParamErrorMsg(r, param, value );
        
        return OPT_BAD_INPUT;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif





int OPT_Mosek::setVariableType( const int index, const OPT_VARTYPE varType )
#if OPT_HAVE_MOSEK
{
    const MSKrescodee r = MSK_putvartype(task, index, varType == OPT_VT_INTEGER ? MSK_VAR_TYPE_INT : MSK_VAR_TYPE_CONT );
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



#if OPT_HAVE_MOSEK
int OPT_Mosek::getSolution( const MSKsoltypee solType, const bool storeSol, const bool storeConstrs, const bool storeDualSol )
{
    int r, nvars, ncons;
    

    //feasSol = false;
    
    MSK_getnumvar(task, &nvars);
    MSK_getnumcon(task, &ncons);
    
    
    
    r = MSK_getprimalobj(task, solType, &objValue);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_UNDEFINED_ERROR;
        //goto termination;
    }
    
    
    if(storeSol)
    {
        r = MSK_getsolutionslice(task, solType, MSK_SOL_ITEM_XX, 0, nvars, sol);
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_UNDEFINED_ERROR;
            //goto termination;
        }
    }
    
    
    if(storeConstrs)
    {
        r = MSK_getsolutionslice(task, solType, MSK_SOL_ITEM_XC, 0, ncons, constr);
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_UNDEFINED_ERROR;
            //goto termination;
        }
    }
    
    
    //feasSol = true;
    
    
    
    
    
    
    
    if( solType != MSK_SOL_ITG && storeDualSol )
    {
        //r = MSK_getsolutionslice(task, solType, MSK_SOL_ITEM_Y, 0, ncons, dualSolC);
        /* if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_UNDEFINED_ERROR;
            goto termination;
        }*/
        
        
        r = MSK_getdualobj( task, solType, &dualObjValue );
        /*
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif 
            
            return OPT_UNDEFINED_ERROR;
            goto termination;
        } */
        
        r = MSK_getsolutionslice(task, solType, MSK_SOL_ITEM_Y, 0, ncons, dualSolC );
        
        #if OPT_DEBUG_MODE
            if( r != MSK_RES_OK )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif 
                
                //return OPT_UNDEFINED_ERROR;
                //goto termination;
            }
        #endif
        
        //we do not change more the signal. Instead of, user can decide if wanna signal changing by calling getDualSolution...
        //for(int i = 0; i < ncons; i++)
            //dualSolC[i] = -dualSolC[i]; //lagrangean definition in Mosek has a minus signal before multipliers. So, I think (I SAD THINK) we need multiply the lagrange multipliers by -1;
        
        
        r = MSK_getsolutionslice(task, solType, MSK_SOL_ITEM_SLX, 0, nvars, dualSolV);
        
        #if OPT_DEBUG_MODE
            if( r != MSK_RES_OK )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif 
                
                //return OPT_UNDEFINED_ERROR;
                //goto termination;
            }
        #endif
        
        
        r = MSK_getsolutionslice(task, solType, MSK_SOL_ITEM_SUX, 0, nvars, &dualSolV[nvars]);
        
        #if OPT_DEBUG_MODE
            if( r != MSK_RES_OK )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif 
                
                //return OPT_UNDEFINED_ERROR;
                //goto termination;
            }
        #endif
        
    }
    
    
    
    
    
    
    return 0;
}
#endif





#if OPT_HAVE_MOSEK

int MSKAPI optsolvers::OPT_mosekStruc(void    *nlhandle, MSKintt  *numgrdobjnz, MSKidxt  *grdobjsub, MSKidxt  i, int      *convali, MSKidxt  *grdconinz, MSKidxt  *grdconisub, MSKintt  yo, MSKintt  numycnz, MSKCONST MSKidxt  *ycsub, MSKintt  maxnumhesnz, MSKintt  *numhesnz, MSKidxt  *hessubi, MSKidxt  *hessubj)
{
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        OPT_createFileToThread( tid );
        std::ostream &thOut = *OPT_thsOut[tid];
        
        thOut << "\tEntering at optsolvers::OPT_mosekStruc" << std::endl;
    }
    #endif
    
    OPT_Mosek *myMosek = (OPT_Mosek *) nlhandle;
    const OPT_SparseMatrix &J = myMosek->J;
    const OPT_SparseMatrix &lagH = myMosek->lagH;
    int n, m; //number of variables, number of constraints
    
    //printf("Entrei em OPT_mosekStruc Thread: %d\n", myMosek->threadNumber);
    //fflush(stdout);
    
    myMosek->getNumberOfVars(n);
    myMosek->getNumberOfConstraints(m);
    
    
    if( numgrdobjnz )
    {
        //all original variables can appear in nonlinear objective function
        numgrdobjnz[0] = myMosek->nlObj ? n : 0;
        
        //std::cout << "numgrdobjnz[0]: " << numgrdobjnz[0] << "\n";
    }
    
    if( grdobjsub )
    {
        if( myMosek->nlObj )
        {
            #pragma ivdep
            #pragma GCC ivdep
            for(int j = 0; j < n; j++)
                grdobjsub[j] = j;
            
            //for(int j = 0; j < n; j++)
                //std::cout << "grdobjsub["<<j<<"]: " << grdobjsub[j] << "\n";
        }
    }
    
    if( convali )
    {
        //I think I should put true if the constraint i is nonlinear...
        convali[0] = (int) myMosek->nlConstr[i];
        
        //std::cout << "i: " << i << " convali: " << convali[0] << "\n";
    }
    
    if( grdconinz )
    {
        //grdconinz[0] = myMosek->J[i].getNumberOfElements();
        grdconinz[0] = J.getNumberOfElementsAtRow(i);
        //std::cout << "i: " << i << " grdconinz: " << grdconinz[0] << "\n";
    }
    
    if( grdconisub )
    {
        //myMosek->J[i].getStructure(grdconisub);
        J.getRowStructure(i, grdconisub, NULL);
        
        //for(int j = 0; j < J.getNumberOfElementsAtRow(i); j++)
            //std::cout << "grdconisub["<<j<<"]: " << grdconisub[j] << "\n";
    }
    
    
    if( numhesnz )
    {
        
        if( myMosek->hasNLConstrs == false && yo == 0.0 )
        {
            numhesnz[0] = 0;
        }
        else
        {
            //we do not have hessian of constraints and lagrangean separately
            numhesnz[0] = lagH.getNumberOfElements();
        }
        
        //std::cout << "numhesnz[0]: " << numhesnz[0] << "\n";
        
        if( maxnumhesnz < numhesnz[0] )
        {
            //OPT_getchar();
            //Not enough space have been allocated for storing the Hessian.
            return 0;
        }
    }
    
    if( hessubi && hessubj )
    {
        int j = lagH.getStructure(hessubi, hessubj);
        
        #if OPT_DEBUG_MODE
            assert(j == numhesnz[0]);
        #endif
            
        //for(int w = 0; w < j; w++)
            //std::cout << "w: " << w << " hessubi: " << hessubi[w] << " hessubj: " << hessubj[w] << "\n";
    }
    
    //lagH.printSparseMatrix();
    //OPT_getchar();
    
    //printf("Sai de OPT_mosekStruc Thread: %d\n", myMosek->threadNumber);
    //fflush(stdout);
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        OPT_createFileToThread( tid );
        std::ostream &thOut = *OPT_thsOut[tid];
        
        thOut << "\tLeaving at optsolvers::OPT_mosekStruc" << std::endl;
        
        OPT_closeFileToThread(tid);
    }
    #endif
    
    return 0;
}


int MSKAPI optsolvers::OPT_mosekEval(void      *nlhandle, MSKCONST double    *xx, double    yo, MSKCONST double    *yc, double    *objval, MSKintt   *numgrdobjnz, MSKidxt   *grdobjsub, double    *grdobjval, MSKintt   numi, MSKCONST MSKidxt   *subi, double    *conval, MSKCONST MSKintt   *grdconptrb, MSKCONST MSKintt   *grdconptre, MSKCONST MSKidxt   *grdconsub, double    *grdconval, double    *grdlag, MSKintt   maxnumhesnz, MSKintt   *numhesnz, MSKidxt   *hessubi, MSKidxt   *hessubj, double    *hesval)
{
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        OPT_createFileToThread( tid );
        std::ostream &thOut = *OPT_thsOut[tid];
        
        thOut << "\tEntering at optsolvers::OPT_mosekEval" << std::endl;
    }
    #endif
    
    OPT_Mosek *myMosek = (OPT_Mosek *) nlhandle;
    OPT_NonLinearEval *nlEval = myMosek->nlEval;
    const unsigned int thnumber = myMosek->threadNumber;
    const bool *nlConstr = myMosek->nlConstr;
    bool *auxCEval = myMosek->auxCEval;
    
    
    int n, m; //number of variables, number of constraints
    
    
    int code;
    bool new_x = true;
    
    
    double *auxConstr = myMosek->auxValues;
    
    //printf("Entrei em OPT_mosekEval Thread: %d\n", myMosek->threadNumber);
    //std::cout << "Id da Thread: " << std::this_thread::get_id() << "\n";
    //fflush(stdout);
    
    myMosek->getNumberOfVars(n);
    myMosek->getNumberOfConstraints(m);
    
    
    //for(int w = 0; w < n; w++)
        //std::cout << "x["<<w<<"]: " << xx[w] << "\n";
    
    if( objval )
    {
        //we eval just nonlinear part
        
        if( myMosek->nlObj )
        {
            int r = nlEval->eval_nl_obj_part( thnumber, n, new_x, xx, objval[0] );
            
            if( r != 0 )
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    OPT_PRINTCALLBACKERRORNUMBER(r);
                #endif
                
                code = r;
                goto termination;
            }
            
            objval[0] *= myMosek->nlObjFactor;
            new_x = false;
        }
        else
        {
            objval[0] = 0.0;
        }
        
        //std::cout << "objval: " << objval[0] << "\n";
    }
    
    
    if( numgrdobjnz )
    {
        #if OPT_DEBUG_MODE 
        //just check if some day we can omit line below (no we cant)
            //assert( numgrdobjnz[0] == (myMosek->nlObj ? n : 0) );
        #endif
        numgrdobjnz[0] = myMosek->nlObj ? n : 0;
        
        //std::cout << "numgrdobjnz: " << numgrdobjnz[0] << "\n";
    }
    
    
    if( grdobjsub && myMosek->nlObj)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(int j = 0; j < n; j++)
        {
            #if OPT_DEBUG_MODE
                //assert( grdobjsub[j] == j );
            #endif
            grdobjsub[j] = j; //that is ridiculous, but we really need do it..
        }
        
        //for(int j = 0; j < n; j++)
            //std::cout << "grdobjsub["<<j<<"]: " << grdobjsub[j] << "\n";
    }
    
    
    if( grdobjval )
    {
        if( myMosek->nlObj )
        {
            const double of = myMosek->nlObjFactor;
            
            int r = nlEval->eval_grad_nl_obj_part( thnumber, n, new_x, xx, grdobjval );
            
            if(r != 0)
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    OPT_PRINTCALLBACKERRORNUMBER(r);
                #endif
                
                code = r;
                goto termination;
            }
            
            new_x = false;
            
            if( of != 1.0 )
            {
                OPT_multiplyAllArray(n, of, grdobjval);
                //for(int j = 0; j < n; j++)
                    //grdobjval[j] *= of;
            }
            
            //for(int j = 0; j < n; j++)
                //std::cout << "grdobjval["<<j<<"]: " << grdobjval[j] << "\n";
        }
    }
    
    
    if( (conval || grdconval || grdlag) && myMosek->hasNLConstrs )
    {
        for(decltype(numi) j = 0; j < numi; j++)
        {
            auto k = subi[j];
            
            if( nlConstr[k] )
                auxCEval[k] = true; //we assume mosek would not ask to evaluate a linear or quadratic consraint... (but it seems this is not true...)
        }
    }
    
    
    if(conval)
    {
        //I think that is the value of the constraints...
        
        if( myMosek->hasNLConstrs )
        {
            int r = nlEval->eval_nl_constrs_part( thnumber, n, m, new_x, auxCEval, xx, auxConstr );
            
            if(r != 0)
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    OPT_PRINTCALLBACKERRORNUMBER(r);
                #endif
                
                code = r;
                goto termination;
            }
            
            new_x = false;
        }
        
        for(decltype(numi) j = 0; j < numi; j++)
        {
            auto k = subi[j];
            
            if( nlConstr[k] )
                conval[j] = auxConstr[k];
            else
                conval[j] = 0.0; //mosek should not ask to evaluate a constraint having   no nonlinear terms, but it does...
        }
        
        //for(int j = 0; j < numi; j++)
            //std::cout << "conval[" << j << "]: " << conval[j] << "\n";
    }
    
    
    
    if( grdconval )
    {
        if( myMosek->hasNLConstrs )
        {
            OPT_SparseMatrix &J = myMosek->J;
            
            int r = nlEval->eval_grad_nl_constrs_part( thnumber, n, m, J.getNumberOfElements(), new_x, auxCEval, xx, J);
            
            if(r != 0)
            {
                #if OPT_PRINT_CALLBACK_ERROR_MSG
                    OPT_PRINTCALLBACKERRORNUMBER(r);
                #endif
                
                code = r;
                goto termination;
            }
            
            new_x = false;
            
            //TODO: Check if it is possible changing valus pointer like in hessian eval. If grdconval and grdconptrb are always the same, we can perform this change to avoid this copy
            
            for(decltype(numi) j = 0; j < numi; j++)
                J.getRowValues( subi[j],  &grdconval[ grdconptrb[j] ], NULL );
            
            /*for(int j = 0; j < numi; j++)
            {
                const int nz = myMosek->J.getNumberOfElementsAtRow( subi[j] );
                
                for(int k = 0; k < nz; k++)
                    std::cout << "subi["<<j<<"]: " <<  subi[j] << " grdconval["<<k<<"]: " << grdconval[ grdconptrb[j] + k ] << "\n";
            }*/
        }
        
    }
    
    if( grdlag )
    {
        //Compute and store the gradient of the Lagrangian. Note that it is stored as a dense vector.
        
        if( myMosek->nlObj && yo != 0.0 )
        {
            if( grdobjval )
            { //if grdobjval is not NULL, we already evaluate the gradient...
                
                if( yo == 1.0 )
                    OPT_copyArray(n, grdobjval, grdlag);
                else
                {
                    //for(int j = 0; j < n; j++)
                        //grdlag[j] = yo*grdobjval[j];
                    OPT_copyArrayTimesFactor(n, yo, grdobjval, grdlag);
                }
            }
            else
            {
                const double f = yo * myMosek->nlObjFactor;
                
                int r = nlEval->eval_grad_nl_obj_part(thnumber, n, new_x, xx, grdlag);
                
                if(r != 0)
                {
                    #if OPT_PRINT_CALLBACK_ERROR_MSG
                        OPT_PRINTCALLBACKERRORNUMBER(r);
                    #endif
                    
                    code = r;
                    goto termination;
                }
                
                new_x = false;
                
                if( f != 1.0 )
                    OPT_multiplyAllArray(n, f, grdlag);
            }
        }
        else
        {
            OPT_setAllArray(n, grdlag, 0.0);
        }
        
        
        if( myMosek->hasNLConstrs )
        {
            double factor;
            OPT_SparseMatrix &J = myMosek->J;
            
            if( grdconval == NULL )
            {//if gradconval != NULL, we have already set auxCEval and evaluate;
                
                int r = nlEval->eval_grad_nl_constrs_part(thnumber, n, m, J.getNumberOfElements(), new_x, auxCEval, xx, J);
                
                if(r != 0)
                {
                    #if OPT_PRINT_CALLBACK_ERROR_MSG
                        OPT_PRINTCALLBACKERRORNUMBER(r);
                    #endif
                    
                    code = r;
                    goto termination;
                }
                
                new_x = false;
            }
            
            for( decltype(numi) j = 0; j < numi; j++ )
            {
                const auto k = subi[j];
                
                factor = yc[ k ];
                if( factor != 0.0 )
                {
                    J.accumulateRowInArray(k, grdlag, -factor);
                }
            }
        }
    }
    
    
    
    if( numhesnz )
    {
        if( !myMosek->hasNLConstrs && yo == 0.0)
        {
            numhesnz[0] = 0;
            hesval = NULL; //we do not do any evaluation of hessian, so w eput NULL to avoid if below...
        }
        else
            numhesnz[0] = myMosek->lagH.getNumberOfElements();
    }
    
    
    if( hesval )
    {
        OPT_SparseMatrix &lagH = myMosek->lagH;
        double *lambda = auxConstr;
        
        
        OPT_setAllArray(m, lambda, 0.0);
        
        for(decltype(numi) j = 0; j < numi; j++)
        {
            auto k = subi[j];
            
            if( nlConstr[k] )
                lambda[k] = -yc[k]; ////in the lagrangian definition, we have a minus signal before the sum
        }
        
        
        if( myMosek->nlObj == false )
            yo = 0.0;
        
        //we try avoid copy values from lagH changing the values pointer directly to hesVal
        double *oldValues = lagH.baseStructure.values;
        lagH.baseStructure.values = hesval;
        
        int r = nlEval->eval_hessian_nl_lagran_part( thnumber, n, m, lagH.getNumberOfElements(), new_x, xx, yo*myMosek->nlObjFactor, lambda, lagH );
        
        lagH.baseStructure.values = oldValues; //restoring the original pointer
        
        if(r != 0)
        {
            #if OPT_PRINT_CALLBACK_ERROR_MSG
                OPT_PRINTCALLBACKERRORNUMBER(r);
            #endif
            
            code = r;
            goto termination;
        }
        
        lagH.getStructure(hessubi, hessubj);
        
        //myMosek->lagH.getStructureAndValues(hessubi, hessubj, hesval);
    }
    
    
    
    code = 0;
    
    
termination:
    
    if( (conval || grdconval || grdlag) && myMosek->hasNLConstrs )
    {
        for( decltype(numi) j = 0; j < numi; j++ )
        {
            auto k = subi[j];
            
            if( nlConstr[k] )
                auxCEval[ k ] = false;
        }
    }
    
    //OPT_getchar();
    
    //printf("Sai de OPT_mosekEval Thread: %d\n", myMosek->threadNumber);
    //fflush(stdout);
    
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        OPT_createFileToThread( tid );
        std::ostream &thOut = *OPT_thsOut[tid];
        
        thOut << "\tLeaving at optsolvers::OPT_mosekEval" << std::endl;
        
        OPT_closeFileToThread(tid);
    }
    #endif
    
    
    
    return code;
}


#endif














int OPT_Mosek::solve(const bool resetSol, const bool storeSol, const bool storeConstrs, const bool storeDualSol )
#if OPT_HAVE_MOSEK
{
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        OPT_createFileToThread( tid );
        std::ostream &thOut = *OPT_thsOut[tid];
        
        thOut << "Entering at OPT_Mosek::solve" << std::endl;
    }
    #endif
    
    int nint = 0;
    
    MSKrescodee  r;
    MSKrescodee trmcode;
    MSKsolstae solsta;
    
    /*MSK_writeparamfile(task, "parametros2.txt");
    std::cout << "Escrevi parametros2.txt";
    OPT_getchar(); */
    
    
    if(resetSol)
    {
        origSolverRetCode = MSK_SOL_STA_UNKNOWN;
        this->resetSol();
    }
    
    feasSol = false;
    
    /*//TODO:remove this if
    if( nlObj || hasNLConstrs )
    {
        MSKuserhandle_t data;
        MSKnlgetspfunc *f1 = NULL;
        MSKnlgetvafunc *f2 = NULL;
        
        
        MSK_getnlfunc(task, &data, f1, f2);
        
        if( f1 == NULL )
        {
            r = MSK_putnlfunc(task, this, OPT_mosekStruc, OPT_mosekEval);
            
            
            if( r != MSK_RES_OK )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                retCode = OPT_UNDEFINED_ERROR;
                goto termination;
            }
        }
    } */
    
    //printf("Entrando em MSK_optimizetrm Thread: %d\n", threadNumber);
    //fflush( stdout);
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        std::ostream &thOut = *OPT_thsOut[std::this_thread::get_id()];
        thOut << "\tOPT_Mosek::solve 1" << std::endl;
    }
    #endif
    
    /*#if OPT_SUPER_THREAD_DEBUG_MODE
    {
        std::stringstream dumpFileName;
        
        dumpFileName << "mosek_dump_file_" << std::this_thread::get_id() << ".tsk";
        
        MSK_writetask( task, dumpFileName.str().c_str() );
    }
    #endif */
    
    
    r = MSK_optimizetrm(task, &trmcode);
    
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        std::ostream &thOut = *OPT_thsOut[std::this_thread::get_id()];
        thOut << "\tOPT_Mosek::solve 2" << std::endl;
    }
    #endif
    
    
    //printf("Saindo de MSK_optimizetrm Thread: %d\n", threadNumber);
    //fflush( stdout);
    
    if( r != MSK_RES_OK )
    { 
        #if OPT_DEBUG_MODE
            #if OPT_PRINT_ERROR_RETURN_CODE_ON_SOLVE
                OPT_PRINTERRORNUMBER(r);
            #endif
        #endif
        
        
        if( OPT_CHECKMOSEKLICENSEERROR(r) )
        {
            OPT_PRINTERRORMSG("License error in mosek library! Check your license file."); ////we print the error even if we are not in debug mode. 
            retCode = OPT_LICENSE_ERROR;
        }
        else
        {
            retCode = OPT_UNDEFINED_ERROR;
        }
        
        origSolverRetCode = r;
        goto termination;
    }
    
    //termcode is checked below
    
    
    MSK_getnumintvar(task, &nint);
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        std::ostream &thOut = *OPT_thsOut[std::this_thread::get_id()];
        thOut << "\tOPT_Mosek::solve 3" << std::endl;
    }
    #endif
    
    
    
    if( nint > 0 )
    {
        MSK_getsolsta(task, MSK_SOL_ITG, &solsta);
        
        origSolverRetCode = solsta;
        
        getSolution(MSK_SOL_ITG, storeSol, storeConstrs, storeDualSol);
        
        
        switch (solsta)
        {
            case MSK_SOL_STA_INTEGER_OPTIMAL:
            case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
                
                
                
                retCode = OPT_OPTIMAL_SOLUTION;
                break;
                
            default:
                
                switch(trmcode)
                {
                    case MSK_RES_OK:
                        break;
                    
                    
                    case MSK_RES_TRM_MAX_TIME:
                        retCode = OPT_MAX_TIME;
                        goto termination;
                    
                    
                    case MSK_RES_TRM_MAX_ITERATIONS:
                        #if OPT_PRINT_MAX_ITER_WARNING
                            std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Mosek solving!\n";
                        #endif
                        
                        retCode = OPT_MAX_ITERATIONS;
                        goto termination;
                        
                        
                    default:
                        #if OPT_DEBUG_MODE
                            std::cout << OPT_PREPRINT "Problem " << trmcode << " on Mosek termination code at OPT_Mosek::solve\n";
                            //getchar();
                        #endif
                        //we do not go to because the solution can still be near to optimal
                        break;
                }
                
                
                retCode = OPT_UNDEFINED_ERROR;
        }
        
    }
    else
    {
        MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
        
        //std::cout << "trmcode: " << trmcode << " solsta: " << solsta << "\n";
        
        origSolverRetCode = solsta;
        
        getSolution(MSK_SOL_ITR, storeSol, storeConstrs, storeDualSol);
        
        switch(solsta)
        {
            case MSK_SOL_STA_OPTIMAL:
            case MSK_SOL_STA_NEAR_OPTIMAL:
                
                retCode = OPT_OPTIMAL_SOLUTION;
                feasSol = true;
                break;
                
                
            case MSK_SOL_STA_PRIM_INFEAS_CER:
                
                retCode = OPT_INFEASIBLE_PROBLEM;
                break;
                
            case MSK_SOL_STA_PRIM_FEAS:
            case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
                
                retCode = OPT_FEASIBLE_SOLUTION;
                feasSol = true;
                break;
                
            case MSK_SOL_STA_DUAL_INFEAS_CER:
                
                /*It is not working... :( Mosek sometimes return a array of constraints having feasible values, but problem is unbounded (??????)
                * {
                    feasSol = true;
                    int m;
                    double *b;
                    
                    MSK_getnumcon(task, &m);
                    
                    if( storeConstrs )
                    {
                        b = constr;
                    }
                    else
                    {
                        b = auxValues;
                        
                        const int r = MSK_getsolutionslice(task, MSK_SOL_ITR, MSK_SOL_ITEM_XC, 0, m, b);
                        if( r != MSK_RES_OK )
                        {
                            #if OPT_DEBUG_MODE
                                cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
                            #endif
                            feasSol = false;
                        }
                    }
                    
                    if( feasSol )
                    {
                        for(int i  = 0; i < m; i++)
                        {
                            double cilb, ciub;
                            
                            getConstraintBounds(i, cilb, ciub);
                            
                            if( b[i] < cilb || b[i] > ciub )
                            {
                                feasSol = false;
                                break;
                            }
                        }
                    }
                    
                    retCode = feasSol ? OPT_UNBOUNDED_PROBLEM : OPT_INFEASIBLE_PROBLEM; //maybe we have a unbounded problem and no feasible solution given by the solver 
                } */
                
                //retCode = OPT_UNBOUNDED_PROBLEM;
                retCode = OPT_INFEASIBLE_PROBLEM; //maybe the problem is unbounded, but we declare infeasibility...
                break;
                
            
            default:
                
                switch(trmcode)
                {
                    case MSK_RES_OK:
                        break;
                    
                    
                    case MSK_RES_TRM_MAX_TIME:
                        retCode = OPT_MAX_TIME;
                        goto termination;
                    
                    
                    case MSK_RES_TRM_MAX_ITERATIONS:
                        #if OPT_PRINT_MAX_ITER_WARNING
                            std::cerr << OPT_PREPRINT "Warning: Maximum iteration achieved on Mosek solving!\n";
                        #endif
                        
                        retCode = OPT_MAX_ITERATIONS;
                        goto termination;
                        
                        
                    default:
                        #if OPT_DEBUG_MODE
                            std::cout << OPT_PREPRINT "Problem " << trmcode << " on Mosek termination code at OPT_Mosek::solve\n";
                            //getchar();
                        #endif
                        //we do not go to because the solution can still be near to optimal
                        break;
                }
                
                
                retCode = OPT_UNDEFINED_ERROR;
        }
        
    }
    
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        std::ostream &thOut = *OPT_thsOut[std::this_thread::get_id()];
        thOut << "\tOPT_Mosek::solve 4" << std::endl;
    }
    #endif
    
    
termination:
    
    #if OPT_SUPER_THREAD_DEBUG_MODE
    {
        const std::thread::id &tid = std::this_thread::get_id();
        OPT_createFileToThread( tid );
        std::ofstream &thOut = * (std::ofstream *) OPT_thsOut[tid];
        
        thOut << "Leaving at OPT_Mosek::solve" << std::endl;
        
        //sometimes, we have to execute the so many times generating very large outuput files. So, if we reach this point, it is because this execution of solve was successfull and we can close the debug file to reopen it as an empty file in the next execution. In this way, the debug will not be large, and we will have the information that we want! ;)
        
        OPT_closeFileToThread(tid);
    }
    #endif
    
    
    return retCode;
    
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


#if OPT_HAVE_MOSEK
MSKboundkeye OPT_Mosek::boundKeye(const double lb, const double ub)
{
    MSKboundkeye auxBKeye;
    
    if( lb > -OPT_INFINITY )
    {
        if( ub < OPT_INFINITY )
            auxBKeye = lb == ub ? MSK_BK_FX : MSK_BK_RA;
        else
            auxBKeye = MSK_BK_LO;
    }
    else
    {
        auxBKeye = MSK_BK_UP;
    }
    
    return auxBKeye;
}
#endif


int OPT_Mosek::resetConstraintLinearPart( const int constrIndex, const int nzs, const int* cols, const double* values )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    
    
    /* r = MSK_putconbound(task, constrIndex, boundKeye(lb, ub), lb, ub);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    } */
    
    
    
    r = MSK_putarow(task, constrIndex, nzs, cols, values);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::__addConstraints(const int nm)
#if OPT_HAVE_MOSEK
{
    int code;
    MSKrescodee r;
    
    
    if(hasNLConstrs || nlObj)
    {
        //Mosek does not allow add constraint in a nonlinear problem... ridicuclous... 
        int r = MSK_putnlfunc(task, NULL, NULL, NULL);
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            code = OPT_UNDEFINED_ERROR;
            goto termination;
        }
    }
    
    
    r = MSK_appendcons(task, nm);
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code =  OPT_SOLVER_ERROR;
        goto termination;
    }
    
    
    code = 0;
    
termination:
    
    
    if(hasNLConstrs || nlObj)
    {
        int r = MSK_putnlfunc(task, this, OPT_mosekStruc, OPT_mosekEval);
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            code = OPT_UNDEFINED_ERROR;
        }
    }
    
    
    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::__addVariables(const int nn, const bool initFree)
#if OPT_HAVE_MOSEK
{
    const bool changeFunctionPointers = nlObj || hasNLConstrs;
    
    int code = 0;
    
    
    if( changeFunctionPointers )
    {
        //unfortunatelly, some mosek versins do not allow add variable if we have nonlinear functions...
        MSK_putnlfunc(task, NULL, NULL, NULL);
    }
    
    
    
    MSKrescodee r = MSK_appendvars(task, nn);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
        goto termination;
    }
    
    
    if( initFree )
    {
        int nvars = 0;
        MSK_getnumvar(task, &nvars);
        
        
        for(int i = nvars - nn; i < nvars; i++)
        {
            const MSKrescodee r = MSK_putvarbound(task, i, MSK_BK_FR, 0.0, 0.0);
            
            if( r != MSK_RES_OK )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR;
                //do not put goto here
            }
        }
        
    }
    
    
termination:
    
    if( changeFunctionPointers )
    {
        // restoring function pointers
        const int r = MSK_putnlfunc(task, this, OPT_mosekStruc, OPT_mosekEval);
        
        assert(r == MSK_RES_OK);
    }
    
    
    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::generateModelFile(const char *fileName)
#if OPT_HAVE_MOSEK
{
    const bool changeFunctionPointers = nlObj || hasNLConstrs;
    
    MSK_putobjname(task, "MosekOptSolvers");
    
    
    if( changeFunctionPointers )
    {
        std::cerr << OPT_PREPRINT "Warning: nonlinear functions set in mosek model. Ignoring these functions to gerenate model file " << fileName <<  ".\n" ;
        
        MSK_putnlfunc(task, NULL, NULL, NULL);
    }
    
    
    MSK_commitchanges(task);
    
    MSKrescodee r = MSK_writedata(task, fileName);
    
    if( changeFunctionPointers )
    {
        // restoring function pointers
        const int r = MSK_putnlfunc(task, this, OPT_mosekStruc, OPT_mosekEval);
        
        assert(r == MSK_RES_OK);
    }
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getConstraintBounds( const int index, double &lb, double &ub )
#if OPT_HAVE_MOSEK
{
    MSKboundkeye bk;
    
    lb = -OPT_INFINITY;
    ub = OPT_INFINITY;
    
    MSKrescodee r = MSK_getconbound(task, index, &bk, &lb, &ub);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getConstraintLinearCoef( const int constrIndex, const int varIndex, double &value)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getaij(task, constrIndex, varIndex, &value);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Mosek::getConstraintLinearPart(const int constrIndex, int &nzs, int *cols, double *values)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getarow( task, constrIndex, &nzs, cols, values );
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getObjLinearCoef( const int index, double &value )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getcj(task, index, &value);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getNumberOfConstraints(int& m)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getnumcon(task, &m);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Mosek::getNumberOfConstraintLinearCoefs( const int constrIndex, int &nzs)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getarownumnz(task, constrIndex, &nzs);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Mosek::getNumberOfVars(int& n)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getnumvar(task, &n);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getNumberOfIntVars(int &nI)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getnumintvar(task, &nI);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Mosek::getObjConstant(double &objConstant)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getcfix(task, &objConstant);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Mosek::getObjSense(OPT_OPTSENSE &sense)
#if OPT_HAVE_MOSEK
{
    MSKobjsensee s;
    MSKrescodee r = MSK_getobjsense(task, &s);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    sense = s == MSK_OBJECTIVE_SENSE_MINIMIZE ? OPT_MINIMIZE : OPT_MAXIMIZE ;
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Mosek::getVariableBounds(const int index, double &lb, double &ub)
#if OPT_HAVE_MOSEK
{
    MSKboundkeye bk;
    
    lb = -OPT_INFINITY;
    ub = OPT_INFINITY;
    
    
    //int n;
    //getNumberOfVars(n);
    //cout << "n: " << n << " index: " << index << endl;
    
    
    MSKrescodee r = MSK_getvarbound( task, index, &bk, &lb, &ub );
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::__removeConstraints(const int ninds, const int* indices )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_removecons(task, ninds, indices);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    
    if( hasNLConstrs )
        updatehasNLConstr();
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::removeVars(const int ninds, const int* indices )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_removevars(task, ninds, indices);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



void OPT_Mosek::setObjConstant(const double value)
{
#if OPT_HAVE_MOSEK
    MSKrescodee r = MSK_putcfix(task , value);
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            OPT_PRINTERRORNUMBER(r);
            //return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    //return 0;
#endif
}


int OPT_Mosek::setLinearColumn( const int varIndex, const int nzs, const int *rows, const double *values)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putacol(task, varIndex, nzs, rows, values);
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
    
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setConstraintBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_MOSEK
{
    MSKboundkeye auxBKeye;
    MSKrescodee r;
    
    
    if( lb > -OPT_INFINITY )
    {
        if( ub < OPT_INFINITY )
            auxBKeye = lb == ub ? MSK_BK_FX : MSK_BK_RA;
        else
            auxBKeye = MSK_BK_LO;
    }
    else
    {
        auxBKeye = MSK_BK_UP;
    }
    
    
    r = MSK_putconbound(task, index, auxBKeye, lb, ub);
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


/*
int OPT_Mosek::setConstraintLowerBound( const int index, const double lb )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putconbound(task, index, MSK_BK_LO, lb, OPT_INFINITY);
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::setConstraintUpperBound( const int index, const double ub )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putconbound(task, index, MSK_BK_UP, -OPT_INFINITY, ub);
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif
*/




int OPT_Mosek::setConstraintsLinearCoefs( const int nzs, const int *rows, const int *cols, const double *values )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putaijlist(task, nzs, rows, cols, values);
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::setConstraintLinearCoefs( const int constrIndex, const int nzs, const int* cols, const double* values)
#if OPT_HAVE_MOSEK
{
    OPT_setAllArray(nzs, auxIndex, constrIndex);
    
    
    //MSKrescodee r = MSK_putarow(task, constrIndex, nzs, cols, values);
    
    MSKrescodee r = MSK_putaijlist(task, nzs, auxIndex, cols, values);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setConstraintLinearCoef( const int constrIndex, const int varIndex, const double value)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putaij(task, constrIndex, varIndex, value);
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setObjLinearCoef( const int index, const double value )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putcj(task, index, value);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
    
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::setObjLinearCoefs( const int nzs, const int *cols, const double *values )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putclist(task, nzs, cols, values);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setObjLinearPart( const int n, const double *values )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    int code = 0, on;
    
    r = MSK_getnumvar(task, &on);
    
    if( on < n )
        return OPT_BAD_INPUT;
    
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    for(int i = 0; i < on; i++)
        auxIndex[i] = i;
    
    
    r = MSK_putclist( task, n, auxIndex, values );
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        code = OPT_SOLVER_ERROR;
    }
    
    
    if(n != on)
    {
        for(int i = n; i < on; i++)
            auxValues[i] = 0.0;
        
        r = MSK_putclist( task, on-n, &auxIndex[n], &auxValues[n] );
    
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            code = OPT_SOLVER_ERROR;
        }
    }
    
    
    
    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setObjSense( const OPT_OPTSENSE sense )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putobjsense(task, sense == OPT_MINIMIZE ? MSK_OBJECTIVE_SENSE_MINIMIZE : MSK_OBJECTIVE_SENSE_MAXIMIZE );
    
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            OPT_PRINTERRORNUMBER(r);
            
            return OPT_SOLVER_ERROR;
        }
    #endif
    
    
    return 0;
    
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setnVariablesBounds( const int n, const double *lb, const double *ub )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    for(int i = 0; i < n; i++)
        auxIndex[i] = boundKeye( lb[i], ub[i] );
    
    
    r = MSK_putvarboundslice( task, 0, n, (MSKboundkeye*) auxIndex, lb, ub);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setVariableBounds( const int index, const double lb, const double ub )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putvarbound( task, index, boundKeye(lb,ub), lb, ub );
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setVariablesBounds( const int ninds, const int *inds, const double *lb, const double *ub )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r;
    
    for(int i = 0; i < ninds; i++)
        auxIndex[i] = boundKeye(lb[i], ub[i]);
    
    r = MSK_putvarboundlist( task, ninds, inds, (MSKboundkeye*) auxIndex, lb, ub );
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



/*
int OPT_Mosek::setVariableLowerBound( const int index, const double lb )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putvarbound(task, index, MSK_BK_LO, lb, OPT_INFINITY);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::setVariableUpperBound( const int index, const double ub )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putvarbound(task, index, MSK_BK_UP, -OPT_INFINITY, ub);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            cerr << "optsolvers: Error " << r << OPT_GETFILELINE << endl;
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif
*/


int OPT_Mosek::getNumberOfQuadObjTerms(int &nzs)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getnumqobjnz( task, &nzs );
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getObjQuadTerm( const int row, const int col, double &value)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getqobjij( task, row, col, &value );
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getObjQuadPart( int &nzs, int *rows, int *cols, double *values )
#if OPT_HAVE_MOSEK
{
    int n;
    
    MSK_getnumvar(task, &n);
    
    int size = (n*(n+1))/2;
    
    MSKrescodee r = MSK_getqobj( task, size, &size, &nzs, rows, cols, values );
    
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setObjQuadCoef( const int row, const int col, const double value )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putqobjij(task, row, col, value);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::setObjQuadMatrix( const int nzs, const int* rows, const int* cols, const double* values )
#if OPT_HAVE_MOSEK
{
    //int nz = 0;
    MSKrescodee r;
    
    /*r = MSK_getnumqobjnz(task, &nz);
    
    #if OPT_DEBUG_MODE
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    #endif*/
    
    //cout << "mosek nz: " << nz << endl;
    
    //if( nz <= 0)
    {
        r = MSK_putqobj(task, nzs, rows, cols, values);
        
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
            
            return OPT_SOLVER_ERROR;
        }
    }
    /*else
    {
        int code = 0;
        
        for(int i = 0; i < nzs; i++)
        {
            r = MSK_putqobjij(task, rows[i], cols[i], values[i]);
            
            if( r != MSK_RES_OK )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                
                code = OPT_SOLVER_ERROR;
            }
        }
        
        if( code != 0 )
            return code;
    }*/
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif


int OPT_Mosek::getNumberOfConstraintQuadTerms( const int index, int &nzs)
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_getnumqconknz(task, index, &nzs);
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::getConstraintQuadMatrix( const int index, int &nzs, int *rows, int *cols, double *values )
#if OPT_HAVE_MOSEK
{
    int n;
    
    MSK_getnumvar( task, &n );
    
    
    int maxnzs = (n*(n+1))/2;
    
    MSKrescodee r = MSK_getqconk( task, index, maxnzs, &maxnzs, &nzs, rows, cols, values );
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::setConstraintQuadMatrix( const int index, const int nzs, const int* qrows, const int* qcols, const double* qvalues )
#if OPT_HAVE_MOSEK
{
    MSKrescodee r = MSK_putqconk( task, index, nzs, qrows, qcols, qvalues );
    
    if( r != MSK_RES_OK )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif
        
        return OPT_SOLVER_ERROR;
    }
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif






int OPT_Mosek::getConstrNLFlag(const int index, bool &flag)
{
#if OPT_DEBUG_MODE
    int m;
    getNumberOfConstraints(m);

    if( index < 0 || index >= m )
        return OPT_BAD_INPUT;
#endif

    flag = nlConstr[index];

    return 0;
}




int OPT_Mosek::getNumberOfNLConstraints(int &mnl)
{
    int r, m;

    mnl = 0;

    r = getNumberOfConstraints(m);
    if( r != 0 )
    {
        #if OPT_DEBUG_MODE
            OPT_PRINTERRORNUMBER(r);
        #endif

        return OPT_SOLVER_ERROR;
    }


    for(int i = 0; i < m; i++)
        if( nlConstr[i] )
            mnl++;

    return 0;
}



int OPT_Mosek::getObjNLFlag(bool &flag)
{
    flag = nlObj;
    return 0;
}



int OPT_Mosek::setConstrNLFlag(const int index, const bool flag)
#if OPT_HAVE_MOSEK
{
    const bool oldHasNLConstr = hasNLConstrs;
    bool check = false;
    
#if OPT_DEBUG_MODE
    int m;
    getNumberOfConstraints(m);

    if( index < 0 || index >= m )
        return OPT_BAD_INPUT;
#endif
    
    if( flag )
        hasNLConstrs = true;
    else if( nlConstr[index] ) //flag == false
        check = true;

    nlConstr[index] = flag;
    
    if( check )
        updatehasNLConstr();
    
    
    //if oldHasNLConstr is equal to oldHasNLConstr, we do not update NL function pointers in mosek
    if(hasNLConstrs != oldHasNLConstr)
    {
        if(!nlObj)
        {
            MSKnlgetspfunc pStruc;
            MSKnlgetvafunc pEval;
            
            if(hasNLConstrs)
            {
                pStruc = OPT_mosekStruc;
                pEval = OPT_mosekEval;
            }
            else
            {
                pStruc = NULL;
                pEval = NULL;
            }
            
            const int r = MSK_putnlfunc(task, this, pStruc, pEval);
            
            if( r != MSK_RES_OK )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif
                    
                return OPT_UNDEFINED_ERROR;
            }
        }
    }
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif



int OPT_Mosek::setInitialSolution(const double *x, const double *dualConstrs, const double *dualVars)
#if OPT_HAVE_MOSEK
{
    int n, m, code = 0;
    MSKrescodee r;
    
    
    if( x && dualVars )
    {
        
        r = MSK_getnumvar(task, &n);
        if( r != 0 )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif

            return OPT_SOLVER_ERROR;
        }
        
        
        
        const double *dualvlb = dualVars;
        const double *dualvub = &dualVars[n];
        
        for(int i = 0; i < n; i++)
        {
            r = MSK_putsolutioni( task, MSK_ACC_VAR, i, MSK_SOL_ITR, MSK_SK_UNK, x[i], dualvlb[i], dualvub[i], 0.0 );
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif

                code = OPT_SOLVER_ERROR;
            }
        }
        
    }
    
    if(dualConstrs)
    {
        r = MSK_getnumcon(task, &m);
        
        for(int i = 0; i < m; i++)
        {
            r = MSK_putsolutionyi( task, i, MSK_SOL_ITR, dualConstrs[i] );
            
            if( r != 0 )
            {
                #if OPT_DEBUG_MODE
                    OPT_PRINTERRORNUMBER(r);
                #endif

                code = OPT_SOLVER_ERROR;
            }
        }
    }
    
    
    
    return code;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif




int OPT_Mosek::setProblemFrom(const minlpproblem::MIP_MINLPProb &prob, const bool setObj, const bool setConstrs, const bool setVarBounds, const bool setVarType, const int naddvars, const int naddconstrs )
{
    const int r = OPT_NLPSolver::setProblemFrom( prob, setObj, setConstrs, setVarBounds, setVarType, naddvars, naddconstrs);
    
    nlObjFactor = prob.objFactor;
    
    return r;
}




int OPT_Mosek::setObjNLFlag(const bool flag)
#if OPT_HAVE_MOSEK
{
    //if hasNLConstrs is true, we do not need worry about MSK_putnlfunc (it is already non NULL)
    if( !hasNLConstrs )
    {
        int r;
        MSKnlgetspfunc pStruc;
        MSKnlgetvafunc pEval;
        
        if(flag)
        {
            pStruc = OPT_mosekStruc;
            pEval = OPT_mosekEval;
        }
        else
        {
            pStruc = NULL;
            pEval = NULL;
        }
        
        r = MSK_putnlfunc(task, this, pStruc, pEval);
        if( r != MSK_RES_OK )
        {
            #if OPT_DEBUG_MODE
                OPT_PRINTERRORNUMBER(r);
            #endif
                
            return OPT_UNDEFINED_ERROR;
        }
    }
    
    nlObj = flag;
    
    
    return 0;
}
#else
{
    OPT_LIBNOTAVAILABLERET(getSolverCode());
}
#endif









