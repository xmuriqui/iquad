/*
* That library defines functions and data structures used in the development of iquad. 
* In general, users do not need include it in their applications.
* 
*/

#ifndef _IQD_IQUAD_TOOLS_HPP
#define _IQD_IQUAD_TOOLS_HPP

#include "iquad.hpp"
#include <random>






namespace iquad
{
    
    #define IQD_STRATT(att)   #att, att
    
    
    #define IQD_PREPRINT "iquad: "
    
    
    #ifdef __FILE__
        #ifdef __LINE__
            #define IQD_DEF_GETFILELINE 1
        #endif
    #endif
    
    
    #ifdef IQD_DEF_GETFILELINE

        #define IQD_GETFILELINE  \
            " on file: " << __FILE__ << " line: " << __LINE__
    #else
        #define IQD_GETFILELINE ""
    #endif
    
    
    #define IQD_PRINTMEMERROR std::cerr << IQD_PREPRINT << "Memory error" << IQD_GETFILELINE << std::endl


    #define IQD_PRINTERROR std::cerr << IQD_PREPRINT << "Error" << IQD_GETFILELINE << std::endl


    #define IQD_PRINTERRORNUMBER(number) std::cerr << IQD_PREPRINT << "Error " << number << IQD_GETFILELINE << std::endl


    #define IQD_PRINTERRORMSG(msg) std::cerr << IQD_PREPRINT << msg << IQD_GETFILELINE << std::endl


    #define IQD_getchar()  IQD_PRINTERRORMSG("stopped in a getchar"), getchar() 
    
    
    
    
    
    
    
    #define IQD_SAVE_OUTPUT_FILE    1

    #if IQD_SAVE_OUTPUT_FILE
        #define IQD_OUT_FILE_NAME	"iquad_output.txt"
        #define IQD_CHAR_SEP	"#"
    #endif
    
    
    
    
    #define IQD_SAVE_OUTPUT_FILE_BY_HOUR 1
    
    #if IQD_SAVE_OUTPUT_FILE_BY_HOUR
        #define IQD_OUT_FILE_BY_HOUR_NAME "iquad_output_hours.txt"
    #endif
    
    
    
    
    class IQD_OutFile
    {
        FILE *outFile;
        
    public:
        
        double lastHoursToSaveFile;
        
        IQD_OutFile();
        
        ~IQD_OutFile();
        
        int open(const char* fileName);
        
        int writeProblemInfo(const char* name, IQD_IQuadProb &prob);
        
        int writeSolInfo(const int code, const double lb, const double objF, const double cpu_time, const double time, const long unsigned int niters, const double split_time, const double root_node_lb, const double root_node_ub);
        
        int writeDouble(const double value);
        
        int writeInt(const int value);
        
        int writeString(const char *value);
        
        int writeBreakLine();
        
        int writeToCompleteHours( const int ndata, const double maxCpuTime, const double hoursToSaveFile );
        
        int flush();
        
        void close();
    
    };
    
    
    
    class IQD_SDPLambda
    {
    
    public:
        
        int calculateLambda(const int sdpSolver, IQD_SparseMatrix& Q, const double Qfactor,  const int nv, double** v, double* lambda, const double* objCoeflambda, IQD_GeneralSolverParams* sdpSolverParams, const double zero_tol);//, const double eps_diag);
        
    private:
        
        int calculateLambdaMosek(IQD_SparseMatrix& Q, const double Qfactor, const int nv, double** v, double* lambda, const double* objCoefLambda, iquad::IQD_GeneralSolverParams* sdpSolverParams);
        
        int calculateLambdaCSDP(IQD_SparseMatrix& Q, const double Qfactor, const int nv, double** v, double* lambda, const double* objCoefLambda);
    };
    
    
    
    class IQD_AuxBoundsNLEval : public optsolvers::OPT_NonLinearEval
    {
        
    public:
        
        int indobjcut;
        //double *values; //commented this code in v0.2.2b
        
        optsolvers::OPT_NonLinearEval *reval;
        IQD_IQuadProb *rprob;
        
        
        virtual int initialize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess);
        
        
        virtual int eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value);
        
        
        virtual int eval_nl_constrs_part(const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values);
        
        
        virtual int eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values);
        
        
        virtual int eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, minlpproblem::MIP_SparseMatrix& jacobian);
        
        
        virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, minlpproblem::MIP_SparseMatrix& hessian);
        
        
        virtual void finalize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess);
        
        
        virtual ~IQD_AuxBoundsNLEval();
        
    };
    
    
    
    
    
    
    class IQD_AuxBoundsUpdate
    {
        int indobjcut;
        int indsecconstrs;
        
        int *cols;
        double *vals;
        
        IQD_AuxBoundsNLEval *nleval;
        
    public:
        
        optsolvers::OPT_Solver *solver;
        
        IQD_AuxBoundsUpdate();
        
        ~IQD_AuxBoundsUpdate();
        
        int allocateAuxArrays( const int size );
        
        
        //note, rprob has the reformulated constraint, auxiliary bounds constraints, but does not have the secant constraints...
        int buildPreProblem(iquad::IQD_IQuadProb &rprob, const int dimv, const int qpSolver, const int qcpSolver, const int nlpSolver);
        
        int calculateNewBounds( const int norigvars, const int dimv, const int indvconstrs, iquad::IQD_SparseMatrix& v, const double zu, const double* ly, const double* uy, double* newly, double* newuy );
        
        void desallocateAuxArrays();
        
        void desallocateSolver();
        
        void initialize();
    };
    
    
    
    
    class IQD_RefProb2NonLinearEval : public minlpproblem::MIP_NonLinearEval
    {
        
    public:
        
        int nthreads;
        int nChgConstr;
        
        bool **constrEvals;
        int *indChgConstr;
        double **lambdas;
        IQD_RefProb2 *refProb;
        
        
        IQD_RefProb2NonLinearEval(IQD_RefProb2 *refProb = NULL);
        
        virtual ~IQD_RefProb2NonLinearEval();
        
        int allocate(const int nthreads);
        
        void desallocate();
        
        void initialize(IQD_RefProb2 *refProb);
        
    public:
        
        virtual int eval_nl_obj_part(const int thnumber, const int n, const bool newx, const double* x, double& value);
        
        virtual int eval_nl_constrs_part(const int thnumber, const int n, const int m, const bool newx, const bool* constrEval, const double* x, double* values);
        
        virtual int eval_grad_nl_obj_part(const int thnumber, const int n, const bool newx, const double* x, double* values);
        
        virtual int eval_grad_nl_constrs_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, minlpproblem::MIP_SparseMatrix& jacobian);
        
        virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, minlpproblem::MIP_SparseMatrix& hessian);
        
    };
    
    
    
    
    class IQD_RefProb2
    {
        
    protected:
        
        //double *auxv;
        int dimLambdaConstr;
        
    public:
        
        IQD_SPLITQUAD splitWay;
        
        IQD_IQuadProb *oprob;
        IQD_IQuadProb rprob; //reformulated problem
        
        bool nonConvexConstrs; //true if original problem has noncovex quadratic constraints...
        
        int dimv;
        
        int indvconstrs; //index where directions bound constraints on v start.
        
        int mdis; //number of dismembered constraints, i.e., double bound constraints that needed be dismembered in two.
        
        bool *chgSignal; //to 0 <= i < m+mdis, chgSignal[i] is true if constraint i had its signal changed. (It is usefull to know for waht constraints we should multiply nonlinear term by (-1) )
        
        int *origConstr; //to associate the new constraints with the original constraints... to 0 <= i < m, origConstr[i] has the index of new constraint associated to original constraint i. To m <= i < m+mdis, origConstr[i] it has the index of original constraint associated to new constraint i
        
        double *lambdaObj; //vector used to give weights to vectors v in decomposition of R. Dimension is dimv
        double *lambdaObjNoInt; //vector used to give weights to vector v in dcecompisiton of R. Here we removed positions about integer vars...
        double **lambdaConstr; //vector lambda for each constraint used to decompose Q_i ;
        IQD_SparseMatrix v; //double **v; //n-vectors used to decomposed R matrices in the problem. The number of vectors is dimv.
        
        IQD_SparseMatrix *PnoInt;
        
        double *ly, *uy; //bounds to new y variables to do branching
        
        double eps_diag_matrix_from_sdp;
        double abs_eps_diag_P;
        double rel_eps_diag_P;
        double zero_tol;
        double max_cpu_time, max_time;
        
        
        IQD_RefProb2();
        
        ~IQD_RefProb2();
        
        
        
        int allocateOrigConstrArrays(const int size);
        
        //allocate v, lambdaObj, lambdaConstr
        int allocateVAndLambda(const int nconstr, const int dimv, const int cols );
        
        void calculatey(const double *x, double *y);
        
        void desallocate();
        
        unsigned int getNDual();
        
        void initialize();
        
        //that method tries to calculate weak bounds to auxiliary variables y to Branch-And-Bound
        int calculateAuxBoundsToBB();
        
        
        //we perform the split on indefinite quadratic matrices and add auxiliary variables in the constraints, but we do not set the auxiliary constraints having bounds to y and secant inequalities...
        int buildRefProb( bool splitWithouIntVars, const bool certainOfNonconvexQ, const int sdpSolver, const int qpSolver, const int qcpSolver, const int nlpSolver, iquad::IQD_GeneralSolverParams* sdpSolverParams, iquad::IQD_GeneralSolverParams* solverParams );
        
        
        //auxColumns and auxValues are arrays of size n to help computations..
        //int setSecantConstraints( const int indsecconstrs, optsolvers::OPT_LPSolver *solver, const double *ly, const double *uy, const bool equality, int *auxColumns, double *auxValues );
        
        
        
    
    protected:
        
        
        int calculateBoundsOnDirections( const double max_cpu_time, const double max_time, const int qpSolver, const int qcpSolver, const int nlpSolver, IQD_GeneralSolverParams *solverParams, const int ndirs, double **dirs, bool *lbconc, bool *ubconv, double *lbdir, double *ubdir );
        
        int calculateLambdaByDiagDominant( IQD_SparseMatrix& Q, const double factor, const double epsilon, double *lambda);
        
        int calculateLambdaByIdentity( IQD_SparseMatrix& Q, const double factor, double *lambda );
        
        void calculateRbyLambdaAndv(const int n, const int dimv, double **v, const double *lambda, double *R, const bool initializeRwithZeros);
        
        //lambda is a output array
        int splitMatrix(const int splitWay, IQD_SparseMatrix& Q, const double factor, const int nvecs, double** vecs, double* vecsWeight, int sdpSolver, IQD_GeneralSolverParams* sdpSolverParams,  double* lambda, double* R, const bool initializeRwithZeros = true);
        
        
        
    };
    
    
    
    
    
    class IQD_Random
    {
        std::mt19937 gen;
    
        long int seed;
        
    public:
        
        IQD_Random();
        
        IQD_Random(const int long seed);
        
        //~IQD_Random();
        
        //if NULL, current time is used like seed...
        long int setSeed(const long int* seed = NULL);
        
        //generates a random integer in a interval [begin  end]
        int randInt(const int begin, const int end);
        
        
        //generates true with probability prob
        bool randBool(const double prob);
        
        //generates a uniform random real in the interval [0 1)
        double random();
        
        //generates a normal random real in the interval [0 1)
        double randomNormal();
        
        //generates a random real in the interval [begin end)
        double random(const double begin, const double end);
        
        //generates a random real in the interval [begin end)
        double randomNormal(const double begin, const double end);
        
    };
    
    
    struct IQD_itemToSort
    {
        unsigned int index;
        double value;
    };
    
    
    
    template <class myClass>
    inline void IQD_copyArray(const unsigned int n, const myClass *source, myClass* const destination)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < n; i++)
            destination[i] = source[i];
    };
    
    
    template <class myClass>
    inline void IQD_setAllarray( const unsigned int size, myClass *a, myClass value)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] = value;
    }
    
    
    template <class myClass>
    inline void IQD_multiplyAllArray( const unsigned int size, myClass *a, myClass value)
    {
        #pragma ivdep
        #pragma GCC ivdep
        for(unsigned int i = 0; i < size; i++)
            a[i] *= value;
    }
    
    
    
    template <class myClass>
    inline void IQD_secFree(myClass* &p)
    {
        if(p)
        {
            free(p);
            p = nullptr;
        }
    }


    template <class myClass>
    inline void IQD_secDelete(myClass* &p)
    {
        if(p)
        {
            delete p;
            p = nullptr;
        }
    }
    
    template <class myClass>
    inline void IQD_secDeleteArray(myClass* &p)
    {
        if(p)
        {
            delete[] p;
            p = nullptr;
        }
    }
    
    
    template <class myClass1, class myClass2>
    inline void IQD_swap( myClass1 &v1, myClass2 &v2 )
    {
        myClass1 vt = v1;
        v1 = v2;
        v2 = vt;
    }
    
    
    template <class myClass>
    inline void IQD_printMatrix( const unsigned int m, const unsigned int n, const myClass *M )
    {
        for(unsigned int i = 0; i < m; i++)
        {
            for(unsigned int j = 0; j < n; j++)
                std::cout << M[i*n + j] << " ";
            
            std::cout << ";" << std::endl;
        }
    }
    
    
    
    //set item if itemName matches with userName
    template <class myClass>
    inline int IQD_setEnum( const char *itemName, const myClass itemValue, const char *userName, myClass &item )
    {
        if( strcmp(itemName, userName) == 0 )
        {
            item = itemValue;
            return 0;
        }
        
        return IQD_NAME_ERROR;
    }
    
    
    
    //set att to userValue if attName matches with userName
    template <class myClass>
    inline int IQD_setAtt( const char *attName, myClass &att, const char *userName, myClass  userValue  )
    {
        if( strcmp(attName, userName) == 0)
        {
            att = userValue;
            return 0;
        }
        
        return IQD_NAME_ERROR;
    }


    template <class myClass>
    inline int IQD_setStrAtt( const char *attName, myClass &att, const char *userName, const char *userValue )
    {
        if( strcmp(attName, userName) == 0 )
        {
            const int r = IQD_getStrEnumValue(userValue, att);
            
            return r == 0 ? 0 : 1; //if we got a error in the value, we return a positive value to diferentiate of the case where attribute name is wrong
        }
        
        return IQD_NAME_ERROR;
    }
    
    
    
    
    
    
    
    inline bool IQD_isIntVarsFixed( const int nI, const int *intVars, const double *lb, const double *ub ) 
    {
        
        for(int i = 0; i < nI; i++)
        {
            const int ind = intVars[i];
            
            if( lb[ind] != ub[ind] )
                return false;
        }
        
        return true;
    }
    
    
    
    
    inline void IQD_addEpsOnNonZeroDiag( unsigned int n, double *M, const double zero_tol, const double absEps, const double relEps )
    {
        
        for(unsigned int i = 0; i < n; i++)
        {
            const int diag = i*n + i;
            
            if( IQD_abs( M[diag] ) > zero_tol )
                M[diag] += absEps + M[diag]*relEps; //we use M[diag] because in some cases we are subtrating for -P, i.e., subtrating from a negative semidefinite matrix
        }
    }
    
    

    int IQD_itemAscComparison(const void *v1, const void*v2);
    
    int IQD_itemDescComparison(const void *v1, const void*v2);
    
    inline void IQD_sortItens(const bool ascendent, const int nitens, struct IQD_itemToSort *itens )
    {
        qsort(itens, nitens, sizeof(struct IQD_itemToSort), ascendent ? IQD_itemAscComparison : IQD_itemDescComparison );
    }
    
    
    
    inline bool IQD_isDiagonalSplitQWay(const int way)
    {
        return way == IQD_SQ_DIAG_SDP || way == IQD_SQ_DIAG_DOM || way == IQD_SQ_IDENTITY ; //|| way == IQD_SQ_USER_DIAGONAL;
        
        //return way != IQD_SQ_SDP && way != IQD_SQ_SCHUR && way != IQD_SQ_USER_MATRIX && way != IQD_SQ_2_BLOCK_SDP && way != IQD_SQ_MILP_BASIS && way !=IQD_SQ_MIN_EIG_SDP;
    }

    inline bool IQD_isInteger(const int type)
    {
        //return type == IQD_VT_INTEGER || type == IQD_VT_BINARY;
        return minlpproblem::MIP_isIntegerType(type);
    }
    
    inline bool IQD_isIntegerSol(const int nI, const int *intVars, const double *sol, const double intTolerance)
    {
        for(int i = 0; i < nI; i++)
        {
            const double sind = sol[ intVars[i] ];
            
            if( IQD_abs(round(sind) - sind) > intTolerance )
                return false;
        }
        
        return true;
    }
    
    //warning: we assume n > 0
    inline double IQD_vectorTimes(const int n, const double *v1, const double *v2)
    {
        int i;
        double r = v1[0]*v2[0];
        
        for(i = 1; i < n; i++)
            r += v1[i]*v2[i];
        
        return r;
    }
    
    
    inline static double IQD_intGap( const double v )
    {
        return IQD_abs(round(v) - v);
    }
    
        
    int IQD_getEigenValues(const int n, const double *M, double* eig);
    
    
    int IQD_getEigenValues(const iquad::IQD_SparseMatrix& Q, double* eig);
    
    
    double IQD_norm2(const int n, const double *v);
    
    int IQD_writeREigInfo(const char* probName, const int n, const double* R, iquad::IQD_OutFile& outFile, const double zero_tol);

    int IQD_writeEigInfo(const char* probName, const iquad::IQD_SparseMatrix& Q, const char* fileName, const double zero_tol);

    void IQD_generateRstatistic(const char* probName, iquad::IQD_IQuadProb& prob, const iquad::IQD_SPLITQUAD splitWay, const int sdpSolver, const int milpSolver, const int qcpSolver, const int nlpSolver, const int minlpSolver, iquad::IQD_GeneralSolverParams* sParams, iquad::IQD_OutFile& outFile, const double zero_tol);

    double IQD_maxDifferenceBetweenMatrices(const int nrows, const int ncols, const double *M1, const double *M2);

    int IQD_solveLinearSystem(const int n, const double* A, double* b, double* x);

    //that function invert negative eigenvalues except the lowest in matrix Q.
    int IQD_inverteQeigenvalues(IQD_IQuadProb& prob, const double zero_tol);
    
    int IQD_changeQeigenvalues( IQD_IQuadProb &prob, const int nNegEig, const double zero_tol);
    
    int IQD_replaceZeroEigValuesOfQByRandom(IQD_IQuadProb &prob, const double zero_tol);


    int IQD_ampl(char* stub, char** argvo, const char* outDirectory);

    int IQD_biqmac(const char* fileName, const char* outDirectory, const bool considerIntegrality, const bool maxcut = false);

    int IQD_boxqp(const char* fileName, const char* outDirectory);
    
    int IQD_myqp(const char* fileName, const char* outDirectory);
    
    int IQD_readQuadProg(const char* fileName, const char* outDirectory);

    int IQD_schurDecomposition(const double* A, int dim, double* T, double* Zt, double* eigenValues);

    double IQD_getTime(void);

    void IQD_welcome(void);
    
    
    int IQD_setSecantConstraints( const int n, const int dimv, const IQD_SparseMatrix &v, const int indsecconstrs, optsolvers::OPT_LPSolver *solver, const double *ly, const double *uy, const bool equality, int *auxColumns, double *auxValues );
    
    int IQD_runProblem( const char *probName, const char *outDirectory, iquad::IQD_IQuadProb &prob, iquad::IQD_BranchAndBound &bb, iquad::IQD_GeneralSolverParams *solverParams, iquad::IQD_GeneralSolverParams *sdpSolverParams);
    
    
    //return true if quadratic part of problem is convex... (note, that function does NOT look to general nonlinear part). auxLambda is a optional argument. It is an array n x 1 having space to perform computations.
    int IQD_isConvexProblemAboutQuadratics( const iquad::IQD_IQuadProb& prob, const double zero_tol, bool& convexObj, bool& convexConstr );
    
    
    /*
    * Function to transform iquad problems second Jon's sugestions:
    * 
    * "i believe was to scale the box to substitute (1/10)x_j for x_j (everywhere), so that the problems are redefined on [0,10]^n rather than [0,1]^n. And then take the first n/2 variables and constrain them to be integer."
    * 
    * Here, we ignore nonlinear terms. Scale factor is a factor to increase (or reduce) variable box sizes...
    * 
    * 
    */
    int IQD_scaleProblemAndTurnIntegerVars( IQD_IQuadProb &prob, const double scaleFactor );
    
    
    
}


#endif

