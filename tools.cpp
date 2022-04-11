
#include <math.h>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <climits>


#include <new>

#include "iquad.hpp"
#include "IQD_tools.hpp"
#include "BBL_tools.hpp"


#if IQD_HAVE_LAPACK
extern "C" {
    
        
    /*
    * *
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
    CHARACTER          JOBVS, SORT
    INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
    LOGICAL            BWORK( * )
    DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
    $                   WR( * )
*     ..
*     .. Function Arguments ..
    LOGICAL            SELECT
    EXTERNAL           SELECT
*     ..
*
*  Purpose
*  =======
*
*  DGEES computes for an N-by-N real nonsymmetric matrix A, the
*  eigenvalues, the real Schur form T, and, optionally, the matrix of
*  Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
*
*  Optionally, it also orders the eigenvalues on the diagonal of the
*  real Schur form so that selected eigenvalues are at the top left.
*  The leading columns of Z then form an orthonormal basis for the
*  invariant subspace corresponding to the selected eigenvalues.
*
*  A matrix is in real Schur form if it is upper quasi-triangular with
*  1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
*  form
*          [  a  b  ]
*          [  c  a  ]
*
*  where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
*
*  Arguments
*  =========
*
*  JOBVS   (input) CHARACTER*1
*          = 'N': Schur vectors are not computed;
*          = 'V': Schur vectors are computed.
*
*  SORT    (input) CHARACTER*1
*          Specifies whether or not to order the eigenvalues on the
*          diagonal of the Schur form.
*          = 'N': Eigenvalues are not ordered;
*          = 'S': Eigenvalues are ordered (see SELECT).
*
*  SELECT  (external procedure) LOGICAL FUNCTION of two DOUBLE PRECISION arguments
*          SELECT must be declared EXTERNAL in the calling subroutine.
*          If SORT = 'S', SELECT is used to select eigenvalues to sort
*          to the top left of the Schur form.
*          If SORT = 'N', SELECT is not referenced.
*          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
*          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
*          conjugate pair of eigenvalues is selected, then both complex
*          eigenvalues are selected.
*          Note that a selected complex eigenvalue may no longer
*          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
*          ordering may change the value of complex eigenvalues
*          (especially if the eigenvalue is ill-conditioned); in this
*          case INFO is set to N+2 (see INFO below).
*
*  N       (input) INTEGER
*          The order of the matrix A. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.
*          On exit, A has been overwritten by its real Schur form T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  SDIM    (output) INTEGER
*          If SORT = 'N', SDIM = 0.
*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*                         for which SELECT is true. (Complex conjugate
*                         pairs for which SELECT is true for either
*                         eigenvalue count as 2.)
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          WR and WI contain the real and imaginary parts,
*          respectively, of the computed eigenvalues in the same order
*          that they appear on the diagonal of the output Schur form T.
*          Complex conjugate pairs of eigenvalues will appear
*          consecutively with the eigenvalue having the positive
*          imaginary part first.
*
*  VS      (output) DOUBLE PRECISION array, dimension (LDVS,N)
*          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
*          vectors.
*          If JOBVS = 'N', VS is not referenced.
*
*  LDVS    (input) INTEGER
*          The leading dimension of the array VS.  LDVS >= 1; if
*          JOBVS = 'V', LDVS >= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,3*N).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  BWORK   (workspace) LOGICAL array, dimension (N)
*          Not referenced if SORT = 'N'.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*          > 0: if INFO = i, and i is
*             <= N: the QR algorithm failed to compute all the
*                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
*                   contain those eigenvalues which have converged; if
*                   JOBVS = 'V', VS contains the matrix which reduces A
*                   to its partially converged Schur form.
*             = N+1: the eigenvalues could not be reordered because some
*                   eigenvalues were too close to separate (the problem
*                   is very ill-conditioned);
*             = N+2: after reordering, roundoff changed values of some
*                   complex eigenvalues so that leading eigenvalues in
*                   the Schur form no longer satisfy SELECT=.TRUE.  This
*                   could also be caused by underflow due to scaling.
*
*  =====================================================================
    */ 
    
    void dgees_(char* JOBVS, char* SORT, int (*SELECT)(double *a, double *b), int *N, double *A, int *LDA, int *SDIM, double *WR, double *WI, double *VS, int *LDVS, double *WORK, int *LWORK, int *BWORK, int *INFO );
    
    
    /*
    *
    **
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
    DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,3*N-1).
*          For optimal efficiency, LWORK >= (NB+2)*N,
*          where NB is the blocksize for DSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
    * 
    */ 
    

    void dsyev_(char *JOBZ, char *UPLO, int *N, double *A, int *LDA, double *W, double *WORK, int *LWORK, int *INFO );
    
    /*
    *
    * SUBROUTINE DSGESV( N, NRHS, A, LDA, IPIV, B, LDB, X, LDX, WORK,
    $                   SWORK, ITER, INFO )
*
*  -- LAPACK PROTOTYPE driver routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*  -- April 2011                                                      --
*
*     ..
*     .. Scalar Arguments ..
    INTEGER            INFO, ITER, LDA, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
    INTEGER            IPIV( * )
    REAL               SWORK( * )
    DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( N, * ),
    $                   X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DSGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  DSGESV first attempts to factorize the matrix in SINGLE PRECISION
*  and use this factorization within an iterative refinement procedure
*  to produce a solution with DOUBLE PRECISION normwise backward error
*  quality (see below). If the approach fails the method switches to a
*  DOUBLE PRECISION factorization and solve.
*
*  The iterative refinement is not going to be a winning strategy if
*  the ratio SINGLE PRECISION performance over DOUBLE PRECISION
*  performance is too small. A reasonable strategy should take the
*  number of right-hand sides and the size of the matrix into account.
*  This might be done with a call to ILAENV in the future. Up to now, we
*  always try iterative refinement.
*
*  The iterative refinement process is stopped if
*      ITER > ITERMAX
*  or for all the RHS we have:
*      RNRM < SQRT(N)*XNRM*ANRM*EPS*BWDMAX
*  where
*      o ITER is the number of the current iteration in the iterative
*        refinement process
*      o RNRM is the infinity-norm of the residual
*      o XNRM is the infinity-norm of the solution
*      o ANRM is the infinity-operator-norm of the matrix A
*      o EPS is the machine epsilon returned by DLAMCH('Epsilon')
*  The value ITERMAX and BWDMAX are fixed to 30 and 1.0D+00
*  respectively.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array,
*          dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, if iterative refinement has been successfully used
*          (INFO.EQ.0 and ITER.GE.0, see description below), then A is
*          unchanged, if double precision factorization has been used
*          (INFO.EQ.0 and ITER.LT.0, see description below), then the
*          array A contains the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*          Corresponds either to the single precision factorization
*          (if INFO.EQ.0 and ITER.GE.0) or the double precision
*          factorization (if INFO.EQ.0 and ITER.LT.0).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          The N-by-NRHS right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          If INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N,NRHS)
*          This array is used to hold the residual vectors.
*
*  SWORK   (workspace) REAL array, dimension (N*(N+NRHS))
*          This array is used to use the single precision matrix and the
*          right-hand sides or solutions in single precision.
*
*  ITER    (output) INTEGER
*          < 0: iterative refinement has failed, double precision
*               factorization has been performed
*               -1 : the routine fell back to full precision for
*                    implementation- or machine-specific reasons
*               -2 : narrowing the precision induced an overflow,
*                    the routine fell back to full precision
*               -3 : failure of SGETRF
*               -31: stop the iterative refinement after the 30th
*                    iterations
*          > 0: iterative refinement has been sucessfully used.
*               Returns the number of iterations
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) computed in DOUBLE PRECISION is
*                exactly zero.  The factorization has been completed,
*                but the factor U is exactly singular, so the solution
*                could not be computed.
    */
    
    
    void dsgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV, double *B, int *LDB, double *X, int *LDX, double *WORK, float *SWORK, int *ITER, int *INFO );
    
}
#endif




using namespace std;
using namespace minlpproblem;
using namespace optsolvers;
using namespace iquad;






IQD_Random::IQD_Random()
{
}


IQD_Random::IQD_Random(const long int seed)
{
    setSeed(&seed);
}



long int IQD_Random::setSeed(const long int *seed)
{   
    this->seed = seed ? *seed : time(NULL) ;
    
    gen.seed(this->seed);
    //this->seed = s;
    
    return this->seed;
}


//generates a random integer in a interval [begin  end]
int IQD_Random::randInt(const int begin, const int end)
{
    std::uniform_int_distribution<int> dist(begin, end);
    
    
    return dist(gen);
}


//generates true with probability prob
bool IQD_Random::randBool(const double prob)
{
    return random() <= prob;
}


//generates a random real in the interval [0 1)
double IQD_Random::random()
{
    static std::uniform_real_distribution<double> distUni(0.0, 1.0);
    
    return distUni(gen);
}


//generates a random real in the interval [begin end)
double IQD_Random::random(const double begin, const double end)
{
    return begin + (end - begin)*random();
}



//generates a normal random real in the interval [0 1)
double IQD_Random::randomNormal()
{
    static std::normal_distribution<double> distNormal(0.0, 1.0);
    
    return distNormal(gen);
}


//generates a random real in the interval [begin end)
double IQD_Random::randomNormal(const double begin, const double end)
{
    return begin + (end - begin)*randomNormal();
}




//we wanting order in descendent order. Due to it, we inverted the return values
int iquad::IQD_itemDescComparison(const void *v1, const void*v2)
{
    if( ((struct IQD_itemToSort *) v1)->value < ((struct IQD_itemToSort *) v2)->value  )
        return 1;
    else if( ((struct IQD_itemToSort *) v1)->value == ((struct IQD_itemToSort *) v2)->value  )
        return 0;
    else
        return -1;
}


int iquad::IQD_itemAscComparison(const void *v1, const void*v2)
{
    if( ((struct IQD_itemToSort *) v1)->value < ((struct IQD_itemToSort *) v2)->value  )
        return -1;
    else if( ((struct IQD_itemToSort *) v1)->value == ((struct IQD_itemToSort *) v2)->value  )
        return 0;
    else
        return 1;
}



//pay attention: Zt will be returned as the Schur vectors matrix transposed.
int iquad::IQD_schurDecomposition(const double* A, int dim, double* T, double* Zt, double* eigenValues)
#if IQD_HAVE_LAPACK
{
    int i, code;
    
    //input parameters to dgees_ function
    int LWORK= 10*dim;
    
    //output parameters to dgees_ function
    int SDIM;
    double  *WR, *WI, *WORK;

    WR = eigenValues;
    WI = (double *) malloc( dim * sizeof(double) );
    WORK=(double *) malloc(LWORK* sizeof(double) );
    
    if(!WI || !WORK)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    /*for(i = 0; i < dim; i++)
    {
    WR[i] = 0.0;
    }*/
    
    
    //dgees_ uses the input as output too. So, if we have to copy the values...
    for(i = 0; i < dim *dim; i++)
        T[i] = A[i];
    
    dgees_((char *)"V", (char *)"N", NULL, &dim, T, &dim, &SDIM, WR, WI, Zt, &dim, WORK, &LWORK, NULL, &code);
    
    //The output argument Z is a matriz following the fortran convention. (stored by columns). In that columns, we have the schur vectors. So, in that C++ code, note that the schur vectors are stored in the lines of Z, and it is advantageous for us. Due to it, we do not transpose Z and work on it in its actual format.
    
    /*
    printf("Saida dgees_: %d\n", code);
    printf("T:\n");
    for(i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
            printf("%f ", T[i*dim + j] );
        printf(";\n");
    }
    
    printf("Z:\n");
    for(i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
            printf("%f ", Zt[i*dim + j]);
        
        printf(";\n");
    }
    
    for(i = 0; i < dim; i++)
    printf("WR[%d]: %f WI[%d]: %f\n", i, WR[i], i, WI[i]);
    */
    
    
desallocate_memory:

    if(WI)
    free(WI);
    
    if(WORK)
    free(WORK);
    
    return code;
} 
#else
{
    return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif


//calcute the eigenvalues of full symmetric matrix M. The eigenvalues will be in ascending ordar in eig vector.
int iquad::IQD_getEigenValues(const int n, const double *M, double* eig)
#if IQD_HAVE_LAPACK
{
    int aux, i, j, code;
    char JOBZ = 'N', UPLO = 'U'; //we set upper triangle because fortran storage by column. So, our lower triangle turns upper triangle...
    int LWORK = 20*n;
    double *A = NULL, *W = eig, *WORK=NULL;
    
    A = (double *) malloc(  n*n * sizeof(double)  ); //we need a nXn matrix...
    WORK=(double*) malloc( LWORK * sizeof(double) );
    if(!WORK)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    for(i = 0; i < n; i++)
    {
        aux = i*n;
        for(j = 0; j <= i; j++)
            A[aux + j] = M[aux + j];
    }
    
    aux = n;
    dsyev_(&JOBZ, &UPLO, &aux, A, &aux, W, WORK, &LWORK, &code);
    
desallocate_memory:
    
    if(A)	free(A);
    if(WORK)	free(WORK);
    
    return code;
    
}
#else
{
    return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif


//calcute the eigenvalues of symmetric sparse matrix. The eigenvalues will be in ascending order in eig vector.
int iquad::IQD_getEigenValues(const IQD_SparseMatrix& Q, double* eig)
#if IQD_HAVE_LAPACK
{
    int aux, code;
    //IQD_SparseElement *colAux;
    const int n = Q.getNumberOfRows(); // nrows;
    
    char JOBZ = 'N', UPLO = 'U'; //we set upper triangle because fortran storage by column. So, our lower triangle turns upper triangle...
    int LWORK = 20*n;
    double *A = NULL, *W = eig, *WORK=NULL;
    
    A = (double *) calloc(  n*n , sizeof(double)  ); //we need a nXn matrix...
    WORK=(double*) malloc( LWORK * sizeof(double) );
    if(!A || !WORK)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    /*for(i = 0; i < n; i++)
    {
        IQD_SparseRow &myrow = Q[i];
        
        aux = myrow.getNumberOfElements(); //Q.rows[i].nElements;
        //colAux = Q.rows[i].columns;
        
        double *a = &A[i*n];
        
        for(j = 0; j < aux; j++)
            a[ myrow[j].getColumn() ] = myrow[j].getValue(); //we store in A only the lower triangle...
    }*/
    
    
    Q.copyMatrixTo(A, false, false, 1.0);
    
    
    aux = n;
    dsyev_(&JOBZ, &UPLO, &aux, A, &aux, W, WORK, &LWORK, &code);
    
    //W and eig are the same and the eigenvalue are in ascending order.
    
    if( code != 0 )
        code = IQD_UNDEFINED_ERROR;
    
    
desallocate_memory:
    
    if(A)	free(A);
    if(WORK)	free(WORK);
    
    return code;
}
#else
{
    return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif



double iquad::IQD_getTime(void)
{
    return branchAndBound::BBL_getTime();
}
/*#ifdef IQD_HAVE_CLOCK_GETTIME
{
    timespec Time;
    
    clock_gettime(CLOCK_REALTIME, &Time);
    
    return Time.tv_sec + Time.tv_nsec/1.0e9;
    
    //time(NULL);
}
#else
{
    return (double) time(NULL);
}
#endif */



void iquad::IQD_welcome(void)
{
    static bool welcome = false;
    
    if(welcome == false)
    {
        std::cout << "----------------------------------------------------------------------------------\n"
        "iquad - Indefinite Quadratic Solver " IQD_VERSION "\n"
        "idealized by " IQD_IDEALIZER "\n"
        "developed by " IQD_AUTHOR ", " IQD_AUTHOR_FILIATION "\n"
        "collaborator: " IQD_COLLABORATORS "\n"
        "----------------------------------------------------------------------------------\n";
        
        welcome = true;
    }
}



IQD_OutFile::IQD_OutFile()
{
    outFile = NULL;
    lastHoursToSaveFile = 0.0;
}

IQD_OutFile::~IQD_OutFile()
{
    close();
}


int IQD_OutFile::open(const char* fileName)
{
    outFile = fopen(fileName, "a");
    
    if(!outFile)
        return IQD_UNDEFINED_ERROR;
    
    return 0;
}



int IQD_OutFile::writeToCompleteHours( const int ndata, const double maxCpuTime, const double hoursToSaveFile )
{
    const int maxdata = 10;
    const char charsep = IQD_CHAR_SEP[0];
    const double maxCpuTimeHours = maxCpuTime/(60*60);
    
    char text[maxdata + 1];
    int r = 0;
    double next;
    
    
    
    if( ndata > maxdata )
    {
        IQD_PRINTERRORMSG("ndata is greater than reserved space. Increase maxdata ");
        IQD_getchar();
        return IQD_BAD_DEFINITIONS;
    }
    
    for(int i = 0; i < ndata; i++)
        text[i] = charsep;
    text[ndata] = '\0';
    
    
    
    next = lastHoursToSaveFile + hoursToSaveFile;
    
    
    while( next <= maxCpuTimeHours )
    {
        r += fprintf(outFile, "%s", text);
        
        lastHoursToSaveFile = next;
        next += hoursToSaveFile;
    }
    
    
    return r;
}


int IQD_OutFile::writeProblemInfo(const char* name, IQD_IQuadProb &prob)
{
    int i, ml = 0, mq = 0, mnl = 0;
    int nzjac = 0, nzlagH = 0;
    
    
    for(i = 0; i < prob.m; i++)
    {
        if( prob.nlConstr[i] )
        {
            mnl++;
        }
        else
        {
            if(prob.QC[i].getNumberOfElements() == 0)
                ml++;
            else
                mq++;
        }
    }
    
    //if( prob.J )
        nzjac = prob.J.getNumberOfElements();
    
    //if(prob.lagH)
        nzlagH = prob.lagH.getNumberOfElements();
    
    return fprintf(outFile, "%s%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s%d%s", name, IQD_CHAR_SEP, prob.n, IQD_CHAR_SEP, prob.getNumberOfIntegerVars(), IQD_CHAR_SEP, ml, IQD_CHAR_SEP, mq, IQD_CHAR_SEP, prob.m, IQD_CHAR_SEP, prob.Q.getNumberOfElements(), IQD_CHAR_SEP, nzjac, IQD_CHAR_SEP, nzlagH, IQD_CHAR_SEP );
}


int IQD_OutFile::writeSolInfo(const int code, const double lb, const double objF, const double cpu_time, const double time, const unsigned long int niters, const double split_time, const double root_node_lb, const double root_node_ub)
{
    return fprintf(outFile, "%d%s%f%s%f%s%f%s%f%s%lu%s%f%s%f%s%f%s", code, IQD_CHAR_SEP, lb, IQD_CHAR_SEP, objF, IQD_CHAR_SEP, cpu_time, IQD_CHAR_SEP, time, IQD_CHAR_SEP, niters, IQD_CHAR_SEP, split_time, IQD_CHAR_SEP, root_node_lb, IQD_CHAR_SEP, root_node_ub, IQD_CHAR_SEP);
}

int IQD_OutFile::writeDouble(const double value)
{
    return fprintf(outFile, "%f%s", value, IQD_CHAR_SEP);
}

int IQD_OutFile::writeInt(const int value)
{
    return fprintf(outFile, "%d%s", value, IQD_CHAR_SEP);
}

int IQD_OutFile::writeString(const char *value)
{
    return fprintf(outFile, "%s%s", value, IQD_CHAR_SEP);
}

int IQD_OutFile::writeBreakLine()
{
    return fprintf(outFile, "\n");
}

int IQD_OutFile::flush()
{
    return fflush(outFile);
}

void IQD_OutFile::close()
{
    if(outFile)
        fclose(outFile);
    
    outFile = NULL;
}







int IQD_AuxBoundsNLEval::initialize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess)
{ 
    return rprob->nlEval->initialize( nthreads, rprob->n, rprob->m, rprob->J.getNumberOfElements(), rprob->lagH.getNumberOfElements() );
    
}


int IQD_AuxBoundsNLEval::eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
{
    value = 0.0;
    return 0;
}


int IQD_AuxBoundsNLEval::eval_nl_constrs_part( const int threadnumber, const int n, const int m, const bool newx, const bool *constrEval, const double *x, double *values)
{
    bool mynewx = newx;
    int code = 0, ret;
    
    
    if( rprob->hasNLConstraints() )
    {
        ret = rprob->nlEval->eval_nl_constrs_part( threadnumber, rprob->n, rprob->m, true, constrEval, x, values );
        
        if( ret != 0 )
            code = ret;
        
        mynewx = false;
    }
    
    
    if( constrEval[indobjcut] )
    {
        const double of = rprob->getObjFactor();
        
        ret = rprob->nlEval->eval_nl_obj_part( threadnumber, rprob->n, mynewx, x, values[indobjcut] );
        
        if( ret != 0 )
            code = ret;
        
        if( of != 1.0 )
            values[indobjcut] *= of;
    }
    
    
    return code;
}


int IQD_AuxBoundsNLEval::eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values)
{
    for(int i = 0; i < n; i++)
        values[i] = 0.0;
    
    return 0;
}


int IQD_AuxBoundsNLEval::eval_grad_nl_constrs_part( const int threadnumber, const int n, const int m, const int nz, const bool newx, const bool *constrEval, const double *x, MIP_SparseMatrix& jacobian)
{
    bool mynewx = newx;
    int code = 0, ret;
    
    
    if( rprob->hasNLConstraints() )
    {
        ret = rprob->nlEval->eval_grad_nl_constrs_part( threadnumber, rprob->n, rprob->m, rprob->J.getNumberOfElements(), true, constrEval, x, jacobian );
        
        if( ret != 0 )
            code = ret;
        
        mynewx = false;
    }
    
    
    if( constrEval[indobjcut] )
    {
        ret = rprob->nlEval->eval_grad_nl_obj_part( threadnumber, rprob->n, mynewx, x, jacobian(indobjcut) );
        
        if( ret != 0 )
            code = ret;
        
        //jacobian[indobjcut].setElementsByOrder( n, values );
    }
    
    
    return code;
}



//virtual int eval_hessian_nl_obj_part(const int threadnumber, const int n, const int nz, const bool newx, const double *x, MIP_SparseMatrix& hessian) = 0;

int IQD_AuxBoundsNLEval::eval_hessian_nl_lagran_part( const int threadnumber, const int n, const int m, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, MIP_SparseMatrix& hessian)
{
    const double myobjFactor = rprob->getObjFactor() * lambda[indobjcut];
    
    
    return rprob->nlEval->eval_hessian_nl_lagran_part( threadnumber, rprob->n, rprob->m, nz, newx, x, myobjFactor, lambda, hessian );
}


void IQD_AuxBoundsNLEval::finalize(const int nthreads, const int n, const int m, const int nzNLJac, const int nzNLLagHess)
{
    return rprob->nlEval->finalize( nthreads, rprob->n, rprob->m, rprob->J.getNumberOfElements(), rprob->lagH.getNumberOfElements() );
}


IQD_AuxBoundsNLEval::~IQD_AuxBoundsNLEval()
{
    
}






IQD_AuxBoundsUpdate::IQD_AuxBoundsUpdate()
{
    initialize();
}



IQD_AuxBoundsUpdate::~IQD_AuxBoundsUpdate()
{
    desallocateAuxArrays();
    desallocateSolver();
}



int IQD_AuxBoundsUpdate::allocateAuxArrays( const int size )
{
    cols = (int *) malloc( size * sizeof(int) );
    vals = (double *) malloc( size * sizeof(double) );
    
    if( !cols || !vals )
        return IQD_MEMORY_ERROR;
    
    return 0;
}



int IQD_AuxBoundsUpdate::buildPreProblem( IQD_IQuadProb &rprob, const int dimv, const int qpSolver, const int qcpSolver, const int nlpSolver)
{
    const int n = rprob.getNumberOfVars();
    bool nlpProb;
    int code, ret;
    
    int *myrows = NULL, *mycols = NULL;
    double *myvals = NULL;
    
    
    ret = rprob.getProblemType();
    
    nlpProb = MIP_isNonlinearProblemType(ret);
    
    
    if( nlpProb )
        solver = OPT_newNLPSolver(nlpSolver);
    else //if( MIP_isQuadConstrProblemType(ret) )
        solver = OPT_newQCPSolver(qcpSolver); //even if rprob is only quadratiuc, we need a qcp problem becuase the objective cut is a quadratic constraint...
    
    
    if( !solver )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    //we set the problem without the objective function...
    ret = optsolvers::OPT_setMINLPProblem(rprob, solver, false, true, true, false);
    if( ret != 0 )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTERRORNUMBER(ret);
        #endif
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    indobjcut = rprob.getNumberOfConstraints();
    indsecconstrs = indobjcut + 1;
    
    //we already allocate space for objective cut and secant equalities also...
    ret = solver->addConstraints(1 + dimv);
    if( ret != 0 )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTERRORNUMBER(ret);
        #endif
        code = IQD_SOLVER_ERROR;
        goto termination;
    }
    
    
    
    //adding objective cut
    {
        int nzs;
        OPT_QCPSolver *qcp = (OPT_QCPSolver *) solver;
        
        const int nzq =  rprob.getNumberOfObjQuadTerms();
        const double fac = rprob.getObjFactor();
        double oc = rprob.getObjConstant();
        
        const double *c = rprob.c;
        
        const int size = IQD_max(n, nzq);
        
        
        if( nzq > 0 )
        {
            myrows = (int*) malloc(size* sizeof(int) );
            mycols = (int*) malloc(size* sizeof(int) );
            myvals = (double *) malloc( size* sizeof(double) );
            
            if( !myrows || !mycols || !myvals )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTMEMERROR;
                #endif
                code = IQD_MEMORY_ERROR;
                goto termination;
            }
            
            rprob.getObjQuadCoefsMatrix(&nzs, myrows, mycols, myvals);
        }
        
        
        
        nzs = 0;
        for(int i = 0; i < n; i++)
        {
            if( c[i] != 0.0 )
            {
                cols[nzs] = i;
                vals[nzs] = c[i];
                nzs++;
            }
        }
        
        
        if( fac != 1.0 )
        {
            for(int i = 0; i < nzq; i++)
                myvals[i] *= fac;
            
            for(int i = 0; i < nzs; i++)
                vals[i] *= fac;
            
            oc *= fac;
        }
        
        
        ret = qcp->resetConstraintLinearPart( indobjcut, nzs, cols, vals);
        if( ret != 0 )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTERRORNUMBER(ret);
            #endif
            code = IQD_SOLVER_ERROR;
            goto termination;
        }
        
        
        ret = qcp->setConstraintQuadMatrix( indobjcut, nzq, myrows, mycols, myvals);
        if( ret != 0 )
        {
            #if IQD_DEBUG_MODE
                IQD_PRINTERRORNUMBER(ret);
            #endif
            code = IQD_SOLVER_ERROR;
            goto termination;
        }
        
        //we do not set the rhs becuase we do not have the zu value... :D
        
        
        if( nlpProb )
        {
            if( rprob.hasObjNLTerm() )
            {
                OPT_NLPSolver *nlp = (OPT_NLPSolver*) solver;
                
                nleval = new (std::nothrow) IQD_AuxBoundsNLEval;
                if( !nleval )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTMEMERROR;
                    #endif
                    code = OPT_MEMORY_ERROR;
                    goto termination;
                }
                
                nlp->setNonLinearEvalObject( nleval );
                
                
                ret = nlp->setConstrNLFlag( indobjcut, true);
                if( ret != 0 )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTERRORNUMBER(ret);
                    #endif
                    code = IQD_SOLVER_ERROR;
                    goto termination;
                }
                
                
                for(int i = 0; i < n; i++)
                    cols[i] = i;
                
                ret = nlp->setJacobianRowStructure( indobjcut, n, cols );
                
                if( ret != 0 )
                {
                    #if IQD_DEBUG_MODE
                        IQD_PRINTERRORNUMBER(ret);
                    #endif
                    code = IQD_SOLVER_ERROR;
                    goto termination;
                }
            }
            
        }
        
    }
    
    
    code = 0;
    
termination:
    
    if( code != 0 )
        desallocateSolver();
    
    return code;
}



int IQD_AuxBoundsUpdate::calculateNewBounds( const int norigvars, const int dimv, const int indvconstrs, IQD_SparseMatrix &v, const double zu, const double *ly, const double *uy, double *newly, double *newuy )
{
    OPT_LPSolver *lp = (OPT_LPSolver *) solver;
    int ret, code = 0;
    
    
    IQD_setAllarray<double>(dimv, newly, -INFINITY);
    IQD_setAllarray<double>(dimv, newuy, INFINITY);
    
    
    //updating rhs in objective cut...
    ret = lp->setConstraintBounds(indobjcut, -OPT_INFINITY, zu);
    if(ret != 0)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTERRORNUMBER(ret);
        #endif
        code = ret;
        goto termination;
    }
    
    
    //updating bounds from auxiliary variables
    for(int i = 0, k = 0; i < dimv; i++)
    {
        //IQD_SparseRow &vrow = v[i];
        
        const int *rcols = v.getRowColsPointer(i);
        const double *rvals = v.getRowValuesPointer(i);
        
        
        const int nel = v.getNumberOfElementsAtRow(i); //vrow.getNumberOfElements();
        
        if( nel == 1 )
        {
            //we change variable bounds directly instead of set a constraint...
            const double d = rvals[0]; // vrow[0].getValue();
            const int ind = rcols[0]; // vrow[0].getColumn();
            
            double lb = ly[i];
            double ub = uy[i];
            
            if( d != 1.0 )
            {
                lb = lb/d;
                ub = ub/d;
                
                if( d < 0 )
                    IQD_swap(lb, ub);
            }
            
            ret = lp->setVariableBounds(ind, lb, ub);
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                code = IQD_SOLVER_ERROR;
                goto termination;
            }
            
        }
        else if( nel > 1 )
        {
            ret = lp->setConstraintBounds(indvconstrs + k, ly[i], uy[i]);
            if( ret != 0 )
            {
                #if IQD_DEBUG_MODE
                    IQD_PRINTERRORNUMBER(ret);
                #endif
                code = IQD_SOLVER_ERROR;
                goto termination;
            }
            
            k++;
        }
    }
    
    
    //setting secants like equalities instead of inequalities...
    ret = IQD_setSecantConstraints(norigvars, dimv, v, indsecconstrs, lp, ly, uy, true, cols, vals);
    
    if(ret != 0)
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTERRORNUMBER(ret);
        #endif
        code = ret;
        goto termination;
    }
    
    
    
    for(int i = 0; i < dimv; i++)
    {
        const int nzs = v.getNumberOfElementsAtRow(i); //v[i].getStructureAndValues(cols, vals);
        
        const int *rcols = v.getRowColsPointer(i);
        const double *rvalues = v.getRowValuesPointer(i);
        
        ret = lp->setObjLinearCoefs(nzs, rcols, rvalues); //lp->setLinearObjCoefs(nzs, cols, vals);
        
        if( ret == 0 )
        {
            ret = lp->setObjSense( OPT_MINIMIZE );
            if( ret == 0 )
            {
                ret = lp->solve(false);
                
                //cout << "ret minimize: " << ret << endl;
                
                if( ret == OPT_OPTIMAL_SOLUTION )
                {
                    newly[i] = lp->objValue;
                }
                else
                    code = IQD_SOLVER_ERROR;
            }
            else
                code = IQD_SOLVER_ERROR;
            
            
            ret = lp->setObjSense( OPT_MAXIMIZE );
            if( ret == 0 )
            {
                ret = lp->solve(false);
                
                if( ret == OPT_OPTIMAL_SOLUTION )
                {
                    newuy[i] = lp->objValue;
                }
                else
                    code = IQD_SOLVER_ERROR;
            }
            else
                code = IQD_SOLVER_ERROR;
        }
        else
        {
            code = IQD_SOLVER_ERROR;
        }
        
        
        IQD_setAllarray(nzs, vals, 0.0);
        ret = lp->setObjLinearCoefs(nzs, rcols, vals);
        
        if( ret != 0 )
        {
            //unfortunatelly, we cannot continue...
            code = IQD_SOLVER_ERROR;
            break;
        }
        
    }
    
    
termination:
    
    return code;
}



void IQD_AuxBoundsUpdate::desallocateAuxArrays()
{
    IQD_secFree(cols);
    IQD_secFree(vals);
}


void IQD_AuxBoundsUpdate::desallocateSolver()
{
    IQD_secDelete(solver);
    IQD_secDelete(nleval);
}


void IQD_AuxBoundsUpdate::initialize()
{
    #if IQD_DEBUG_MODE
        indobjcut = INT_MAX;
        indsecconstrs = INT_MAX;
    #endif
    
    cols = NULL;
    vals = NULL;
    
    nleval = NULL;
    solver = NULL;
}


















void iquad::IQD_generateRstatistic(const char* probName, IQD_IQuadProb& prob, const IQD_SPLITQUAD splitWay, const int sdpSolver, const int milpSolver, const int qcpSolver, const int nlpSolver, const int minlpSolver, IQD_GeneralSolverParams* sParams, IQD_OutFile& outFile, const double zero_tol)
{
    /*const int n = prob.n;
    IQD_RefProb refProb;
    
    refProb.prob = &prob;
    refProb.setWayOfSplitQ(splitWay);
    
    
    refProb.splitQ(sdpSolver, milpSolver, qcpSolver, nlpSolver, milpSolver, sParams, NULL, NULL, false);
    
    if( IQD_isDiagonalSplitQWay( splitWay ) )
    {
        refProb.allocateR(n);
        
        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
                refProb.R[n*i + j] = 0.0;
            
            refProb.R[n*i + i] = refProb.lambda[i];
        }
    }
    
    IQD_writeREigInfo(probName, n, refProb.R, outFile, zero_tol); */
}


int iquad::IQD_writeREigInfo(const char *probName, const int n, const double* R, IQD_OutFile &outFile, const double zero_tol)
{
    int i, code, nNegEig = 0;
    double *eig = NULL;
    double trace = 0.0, mean = 0.0, median;
    
    eig = (double *) malloc( n * sizeof(double) );
    if(!eig)
    {
        code = IQD_UNDEFINED_ERROR;
        goto desallocate_memory;
    }
    
    i = IQD_getEigenValues(n, R, eig);
    for(i = 0; i < n; i++)
        mean += eig[i];
    
    mean = mean/n;
    
    for(i = 0; i < n; i++)
        trace += R[i*n + i];
    
    for(i = 0; i < n; i++)
        if(eig[i] < -zero_tol)
            nNegEig++;
    
    //getting median
    if(n % 2 == 0)
    {
        i = n / 2;
        median = (eig[i-1] + eig[i])/2.0; 
    }
    else
        median = eig[ (n-1)/2 ];
    
    //outFile.writeString(probName);
    outFile.writeInt(nNegEig);
    outFile.writeDouble( (double) nNegEig/n );
    outFile.writeDouble( eig[0] );
    outFile.writeDouble( eig[n-1] );
    outFile.writeDouble( trace );
    outFile.writeDouble( mean );
    outFile.writeDouble( median );
    
    //fprintf(file, "%d%s%f%s%f%s%f%s%f%s%f%s%f%s", nNegEig, IQD_CHAR_SEP, (double) nNegEig/n, IQD_CHAR_SEP, eig[0], IQD_CHAR_SEP, eig[n-1], IQD_CHAR_SEP, trace, IQD_CHAR_SEP, mean, IQD_CHAR_SEP, median, IQD_CHAR_SEP );
    
desallocate_memory:
    
    //if(file)	fclose(file);
    if(eig)	free(eig);
    
    return code;
}



int iquad::IQD_writeEigInfo(const char *probName, const IQD_SparseMatrix &Q, const char *fileName, const double zero_tol)
{
    const int n = Q.getNumberOfRows();
    int i, code, nNegEig = 0;
    double *eig = NULL;
    double mean = 0.0, median;
    FILE *file = NULL;
    
    eig = (double *) malloc( n * sizeof(double) );
    file = fopen(fileName, "a");
    if(!file)
    {
        code = IQD_UNDEFINED_ERROR;
        goto desallocate_memory;
    }
    
    i = IQD_getEigenValues(Q, eig);
    
    for(i = 0; i < n; i++)
    {
        mean += eig[i];
        
        if(eig[i] < -zero_tol)
        {
            nNegEig++;
        }
        
        printf("%f ", eig[i]);
    }
    
    printf("\n");
    
    
    mean = mean/n;
    
    //getting median
    if(n % 2 == 0)
    {
        i = n / 2;
        median = (eig[i-1] + eig[i])/2.0; 
    }
    else
        median = eig[ (n-1)/2 ];
    
    
    fprintf(file, "\n%s%s%d%s%d%s%f%s%f%s%f%s%f%s%f%s", probName, IQD_CHAR_SEP, n, IQD_CHAR_SEP, nNegEig, IQD_CHAR_SEP, (double) nNegEig/n, IQD_CHAR_SEP, eig[0], IQD_CHAR_SEP, eig[n-1], IQD_CHAR_SEP, mean, IQD_CHAR_SEP, median, IQD_CHAR_SEP );
    
    
    
    
desallocate_memory:
    
    if(file)	fclose(file);
    if(eig)	free(eig);
    
    return code;
    
}



double iquad::IQD_maxDifferenceBetweenMatrices(const int nrows, const int ncols, const double* M1, const double* M2)
{
    int i;
    double dif = 0.0;
    
    for(i = 0; i < nrows*ncols; i++)
    dif = IQD_max(  dif, IQD_abs(M1[i] - M2[i]) );
    
    return dif;
}






//that function invert negative eigenvalues except the lowest in matrix Q.
int iquad::IQD_inverteQeigenvalues( IQD_IQuadProb &prob, const double zero_tol)
{
    const int n = prob.n;
    int i, j, k, code;
    double smallEig;
    double *Q = NULL, *M = NULL, *Zt = NULL, *eig = NULL;
    
    prob.Q.copyMatrixTo(Q, true, true, 1.0); //getFullMatrixCopy();
    M= (double *) malloc( n*n*sizeof(double) );
    Zt= (double *) malloc( n*n*sizeof(double) );
    eig = (double *) malloc( n * sizeof(double) );
    if( !Q || !M || !Zt || !eig )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    IQD_schurDecomposition(Q, n, Q, Zt, eig);
    
    
    printf("eig:\n");
    for(i = 0; i < n; i++)
        printf("%0.16f; ", eig[i]);
    printf("\n");
    
    /*printf("T\n");
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            printf("%0.16f ", Q[i*n + j] );
        printf(";\n");
    }
    
    
    printf("Zt\n");
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            printf("%0.16f ", Zt[i*n + j] );
        printf(";\n");
    } */
    getchar();
    
    
    //Q is a symmetric matrix. So, after schur decomposition, we have Q as a diagonal matrix with eigenvalues.
    

    smallEig = 0.0;
    for(j = 0; j < n; j++)
    {
        if( eig[j] < smallEig )
        {
            smallEig = eig[j];
            k = j;
        }
    }
        
    if( smallEig != 0.0)
    {
        for(j = 0; j < n; j++)
        {
            //inverting the eigenvalues
            if( eig[j] < -zero_tol && k != j )
            eig[j] = -eig[j];
        }
    }
    
    
        //recomposing the new Q matrix Schur factorization A = Z*T*(Z**T).
        
    //M = Z*T  Z is transposed in Zt and T is a diagonal matrix with eigenvalues...W
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            M[i*n + j] = Zt[j*n + i] * eig[j];    
        }
    }

    //Q = M*Zt
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            Q[i*n + j] = 0.0;
            for(k = 0; k < n; k++)
                Q[i*n + j] += M[i*n + k] * Zt[k*n + j];
        }   
    }
    
    
    prob.Q.deleteStructure();
    
    //using M to store the lower triangle of new matrix Q
    
    for(i = 0; i < prob.n; i++)
    {
        for(j = 0; j <= i; j++)
            M[ (i*(i+1))/2 + j ] = Q[i*n + j];
    }
    
    prob.Q.setSymmetricStructureAndValuesByLowerTriangle(M, zero_tol);
    
    
    code = 0;
    
desallocate_memory:
    
    if(Q)	free(Q);
    if(M)	free(M);
    if(Zt)	free(Zt);
    if(eig)	free(eig);
    
    return code;
}


//that function invert eigenvalues on matrix Q until it has nNegEig. The highest eigenvalue in module are keep (or turned) negative
int iquad::IQD_changeQeigenvalues( IQD_IQuadProb &prob, const int nNegEig, const double zero_tol)
{
    bool *flags;
    const int n = prob.n;
    int i, j, k, code;
    double aux, highestEig;
    double *Q = NULL, *M = NULL, *Zt = NULL, *eig = NULL;
    
    
    prob.Q.copyMatrixTo(Q, true, true, 1.0); //getFullMatrixCopy();
    M= (double *) malloc( n*n*sizeof(double) );
    Zt= (double *) malloc( n*n*sizeof(double) );
    eig = (double *) malloc( n * sizeof(double) );
    flags= (bool *) calloc( n , sizeof(bool) );
    if( !Q || !M || !Zt || !eig || !flags )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    IQD_schurDecomposition(Q, n, Q, Zt, eig);
    
    
    for(i = 0; i < nNegEig; i++)
    {
        highestEig = 0.0;
        
        for(j = 0; j < n; j++)
        {
            aux = IQD_abs(eig[j]);
            
            if( aux > highestEig && flags[j] == false )
            {
                highestEig = aux;
                k = j;
            }
        }
        
        if( highestEig == 0.0 )
        {
            code = IQD_BAD_DEFINITIONS;
            goto desallocate_memory;
        }
        
        flags[k] = true;
    }
    
    
    for(i = 0; i < n; i++)
    {
        if( flags[i] )
            eig[i] = -IQD_abs(eig[i]); //that eigenvalue should be negative
        else
            eig[i] =  IQD_abs(eig[i]); //that eigenvalue should be negative
    }
    
    //recomposing the new Q matrix Schur factorization A = Z*T*(Z**T).
        
    //M = Z*T  Z is transposed in Zt and T is a diagonal matrix with eigenvalues...W
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            M[i*n + j] = Zt[j*n + i] * eig[j];    
        }
    }

    //Q = M*Zt
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            Q[i*n + j] = 0.0;
            for(k = 0; k < n; k++)
                Q[i*n + j] += M[i*n + k] * Zt[k*n + j];
        }   
    }
    
    
    prob.Q.deleteStructure();
    
    //using M to store the lower triangle of new matrix Q
    
    for(i = 0; i < prob.n; i++)
    {
        for(j = 0; j <= i; j++)
            M[ (i*(i+1))/2 + j ] = Q[i*n + j];
    }
    
    prob.Q.setSymmetricStructureAndValuesByLowerTriangle(M, zero_tol);
    
    
    
    
    code = 0;
    
desallocate_memory:
    
    return code;
}



//that function replace zero eingevalue by random values beetwen 
int iquad::IQD_replaceZeroEigValuesOfQByRandom(IQD_IQuadProb &prob, const double zero_tol)
{
    const int n = prob.n;
    int i,j, k, code;
    double *Q = NULL, *M = NULL, *Zt = NULL, *eig = NULL;
    double maxeig, mineig;
    IQD_Random random;
    
    
    prob.Q.copyMatrixTo(Q, true, true, 1.0); // getFullMatrixCopy();
    M= (double *) malloc( n*n*sizeof(double) );
    Zt= (double *) malloc( n*n*sizeof(double) );
    eig = (double *) malloc( n * sizeof(double) );
    if( !Q || !M || !Zt || !eig )
    {
        #if IQD_DEBUG_MODE
            IQD_PRINTMEMERROR;
        #endif
        code = IQD_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    IQD_schurDecomposition(Q, n, Q, Zt, eig);
    
    //Q is a symmetric matrix. So, after schur decomposition, we have Q as a diagonal matrix with eigenvalues...
    
    
    maxeig = mineig = eig[0];
    
    for(i = 1; i < n; i++)
    {
        if( eig[i] > maxeig )
            maxeig = eig[i];
        
        if( eig[i] < mineig )
            mineig = eig[i];
    }
    
    random.setSeed(NULL);
    
    for(i = 0; i < n; i++ )
    {
        while( IQD_abs( eig[i] ) <= zero_tol  )
        {
            eig[i] = random.randInt(mineig, maxeig);
        }
    }
    
    
    //M = Z*T  Z is transposed in Zt and T is a diagonal matrix with eigenvalues...W
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            M[i*n + j] = Zt[j*n + i] * eig[j];    
        }
    }

    //Q = M*Zt
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
        {
            Q[i*n + j] = 0.0;
            for(k = 0; k < n; k++)
                Q[i*n + j] += M[i*n + k] * Zt[k*n + j];
        }   
    }
    
    
    prob.Q.deleteStructure();
    
    //using M to store the lower triangle of new matrix Q
    
    for(i = 0; i < prob.n; i++)
    {
        for(j = 0; j <= i; j++)
            M[ (i*(i+1))/2 + j ] = Q[i*n + j];
    }
    
    prob.Q.setSymmetricStructureAndValuesByLowerTriangle(M, zero_tol);
    
    
    
    code = 0;
    
desallocate_memory:
    
    if(Q)	free(Q);
    if(M)	free(M);
    if(Zt)	free(Zt);
    if(eig)	free(eig);
    
    return code;
}




int iquad::IQD_solveLinearSystem(const int n, const double *A, double *b, double *x)
{
    int i, j, INFO, code, iter, *IPIV = NULL;
    float *swork = NULL;
    double *newA = NULL, *work = NULL;
    
    
    newA = (double *) malloc( n*n * sizeof(double) );
    work = (double *) malloc( n * sizeof(double) );
    swork= (float *) malloc( n *(n +1)*sizeof(float) );
    IPIV = (int *) malloc( n * sizeof(int) );
    if(!newA || !work || !swork || !IPIV)
    {
        code = IQD_MEMORY_ERROR;
        goto desallocate_memory;
    }
    
    
    //fortran store matrix by columns
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            newA[j*n + i]  = A[i*n + j];
    }
    
    
    printf("A:\n");
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < n; j++)
            printf("%f ", A[i*n + j]);
        printf(";\n");
    }
    
    
    //for(i = 0; i < n; i++)
        //x[i] = b[i];
    
    i = n;
    j = 1;
    
    dsgesv_(&i, &j, newA, &i, IPIV, b, &i, x, &i, work, swork, &iter, &INFO);
    
    
    printf("dgesv_ INFO: %d\n", INFO);
    
    printf("x:\n");
    for(i = 0; i < n; i++)
        printf("%f ", x[i]);
    
    
    if(INFO == 0)
        code = 0;
    else if(INFO < 0)
        code = IQD_BAD_DEFINITIONS;
    else
        code = IQD_UNDEFINED_ERROR;
    
    
    
    
    
desallocate_memory:

    if(newA) free(newA);
    if(IPIV) free(IPIV);
    
    return code;
}




double iquad::IQD_norm2(const int n, const double *v)
{
    /*int i;
    double norm = 0.0;
    
    for(i = 0; i < n; i++)
        norm += v[i]*v[i]; */
    
    return sqrt( IQD_vectorTimes(n, v, v) );
}





int iquad::IQD_getStrEnumValue(const char *svalue, IQD_SDP_SOLVERS &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_SS_CSDP), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_SS_MOSEK), svalue, value ) == 0)
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}


int iquad::IQD_getStrEnumValue(const char *svalue, IQD_QP_SOLVERS &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_QP_CPLEX), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QP_GUROBI), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QP_XPRESS), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QP_MOSEK), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QP_IPOPT), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QP_WORHP), svalue, value ) == 0)
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}


int iquad::IQD_getStrEnumValue(const char *svalue, IQD_QCP_SOLVERS &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_QCP_XPRESS), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QCP_MOSEK), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QCP_IPOPT), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_QCP_WORHP), svalue, value ) == 0)
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}

int iquad::IQD_getStrEnumValue(const char *svalue, IQD_MILP_SOLVERS &value)
{
    int ret = 0;
    
    
    if( IQD_setEnum( IQD_STRATT(IQD_IP_GLPK), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_IP_CPLEX), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_IP_GUROBI), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_IP_XPRESS), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_IP_MOSEK), svalue, value ) == 0)
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}


int iquad::IQD_getStrEnumValue(const char *svalue, IQD_NLP_SOLVERS &value)
{
    int ret = 0;
    
    
    if( IQD_setEnum( IQD_STRATT(IQD_NLP_MOSEK), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_NLP_IPOPT), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_NLP_WORHP), svalue, value ) == 0)
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}


int iquad::IQD_getStrEnumValue(const char *svalue, IQD_MINLP_SOLVERS &value)
{
    int ret = 0;
    
    
    if( IQD_setEnum( IQD_STRATT(IQD_MS_MURIQUI), svalue, value ) == 0)
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_MS_BONMIN), svalue, value ) == 0)
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}


int iquad::IQD_getStrEnumValue(const char *svalue, IQD_SPLITQUAD &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_SQ_SCHUR), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_SQ_DIAG_SDP), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_SQ_DIAG_DOM), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_SQ_IDENTITY), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_SQ_MIX_SCHUR_DIAG_SDP), svalue, value ) == 0 )
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}



int iquad::IQD_getStrEnumValue(const char *svalue, IQD_INTSTRATEGY &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_IS_ON_BB_ENUMERATION), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_IS_ON_EACH_SUBPROBLEM), svalue, value ) == 0 )
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}



int iquad::IQD_getStrEnumValue(const char *svalue, IQD_EXP_STRATEGY &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_ES_DEPTH), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_ES_WIDTH), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_ES_BEST_LIMIT), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_ES_DEPTH_BEST_LIMIT), svalue, value ) == 0 )
    {}
        ret = IQD_VALUE_ERROR;
    
    return ret;
}


int iquad::IQD_getStrEnumValue(const char *svalue, IQD_BRANCH_STRATEGY &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_BS_MAX_CONSTRS_VIOLATED), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_BS_MOST_VIOLATED_CONSTR), svalue, value ) == 0 )
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    return ret;
}



int iquad::IQD_getStrEnumValue(const char *svalue, IQD_AUX_BOUNDS_CALC_STRATEGY &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_ABCS_SUBPROBLEMS), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_ABCS_ORIG_BOUNDS), svalue, value ) == 0 )
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    
    return ret;
}



int iquad::IQD_getStrEnumValue(const char *svalue, IDQ_RELAX_PROB_WRITING_MODE &value)
{
    int ret = 0;
    
    if( IQD_setEnum( IQD_STRATT(IQD_RPWM_NO_WRITING), svalue, value ) == 0 )
    {}
    else if( IQD_setEnum( IQD_STRATT(IQD_RPWM_WRITING), svalue, value ) == 0 )
    {}
    else
        ret = IQD_VALUE_ERROR;
    
    
    return ret;
}



int iquad::IQD_setSecantConstraints( const int n, const int dimv, const IQD_SparseMatrix& v, const int indsecconstrs, OPT_LPSolver* solver, const double* ly, const double* uy, const bool equality, int* auxColumns, double* auxValues )
{
    //const int dimv = refProb2->dimv;
    //const int n = oprob->n;
    
    //IQD_SparseMatrix &v = refProb2->v;
    
    int code = 0, r, nzs;
    double lb = -OPT_INFINITY, ub, aux;
    //int *auxColumns;// = auxColumns;
    //double *auxValues;// = auxValues;
    
    
    
    for(int i = 0; i < dimv; i++)
    {
        
        if( ly[i] == uy[i] ) 
        {
            //we will fix w_i variable and we need no secant...
            
            r = solver->resetConstraintLinearPart(indsecconstrs+i, 0, NULL, NULL);
            
            r += solver->setConstraintBounds( indsecconstrs+i, -OPT_INFINITY, OPT_INFINITY );
            
            #if IQD_DEBUG_MODE
                assert(r == 0);
            #endif
            
            continue;
        }
        
        aux = -( uy[i] + ly[i] );
        
        ub = -uy[i]*ly[i];
        if( equality )
            lb = ub;
        
        
        //IQD_SparseRow &rowv = v[i];
        //const int *rcols = v.getRowColsPointer(i);
        const double *rvalues = v.getRowValuesPointer(i);
        
        //nzs = v.getNumberOfElementsAtRow(i); //rowv.getNumberOfElements();
        
        v.getRowStructure(i, auxColumns, &nzs); // IQD_copyArray(nzs, rcols, auxColumns);
        
        #pragma ivdep
        #pragma GCC ivdep
        for( int j = 0; j < nzs; j++ )
        {
            //auxColumns[j] = rcols[j]; //rowv[j].getColumn();
            auxValues[j] = aux *rvalues[j]; //aux *rowv[j].getValue() ;
            
            //printf("i: %d auxColumn[%d]: %d rvalues[%d]: %f auxValues[%d]: %f boxsize: %f\n", i, j, auxColumns[j], j, rvalues[j], j, auxValues[j], aux);
        }
        
        /*nzs = v[i].getStructureAndValues( cols[1], vals[1] );
        
        for(int j = 1; j < nzs; j++ )
            vals[j] *= aux; */
        
        
        auxColumns[nzs] = n + i; //w_i
        auxValues[nzs] = -1.0;  //w_i
        
        
        r = solver->resetConstraintLinearPart( indsecconstrs +i, nzs + 1, auxColumns, auxValues ) + 
        solver->setConstraintBounds( indsecconstrs +i, lb, ub);
        
        if( r != 0 )
        {
            IQD_PRINTERRORNUMBER(r);
            code = IQD_SOLVER_ERROR;
        }
        
    }
    
    
    return code;
}



int iquad::IQD_isConvexProblemAboutQuadratics( const IQD_IQuadProb& prob, const double zero_tol, bool& convexObj, bool& convexConstr )
{
    const int n = prob.n;
    const int m = prob.m;
    
    int code = 0, r;
    double *eig = NULL;
    
    const IQD_SparseMatrix &Q = prob.Q;
    const IQD_SparseMatrix *QC = prob.QC;
    
    
    
    convexConstr = false;
    convexObj = false;
    
    
    eig = (double *) malloc( n * sizeof(double) );
    if( !eig )
    {
        code = IQD_MEMORY_ERROR;
        goto termination;
    }
    
    
    convexConstr = true;
    
    //we check first the constraints...
    for(int i = 0; i < m; i++)
    {
        if( QC[i].getNumberOfElements() == 0 )
            continue;
        
        double lci, uci;
        
        r = prob.getConstraintBounds(i, lci, uci);
        #if IQD_DEBUG_MODE
            if( r != 0 )
            {
                IQD_PRINTERRORNUMBER(r);
                code = IQD_UNDEFINED_ERROR;
                goto termination;
            }
        #endif
        
        if( lci <= -MIP_INFINITY && uci >= MIP_INFINITY )
        {
            //free constraint. We do not check
            continue;
        }
        
        
        if( lci > -MIP_INFINITY && uci < MIP_INFINITY )
        {
            //we have a double bounded quadratic constraint. By definiton, it is nonconvex...
            convexConstr = false;
            break;
        }
        
        
        r = IQD_getEigenValues(QC[i], eig);
        if( r != 0 )
        {
            convexConstr = false;
            code = r;
            break;
        }
        
        
        if( uci < MIP_INFINITY )
        {
            //matrix should be Positive SD
            for(int i = 0; i < n; i++)
            {
                if( eig[i] < -zero_tol )
                {
                    convexConstr = false;
                    break;
                }
            }
        }
        else //if (lci > -MIP_INFINITY)
        {
            //matrix should be Negative SD
            for(int i = 0; i < n; i++)
            {
                if( eig[i] > zero_tol )
                {
                    convexConstr = false;
                    break;
                }
            }
        }
        
    }
    
    
    convexObj = true;
    
    if(Q.getNumberOfElements() > 0)
    {
        r = IQD_getEigenValues(Q, eig);
        if( r != 0 )
        {
            convexObj = false;
            code = r;
            goto termination;
        }
        
        for( int i = 0; i < n; i++ )
        {
            if( eig[i] < -zero_tol )
            {
                convexObj = false;
                goto termination;
            }
        }
    }
    
    
    
    
    
termination:
    
    if(eig)		free(eig);
    
    return code;
}






int iquad::IQD_scaleProblemAndTurnIntegerVars( IQD_IQuadProb &prob, const double scaleFactor )
{
    const double revFactor2 = 1.0/( scaleFactor*scaleFactor);
    const int n = prob.n;
    const int m = prob.m;
    
    double *lx = prob.lx, *ux = prob.ux;
    IQD_SparseMatrix *QC = prob.QC;
    
    
    
    //that is not ideal, but we modify MIP_MINLPProblem
    
    if( prob.hasLinCoefObj() )
    {
        double *c = prob.c;
        
        #pragma GCC ivdep
        #pragma ivdep
        for(int i = 0; i < n; i++)
            c[i] /= scaleFactor;
    }
    
    
    prob.Q.multiplyAllElements( revFactor2 ); //maybe would be better run all elements and divide by scaleFactor^2, but in this way is so easy...
    
    
    prob.A.multiplyAllElements( 1.0/scaleFactor );
    
    for( int i = 0; i < m; i++ )
        QC[i].multiplyAllElements( revFactor2 );
    
    
    #pragma GCC ivdep
    #pragma ivdep
    for( int i = 0; i < n; i++ )
        lx[i] *= scaleFactor;
    
    
    #pragma GCC ivdep
    #pragma ivdep
    for( int i = 0; i < n; i++ )
        ux[i] *= scaleFactor;
    
    
    #pragma GCC ivdep
    #pragma ivdep
    for( int i = 0; i < n/2; i++ )
        prob.setVariableType(i, minlpproblem::MIP_VT_INTEGER);
    
    
    
    
    return 0;
}








