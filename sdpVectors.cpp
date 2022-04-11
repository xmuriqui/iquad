/*
 * Implementation of a class that receives a a set of vectors v_i and calculates lambda solving the SDP problem:
 * 
 * Minimize \lambda
 * s.t.:
 * 
 * 		Q + \sum{  \lambda_i V_i   }
 * 
 * 		where V_i = v_i * v_i'
 * 
 */


#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <list>


#include "iquad.hpp"
#include "IQD_tools.hpp"

#include "OPT_solvers.hpp"


#if IQD_HAVE_CSDP
    
    extern "C" {
	#include "declarations.h" // Include CSDP declarations so that we'll know the calling interfaces. 
    }
#endif

#if IQD_HAVE_MOSEK
    #include "mosek.h"
#endif


using namespace iquad;
using namespace std;


#if IQD_HAVE_GUROBI
	//implemented on gurobi.cpp
	int IQD_setGurobiParameters(GRBenv *grbEnv, IQD_GeneralSolverParams *params);
#endif


#if IQD_HAVE_MOSEK
	//implemented on miqpMosek.cpp
	/*void IQD_setMosekParams(MSKtask_t &task, IQD_GeneralSolverParams *params)
	{
		MSKrescodee r;
		list< optsolvers::OPT_SolverParam<int>* >::iterator itInt;
		list< optsolvers::OPT_SolverParam<double>* >::iterator itDbl;
		list< optsolvers::OPT_SolverParam<char*>* >::iterator itStr;
		
		
		if(!params)
			return;
		
		
		for( itInt = params->intParams.begin(); itInt != params->intParams.end(); itInt++ )
		{
			r = MSK_putnaintparam(task, (*itInt)->name,  (*itInt)->value);
			if ( r != MSK_RES_OK )
			{
				fprintf(stderr, "Error %d at setting mosek integer parameter: %s at the value: %d!\n", r, (*itInt)->name, (*itInt)->value );
				//getchar();
			}
		}
		
		for( itDbl = params->dblParams.begin(); itDbl != params->dblParams.end(); itDbl++ )
		{
			r = MSK_putnadouparam(task, (*itDbl)->name,  (*itDbl)->value);
			if ( r != MSK_RES_OK )
			{
				fprintf(stderr, "Error %d at setting mosek double parameter: %s at the value: %f!\n", r, (*itDbl)->name, (*itDbl)->value );
				//getchar();
			}
		}
		
		for( itStr = params->strParams.begin(); itStr != params->strParams.end(); itStr++ )
		{
			r = MSK_putnastrparam(task, (*itStr)->name,  (*itStr)->value);
			if ( r != MSK_RES_OK )
			{
				fprintf(stderr, "Error %d at setting mosek string parameter: %s at the value: %s!\n", r, (*itStr)->name, (*itStr)->value );
				//getchar();
			}
		}
	} */
	
	
	void IQD_setMosekParams(MSKtask_t &task, IQD_GeneralSolverParams *params)
	{
		MSKrescodee r;
		
		for(auto &pair : params->intParams)
		{
			r = MSK_putnaintparam(task, pair.first.c_str(),  pair.second);
			if ( r != MSK_RES_OK )
			{
				fprintf(stderr, "Error %d at setting mosek integer parameter: %s at the value: %ld!\n", r, pair.first.c_str(), pair.second );
				//getchar();
			}
		}

		for(auto &pair : params->dblParams)
		{
			r = MSK_putnadouparam(task, pair.first.c_str(),  pair.second);
			if ( r != MSK_RES_OK )
			{
				fprintf(stderr, "Error %d at setting mosek double parameter: %s at the value: %f!\n", r, pair.first.c_str(), pair.second );
				//getchar();
			}
		}

		for(auto &pair : params->strParams)
		{
			r = MSK_putnastrparam(task, pair.first.c_str(),  pair.second.c_str());

			if ( r != MSK_RES_OK )
			{
				fprintf(stderr, "Error %d at setting mosek string parameter: %s at the value: %s!\n", r, pair.first.c_str(), pair.second.c_str() );
				//getchar();
			}
		}
	}
	
	
#endif








int IQD_SDPLambda::calculateLambda(const int sdpSolver, IQD_SparseMatrix& Q, const double Qfactor, const int nv, double** v, double* lambda, const double* objCoeflambda, IQD_GeneralSolverParams* sdpSolverParams, const double zero_tol)//, const double eps_diag)
{
	//const int n= Q.nrows;
	bool objLambdaAlloc = false;
	int i, code;
	double *objLambda = NULL;
	void *p = &objCoeflambda;
	
	
	
	if( objCoeflambda == NULL)
	{
		objLambda = (double *) malloc( nv * sizeof(double) );
		if(!objLambda)
		{
			#if IQD_DEBUG_MODE
				IQD_PRINTMEMERROR;
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
		}
		
		objLambdaAlloc = true;
		
		for(i = 0; i < nv; i++)
			objLambda[i] = 1.0;
	}
	else
	{
		objLambda =  *((double **) p) ;
	}
	
	
	
	switch(sdpSolver)
	{
		case IQD_SS_MOSEK:
			code = calculateLambdaMosek(Q, Qfactor, nv, v, lambda, objLambda, sdpSolverParams);
			break;
		
		case IQD_SS_CSDP:
			code = calculateLambdaCSDP(Q, Qfactor, nv, v, lambda, objLambda);
			break;
		
		default:
		    std::cerr << IQD_PREPRINT "Solver " << sdpSolver << " is not supported at IQD_SDPLAMBDA::calculateLambda. Please, send that message to " IQD_EMAIL "\n";
		    code = IQD_LIBRARY_NOT_AVAILABLE;
	}
	
	if( code != IQD_OPTIMAL_SOLUTION )
	{
		std::cout << IQD_PREPRINT  "sdp lambda code: " << code << "\n";
		IQD_getchar();
	}
	
	//cout << "zero_tol: " << zero_tol << endl;
	//getchar();
	
	for(i = 0; i < nv; i++)
	{
		if( lambda[i] >= -zero_tol )
			lambda[i] = 0.0;
		//else
			//lambda[i] -= eps_diag;
	}
	
	/*printf("lambda:\n");
	for(i = 0; i < nv; i++)
		printf("%f ", lambda[i]);
	getchar(); */
	
	
	
desallocate_memory:
	
	if(objLambdaAlloc)	free(objLambda);
	
	
	return code;
}




/*
 * Here, we solve the problem:
 * 
 * Max  -tr(Q)
 * s.t.:
 * 		tr(V_i X)		<= objCoefLambda
 * 
 * where V_i = v_i *v_i', i = 1, .., nv
 */



int IQD_SDPLambda::calculateLambdaMosek(IQD_SparseMatrix &Q, const double Qfactor, const int nv, double **v, double *lambda, const double *objCoefLambda, IQD_GeneralSolverParams *sdpSolverParams)
#if IQD_HAVE_MOSEK
{
	const int n = Q.getNumberOfRows(); // nrows;
	const int nsym = (n*(n+1))/2;
	const int nconstr = nv;
	
	int i, j, k, aux, code;
	
	
	MSKrescodee trmcode;
    MSKrescodee  r;
    MSKsolstae solsta;
	
	MSKint32t    DIMBARVAR[] = {n};
	
	MSKint32t *rows = NULL, *cols = NULL;
    double *vals = NULL;
    
    MSKint64t    idx;//, idx2;
    double       falpha = 1.0;
	
	
	MSKenv_t     env = NULL;
    MSKtask_t    task = NULL;
	
	
	/*Q.printSparseMatrix();
	
	for(int w = 0; w < nv; w++)
	{
		for(int x = 0; x < n; x++)
			printf("%0.1f  ", v[w][x]);
		printf("\n");
	} 
	
	printf("objCoefLambda: ");
	for(int w = 0; w < nv; w++)
		printf( "%f  ", objCoefLambda[w] );
	printf("\n"); */
	
	//for(int i = 0; i < nv; i++)
		//printf("objLambda: %f\n", objCoefLambda[i]);
	
	
	rows = (MSKint32t *) malloc( nsym * sizeof(MSKint32t) );
	cols = (MSKint32t *) malloc( nsym * sizeof(MSKint32t) );
	vals = (double *)  malloc( nsym * sizeof(double) );
	
	if( !rows || !cols || !vals )
    {
		#if IQD_DEBUG_MODE
			IQD_PRINTMEMERROR;
		#endif
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
	
	
	/* Create the mosek environment. */
    r = MSK_makeenv(&env,NULL);
    if(r != MSK_RES_OK)
    {
		#if IQD_DEBUG_MODE
			IQD_PRINTERRORNUMBER(r);
		#endif
		code = IQD_SDP_SOLVING_ERROR;
		goto desallocate_memory;
    }
    
    //printf("SDP vector Mosek 2\n");
    
    /* Create the optimization task. */
    r = MSK_maketask(env, nconstr, 0, &task);
    if(r != MSK_RES_OK)
    {
		#if IQD_DEBUG_MODE
			IQD_PRINTERRORNUMBER(r);
		#endif
		code = IQD_SDP_SOLVING_ERROR;
		goto desallocate_memory;
    }
	
	
	IQD_setMosekParams(task, sdpSolverParams);
	
	 /* Append 'nconstr' empty constraints.
       The constraints will initially have no bounds. */
    
    r = MSK_appendcons(task, nconstr);
    if(r != MSK_RES_OK)
    {
		#if IQD_DEBUG_MODE
			IQD_PRINTERRORNUMBER(r);
		#endif
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
	
	//printf("SDP vector Mosek 3\n");
	
	//we have only 1 sdp variable
	r = MSK_appendbarvars(task, 1, DIMBARVAR);
    if(r != MSK_RES_OK)
    {
		#if IQD_DEBUG_MODE
			IQD_PRINTERRORNUMBER(r);
		#endif
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
	
	{
		
		const int *pcols = Q.getRowColsPointer(0);
		const double *pvals;
		//seting objective function. here, we work on a minimization problem, so -tr(QZ) turns tr(QZ)
		
		aux = Q.getStructure(rows, NULL);
		
		
		if( Qfactor != 1.0 )
		{
			Q.getValues(vals);
			
			//for( int w = 0; w < aux; w++ )
				//vals[w] *= Qfactor;
			
			IQD_multiplyAllArray(aux, vals, Qfactor);
			pvals = vals;
		}
		else
		{
			pvals = Q.getRowValuesPointer(0);
		}
		
		
		r = MSK_appendsparsesymmat(task, DIMBARVAR[0], aux, rows, pcols, pvals, &idx);
		if(r != MSK_RES_OK)
		{
			#if IQD_DEBUG_MODE
				IQD_PRINTERRORNUMBER(r);
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
		}
		
		//falpha = -1.0;
		//we take advantage of minimization objective function and do no multiply by -1
		r = MSK_putbarcj(task, 0, 1, &idx, &falpha);
		if(r != MSK_RES_OK)
		{
			#if IQD_DEBUG_MODE
				IQD_PRINTERRORNUMBER(r);
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
		}
		//falpha = 1.0;
	}
	
    
    //setting contraints
    
    /*for(i = aux = 0; i < n; i++)
	{
		for(j = 0; j <= i; j++)
		{
			rows[aux] = i;
			cols[aux] = j;
			aux++;
		}
	} */
	
	//printf("SDP vector Mosek 5\n");
    
    //setting matrices
    for(k = 0; k < nconstr; k++)
	{
		
		/*for(i = aux = 0; i < n; i++)
		{
			for(j = 0; j <= i; j++)
			{
				vals[aux] = v[k][i]*v[k][j] ;
				
				if(vals[aux] != 0.0)
				{
					rows[aux] = i;
					cols[aux] = j;
					aux++;
				}
			}
		}*/
		
		const double *vk = v[k];
		
		for(i = aux = 0; i < n; i++)
		{
			if( vk[i] != 0.0 )
			{
				for(j = 0; j <= i; j++)
				{
					if( vk[j] != 0.0 )
					{
						vals[aux] = vk[i]*vk[j] ;
						rows[aux] = i;
						cols[aux] = j;
						aux++;
					}
				}
			}
		}
		
		//printf("k: %d aux: %d\n", k, aux);
		
		r = MSK_appendsparsesymmat(task, DIMBARVAR[0], aux, rows, cols, vals, &idx);
		if(r != MSK_RES_OK)
	    {
			#if IQD_DEBUG_MODE
				IQD_PRINTERRORNUMBER(r);
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
	    }
		
		
		r = MSK_putbaraij(task, k, 0, 1, &idx, &falpha);
		if(r != MSK_RES_OK)
	    {
			#if IQD_DEBUG_MODE
				IQD_PRINTERRORNUMBER(r);
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
	    }
	    
	    //setting rhs
	    //r = MSK_putconbound(task, k, MSK_BK_UP, 1.0, 1.0);
		r = MSK_putconbound(task, k, MSK_BK_UP, 0.0, objCoefLambda[k]);
		
		#if IQD_DEBUG_MODE
			if(r != MSK_RES_OK)
			{
				#if IQD_DEBUG_MODE
					IQD_PRINTERRORNUMBER(r);
				#endif
				code = IQD_MEMORY_ERROR;
				goto desallocate_memory;
			}
		#endif
	}
    
    //printf("SDP vector Mosek 6\n");
    
    r = MSK_optimizetrm(task,&trmcode);
	
	if( r != MSK_RES_OK )
	{
		#if OPT_DEBUG_MODE
			IQD_PRINTERRORNUMBER(r);
		#endif
			
		code = IQD_SOLVER_ERROR;
		goto desallocate_memory;
	}
	
	
	
	
	MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
	
	//printf("SDP Mosek solution status: %d\n", solsta);
	
	
	
	switch(solsta)
    {
		case MSK_SOL_STA_OPTIMAL:
		case MSK_SOL_STA_NEAR_OPTIMAL:
		    
		case MSK_SOL_STA_DUAL_FEAS:
		case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
		    
		case MSK_SOL_STA_NEAR_DUAL_FEAS:
			
			
			MSK_getsolution(task, MSK_SOL_ITR, NULL, NULL, NULL, NULL, NULL, NULL, NULL, lambda, NULL, NULL, NULL, NULL, NULL );
			
			MSK_getprimalobj(task, MSK_SOL_ITR, &falpha);
		    
		    //printf("Mosek Objective function: %f\n", falpha);
			
			code = IQD_OPTIMAL_SOLUTION;
			
			break;
			
		default:
			std::cerr << IQD_PREPRINT "SDP Mosek solution status: " << solsta <<  " termination code: " << trmcode << "\n";
		    code = IQD_GENERIC_SOLVER_ERROR;
	}
	
	
	
desallocate_memory:
	
	
	if(rows)	free(rows);
	if(cols)	free(cols);
    if(vals)	free(vals);
    
    if(task)	MSK_deletetask(&task);
    if(env)     MSK_deleteenv(&env);
	
	
	
	
	return code;
}
#else
{
    return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif






int IQD_SDPLambda::calculateLambdaCSDP( IQD_SparseMatrix &Q, const double Qfactor, const int nv, double **v, double *lambda, const double *objCoefLambda)
#if IQD_HAVE_CSDP
{
	int i, j, k, aux, code;
	const int n = Q.nrows;
	const int msdp = nv;
	const int nsym = (n*(n+1))/2;
	double objF, dobj; // auxD;
	
	int *cols = NULL;
	double *vals = NULL;
	
	//variables to CSDP problem
    //Auxiliar variable to set block in constraints...
    struct sparseblock *blockptr;
    
    //The problem and solution data.
    struct blockmatrix C;
    double *b;
    struct constraintmatrix *constraints;
    
    // Storage for the initial and final solutions.
    struct blockmatrix X,Z;
    double *lagran;
	
	
	
	cols = (int *) malloc( n * sizeof(int) );
	vals = (double *) malloc( n* sizeof(double));
	b = (double *)malloc( (msdp + 1) *sizeof(double) );
	
	if( !cols || !vals || !b )
	{
		#if IQD_MEMORY_ERROR
			IQD_PRINTMEMERROR;
		#endif
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
	}
	
	
	//first, allocate storage for the C matrix. We have two blocks, but need to allocate for 3;
	
	//The first block is our SDP variable. The second block are slack variable since CSDP only accepts equality constraints
	
	C.nblocks = 2;
    C.blocks = (struct blockrec *) malloc( 3 * sizeof(struct blockrec) );
    if( !C.blocks )
    {
		#if IQD_MEMORY_ERROR
			IQD_PRINTMEMERROR;
		#endif
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
	
	
	C.blocks[1].blockcategory = MATRIX;
    C.blocks[1].blocksize = n;
    C.blocks[2].blockcategory = DIAG;
    C.blocks[2].blocksize = msdp;
	
	
	C.blocks[1].data.mat = (double *) calloc( n*n, sizeof(double) );
	C.blocks[2].data.vec = (double *) calloc( msdp+1, sizeof(double) );
	
	
	if(!C.blocks[2].data.vec || !C.blocks[1].data.mat)
    {
		#if IQD_MEMORY_ERROR
			IQD_PRINTMEMERROR;
		#endif
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
	
	for(i = 1; i <= n; i++)
	{
		const int nzrow = Q.getNumberOfElementsAtRow(i-1);
		
		const int* pcols = Q.getRowColsPointer(i-1);
		const double* pvals;
		
		//Q.getRowStructureAndValues(i-1, aux, cols, vals);
		
		if( Qfactor != 1.0 )
		{
			Q.getRowValues(i-1, vals, NULL);
			for(int w = 0; w < nzrow; w++)
				vals[w] *= Qfactor;
			
			pvals = vals;
		}
		else
		{
			pvals = Q.getRowValuesPointer(i-1);
		}
		
		
		for(int j = 0; j < nzrow; j++)
		{
			C.blocks[1].data.mat[ijtok(i, pcols[j] +1 ,n)] = -pvals[j];
			
			C.blocks[1].data.mat[ijtok(pcols[j] +1, i, n)] = -pvals[j];
		}
	}
	
	
    //setting the constraint matrices
    constraints = (struct constraintmatrix *) malloc( (msdp + 1) * sizeof(struct constraintmatrix) ) ;
    if( !constraints )
    {
		#if IQD_MEMORY_ERROR
			IQD_PRINTMEMERROR;
		#endif
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
	
	
	for(int k = 1; k <= msdp; k++)
		b[k] = objCoefLambda[k-1];//1.0;
	
	
	for(k = 0; k < msdp; k++)
	{
		blockptr = (struct sparseblock *)malloc( sizeof(struct sparseblock) );
	    if(!blockptr)
	    {
			#if IQD_MEMORY_ERROR
				IQD_PRINTMEMERROR;
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
	    }
	    
		//insert the block in the constraint
	    constraints[k+1].blocks = blockptr;
		
		blockptr->blocknum = 1;
		blockptr->blocksize = n;
		blockptr->constraintnum = k + 1;
		blockptr->nextbyblock = NULL;
		
		
		aux = nsym;
		
		blockptr->entries = (double *) malloc( (aux+1) * sizeof(double) );
		blockptr->iindices = (int *) malloc( (aux+1) * sizeof(int) );
		blockptr->jindices = (int *) malloc( (aux+1) * sizeof(int) );
		
		if( !(blockptr->entries) || !(blockptr->iindices) || !(blockptr->jindices) )
	    {
			#if IQD_MEMORY_ERROR
				IQD_PRINTMEMERROR;
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
	    }
		
		
		for(i = aux = 0; i < n; i++)
		{
			const double *vk = v[k];
			if(vk[i] != 0.0) 
			{
				for(j = i; j < n; j++) //upper triangle
				{
					
					if( vk[j] != 0.0 )
					{
						//auxD = v[k][i]*v[k][j];
						
						//if( auxD != 0.0 )
						{
							aux++;
							blockptr->iindices[aux] = i+1;
							blockptr->jindices[aux] = j+1;
							blockptr->entries[aux] = vk[i]*vk[j]; //auxD;
						}
					}
				}
			}
		}
		
		blockptr->numentries = aux; //only the upper triangle
		
		
		blockptr->next = (struct sparseblock *)malloc( sizeof(struct sparseblock) );
	    if(!blockptr)
	    {
			#if IQD_MEMORY_ERROR
				IQD_PRINTMEMERROR;
			#endif
			code = IQD_MEMORY_ERROR;
			goto desallocate_memory;
	    }
		
		blockptr = blockptr->next;
		
		
		// Now, the last blockcategory, i.e., the slack variables
		blockptr->blocknum = 2;
		blockptr->blocksize = msdp;
		blockptr->constraintnum = k+1;
		blockptr->next = NULL;
		blockptr->nextbyblock = NULL;
		blockptr->numentries = 1;
		
		blockptr->entries = (double *) malloc( (1+1)*sizeof(double) );
		blockptr->iindices = (int *) malloc( (1+1)*sizeof(int) );
		blockptr->jindices = (int *) malloc( (1+1)*sizeof(int) );
		
		
		if( !blockptr->entries || !blockptr->iindices || !blockptr->jindices )
		{
			#if IQD_MEMORY_ERROR
				IQD_PRINTMEMERROR;
			#endif
		    code = IQD_MEMORY_ERROR;
		    goto desallocate_memory;
		}
		
		blockptr->iindices[1] = blockptr->jindices[1] = k + 1;
		
		blockptr->entries[1] = 1.0; //first block
	}
	
	
	//writting the prob in the SDPA format
    //write_prob("prob.dat-s", n+msdp, msdp, C, b, constraints);
    
       
    /*
    * Create an initial solution.  This allocates space for X, y, and Z,
    * and sets initial values.
    */
	initsoln(n+msdp,msdp,C,b,constraints,&X,&lagran,&Z);
	
	i = easy_sdp(n+msdp, msdp, C, b, constraints, 0.0, &X, &lagran, &Z, &objF, &dobj);
	
	printf("CSDP return code: %d\n", i);
	
	switch(i)
    {
		case 0:
		case 3:
		    //getting the solution
		    
		    for(i = 0; i < nv; i++)
		    {
				//we must multiply by -1 because p(x) = -r_i x_i ^2
				lambda[i] = -lagran[i+1];
		    }
		    
		    
		    code = IQD_OPTIMAL_SOLUTION;
		    break;
		    
		case 1:
		case 2:
		    code = IQD_INFEASIBLE_PROBLEM;
		    break;
		    
		case 4:
		    code = IQD_MAX_ITERATIONS_STOP;
		    break;
		    
		default:
		    code = IQD_GENERIC_SOLVER_ERROR;
		    break;
    }

    free_prob(n+msdp, msdp, C, b, constraints, X, lagran, Z);
	
	
	
desallocate_memory:
	
	
	if(cols)	free(cols);
	if(vals)	free(vals);
	

	return code;
	
}
#else
{
    return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif






















