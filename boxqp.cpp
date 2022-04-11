#include <cstdio>
#include <new>

#include "iquad.hpp"
#include "IQD_tools.hpp"

using namespace iquad;
using namespace std;



class IQD_BoxqpNlEval:public IQD_NonLinearEval
{
	virtual int eval_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double &value)
	{
		value = 0.0;
		return 0;
	}
	
	virtual int eval_nl_constrs_part(const int threadnumber, const int n, const int mnl, const bool newx, const bool *constrEval, const double *x, double *values)
	{
		return 0;
	}
	
	virtual int eval_grad_nl_obj_part(const int threadnumber, const int n, const bool newx, const double *x, double *values)
	{
		for(int i = 0; i < n; i++)
			values[i] = 0.0;
		
		return 0;
	}
	
	virtual int eval_grad_nl_constrs_part(const int threadnumber, const int n, const int mnl, const int nz, const bool newx, const bool *constrEval, const double *x, IQD_SparseMatrix& jacobian)
	{
		return 0;
	}
	
	virtual int eval_hessian_nl_lagran_part(const int threadnumber, const int n, const int mnl, const int nz, const bool newx, const double *x, const double objFactor, const double *lambda, IQD_SparseMatrix& hessian)
	{
		return 0;
	}
	
};


int iquad::IQD_boxqp(const char* fileName, const char *outDirectory)
{
	//char *gamsModelName = NULL;
	int returnValue, i, j, k, nz, n;
	int *rows = NULL, *cols = NULL;
	double aux, *vals = NULL;
	FILE *file = NULL;
	
	IQD_IQuadProb prob;
	IQD_GeneralSolverParams sParams, mosekSDPParams;
	IQD_BranchAndBound bb;
	IQD_BoxqpNlEval eval;
	
	
	
	file = fopen(fileName, "r");
	if(!file)
	{
		printf("Cannot open %s\n", fileName);
		returnValue = IQD_UNDEFINED_ERROR;
		goto desallocate_memory;
	}

	fscanf(file, "%d", &n);

	nz = ((1 + n)*n)/2;

	rows = (int *) malloc( nz * sizeof(int) );
	cols = (int *) malloc( nz * sizeof(int) );
	vals = (double *) malloc( nz * sizeof(double) );
	
	
	i = prob.setParametersAndAllocate(n, 0);
	if(!rows || !cols || !vals || i != 0)
	{
		returnValue = i;
		goto desallocate_memory;
	}

	for(i = 0; i < n; i++)
	{
		fscanf(file, "%lf", &vals[i]);
		vals[i] = -vals[i];
	}

	prob.setObjLinearCoefficients(vals);

	for(k = 0, i = 0; i < n; i++)
	{
		for(j = 0; j <= i; j++, k++)
		{
			fscanf(file, "%lf", &vals[k] );
			vals[k] = -vals[k]; //maximization problems...
		}

		//fscanf(file, "%lf", &aux); //diagonal position
		//vals[k] = aux;
		//k++;

		for( ;j < n; j++)
			fscanf(file, "%lf", &aux);
	}
	
	fclose(file);
	file = NULL;
	
	
	for(i = k =  0; i < n; i++)
	{
		for(j = 0; j <= i; j++, k++)
		{
			rows[k] = i;
			cols[k] = j;
		}
	}
	
	prob.setObjQuadCoefsMatrix(nz, rows, cols, vals);
	
	for(i = 0; i < n; i++)
	{
		prob.setVariableLowerBound(i, 0.0);
		prob.setVariableUpperBound(i, 1.0);
	}
	
	//prob.printProblem();
	
	
	free(rows);		rows = NULL;
	free(cols);		cols = NULL;
	free(vals);		vals = NULL;
	
	
	
	
	
	
	
	/*
	{
	double dif;
	IQD_RefProb rp1, rp2;

	rp1.prob = &prob;
	rp2.prob = &prob;

	rp1.setWayOfSplitQ(IQD_SQ_SDP);
	rp2.setWayOfSplitQ(IQD_SQ_SCHUR);

	rp1.splitQ();
	rp2.splitQ();

	dif = IQD_maxDifferenceBetweenMatrices(prob.n, prob.n, rp1.R, rp2.R);

	outFile->writeDouble(dif);

	goto desallocate_memory;
	}
	*/
	
	/*
	{
		IQD_OutFile outFile2;
		
		mosekParams.addIntegerParameter("MSK_IPAR_NUM_THREADS", 1 );
		
		outFile2.open("output_eigen.txt");
		outFile2.writeBreakLine();
		outFile2.writeProblemInfo(fileName, prob);
		
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_SCHUR, &mosekParams, outFile2, 1.0e-5);
		IQD_generateRstatistic(fileName, prob, IQD_SQ_MIN_EIG_SDP, &mosekParams,outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_2_BLOCK_SDP, &mosekParams,outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_DIAG_SDP, &mosekParams,outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_DIAG_DOM, &mosekParams,outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_IDENTITY, &mosekParams,outFile2, 1.0e-5);
		
		//goto desallocate_memory;
	} */
	
	
	//IQD_writeEigInfo(fileName, prob.Q, "output_eigen2.txt", 1.0e-5);
	
	//goto desallocate_memory;
	
	
	
    /*
    gamsModelName = (char *) malloc( (strlen(fileName) + 5) * sizeof(char) );
    if(!gamsModelName)
    {
		returnValue = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
    
    strcpy( gamsModelName, fileName);
    strcat( gamsModelName, ".gms" );
    
    
    prob.writeMIQCPModelInGAMSFile(gamsModelName, fileName , "baron", true);
    
    goto desallocate_memory;
    */
    
	
	
	/*{
		int nzQ = prob.Q.nElements;
		int rows[nzQ], cols[nzQ];
		double vals[nzQ];
		double eig[prob.n];
		char fileName[80];
		
		
		IQD_replaceZeroEigValuesOfQByRandom(prob, 1.0e-4);
		
		nzQ = prob.Q.getStructureAndValues(rows, cols, vals);
		
		for(int w = 1; w <= prob.n; w++ )
		{
			IQD_changeQeigenvalues(prob, w, 1.0e-4);
			
			sprintf( fileName, "box3_myqp_%03d_eig%02d.mod", prob.n, w );
			prob.writeMIQCPModelInAMPLFile(fileName);
			
			sprintf( fileName, "box3_myqp_%03d_eig%02d.in", prob.n, w );
			prob.writeMIQCPModelInFile(fileName);
			
			
			//IQD_getEigenValues(prob.Q, eig);
			
			//for(int j = 0; j < prob.n; j++)
				//printf("eig[%d]: %f ", j, eig[j]);
			//printf("\n\n");
			//getchar();
			
			prob.Q.setStructureAndValues(nzQ, rows, cols, vals);
		}
	} */
	
	
	
	
	/*{
		char outAmpl[200];
		
		sprintf(outAmpl, "%s.mod", fileName );
		prob.writeMIQCPModelInAMPLFile(outAmpl, "couenne");
		goto desallocate_memory;
	} */
	
	
	
	
	IQD_runProblem( fileName, outDirectory, prob, bb, &sParams, &mosekSDPParams );
	
	
	returnValue = bb.out_return_code;
	
	
desallocate_memory:
	
	
	if(rows)	free(rows);
	if(cols)	free(cols);
	if(vals)	free(vals);
	if(file)	fclose(file);

	return returnValue;
}


