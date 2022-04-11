


#include <iostream>
#include <new>

#include <cstdio>
#include <cstdlib>

#include "iquad.hpp"
#include "IQD_qpGen.hpp"
#include "IQD_tools.hpp"


using namespace iquad;
using namespace std;



class IQD_MyqpNlEval:public IQD_NonLinearEval
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





int iquad::IQD_myqp(const char* fileName, const char *outDirectory)
{
	int i, returnValue;
	IQD_IQuadProb prob;
	IQD_BranchAndBound bb;
	IQD_GeneralSolverParams sParams, mosekSDPParams;
	IQD_MyqpNlEval eval;
	
	
	
	
	i = prob.readMIQCPModelInFile(fileName);
	
	if(i)
	{
		fprintf(stderr, "Problem %d to read problem in file %s.", i, fileName );
		returnValue = i;
		goto desallocate_memory;
	}
	
	/*{
		char fileOut[100];
		int e = strlen(fileName);
		
		strcpy(fileOut, fileName);
		
		//replacing extension .in by .gms
		
		fileOut[e+1] = '\0';
		fileOut[e] = 's';
		fileOut[e-1] = 'm';
		fileOut[e-2] = 'g';
		
		
		prob.writeMIQCPModelInGAMSFile(fileOut, fileName, "couenne", true);
		
		goto desallocate_memory;
	}*/
	
	
	
	/*{
		char outAmpl[200];
		
		sprintf(outAmpl, "%s.mod", fileName );
		prob.writeMIQCPModelInAMPLFile(outAmpl, "couenne");
		goto desallocate_memory;
	} */
	
	
	
	/*{
		IQD_OutFile outFile2;
		
		mosekParams.addIntegerParameter("MSK_IPAR_NUM_THREADS", 1 );
		
		outFile2.open("output_eigen.txt");
		outFile2.writeBreakLine();
		outFile2.writeProblemInfo(fileName, prob);
		
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_SCHUR, &mosekParams, outFile2, 1.0e-5);
		IQD_generateRstatistic(fileName, prob, IQD_SQ_MIN_EIG_SDP, &mosekParams, outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_2_BLOCK_SDP, &mosekParams,outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_DIAG_SDP, &mosekParams, outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_DIAG_DOM, &mosekParams, outFile2, 1.0e-5);
		//IQD_generateRstatistic(fileName, prob, IQD_SQ_IDENTITY, &mosekParams, outFile2, 1.0e-5);
		
		//goto desallocate_memory;
    }*/
	
	/*{
		int nzQ = prob.Q.nElements;
		int rows[prob.n * prob.n], cols[prob.n * prob.n];
		double vals[prob.n * prob.n];
		double eig[prob.n];
		char fileName[80];
		
		
		IQD_replaceZeroEigValuesOfQByRandom(prob, 1.0e-4);
		
		nzQ = prob.Q.getStructureAndValues(rows, cols, vals);
		
		for(int w = 1; w <= prob.n; w++ )
		{
			IQD_changeQeigenvalues(prob, w, 1.0e-4);
			
			sprintf( fileName, "ctr3_myqp_%03d_eig%02d.mod", prob.n, w );
			prob.writeMIQCPModelInAMPLFile(fileName);
			
			sprintf( fileName, "ctr3_myqp_%03d_eig%02d.in", prob.n, w );
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
		const int n = prob.n;
		int i;
		double eigs[n];
		double traceR, traceQ, sumNegEigQ;
		FILE *out = NULL;
		
		
		IQD_getEigenValues(prob.Q, eigs);
		
		
		traceR = -eigs[0]*n;
		traceQ =  eigs[0];
		
		for(i = 1; i < n; i++)
			traceQ += eigs[i];
		
		
		sumNegEigQ = 0.0;
		for(i = 0; i < n; i++)
		{
			if( eigs[i] < 0.0 )
				sumNegEigQ += eigs[i];
		}
		
		out = fopen("out_eig.txt", "a" );
		
		if(out)
		{
			fprintf(out, "\n%s%s%f%s%f%s%f%s%f%s", fileName, IQD_CHAR_SEP, traceQ, IQD_CHAR_SEP, traceR, IQD_CHAR_SEP, traceR - traceQ, IQD_CHAR_SEP, sumNegEigQ, IQD_CHAR_SEP);
			
			fclose(out);
		}
		
		
		//prob.Q.printAllMatrix();
		
		//getchar();
		
		goto desallocate_memory;
	} */
	
	
	//IQD_writeEigInfo(fileName, prob.Q, "output_eigen2.txt", 1.0e-5);
	//goto desallocate_memory;
	
	
	//bb.in_eps_added_to_diag_P = 1.0e-1;
	
	IQD_runProblem(fileName, outDirectory, prob, bb, &sParams, &mosekSDPParams);
	
	
	
	returnValue = bb.out_return_code;
	
	
desallocate_memory:
	
	return returnValue;
}



