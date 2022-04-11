#include <cstdio>
#include <new>

#include "iquad.hpp"
#include "IQD_tools.hpp"




using namespace minlpproblem;
using namespace iquad;
using namespace std;


int iquad::IQD_biqmac(const char* fileName, const char *outDirectory, const bool considerIntegrality, const bool maxcut)
{
    int returnValue, i, nz, n;
    int *rows = NULL, *cols = NULL;
    double *vals = NULL;
    FILE *file = NULL;
    IQD_IQuadProb prob;
    IQD_GeneralSolverParams sParams, mosekSDPParams;
    IQD_BranchAndBound bb;
	
	
	
	
	file = fopen(fileName, "r");
	if(!file)
	{
		printf("Cannot open %s\n", fileName);
		returnValue = IQD_UNDEFINED_ERROR;
		goto desallocate_memory;
	}
    
    fscanf(file, "%d%d", &n, &nz);
    
    i = IQD_max(n, nz);
    rows = (int *) malloc( i * sizeof(int) );
    cols = (int *) malloc( i * sizeof(int) );
    vals = (double *) malloc( i * sizeof(double) );
    
    i = prob.setParametersAndAllocate(n, 0);
    if(!rows || !cols || !vals || i != 0)
    {
		returnValue = i;
		goto desallocate_memory;
    }
    
    if(considerIntegrality)
	{
	    for(i = 0; i < n; i++)
		{
			prob.setVariableType(i, MIP_VT_INTEGER);
			prob.setVariableLowerBound(i, 0.0);
			prob.setVariableUpperBound(i, 1.0);
		}
    }
    else
	{
		if( maxcut )
		{
			for(i = 0; i < n; i++)
			{
				prob.setVariableLowerBound(i, -1.0);
				prob.setVariableUpperBound(i, 1.0);
			}
		}
		else
		{
			for(i = 0; i < n; i++)
			{
				prob.setVariableLowerBound(i, 0.0);
				prob.setVariableUpperBound(i, 1.0);
			}
		}
	}
    
    
    for(i = 0; i < nz; i++)
    {
		fscanf(file, "%d%d%lf\n", &rows[i], &cols[i], &vals[i]);
		//printf("vals[%d]: %f ", i, vals[i]);
		rows[i]--;
		cols[i]--;
		vals[i] = 2.0*vals[i];
		//printf("rows[%d]: %d cols[%d]: %d vals[%d]: %f\n", i, rows[i], i, cols[i], i, vals[i]);
    }
    
    fclose(file);
    file = NULL;
    
    prob.setObjQuadCoefsMatrix(nz, rows, cols, vals);
    
    //prob.print();
    
    //prob.Q.printAllMatrix();
    
    /*{
		char outAmpl[200];
		
		sprintf(outAmpl, "%s.mod", fileName );
		prob.writeMIQCPModelInAMPLFile(outAmpl, "couenne");
		goto desallocate_memory;
	} */
    
    
    /*
    {
		IQD_OutFile outFile2;
		mosekParams.addIntegerParameter("MSK_IPAR_NUM_THREADS", 1 );
		
		outFile2.open("output_eigen.txt");
		outFile2.writeBreakLine();
		outFile2.writeProblemInfo(fileName, prob);
		
		IQD_generateRstatistic(fileName, prob, IQD_SQ_SCHUR, &mosekParams,outFile2, 1.0e-5);
		IQD_generateRstatistic(fileName, prob, IQD_SQ_MIN_EIG_SDP, &mosekParams,outFile2, 1.0e-5);
		IQD_generateRstatistic(fileName, prob, IQD_SQ_2_BLOCK_SDP, &mosekParams,outFile2, 1.0e-5);
		IQD_generateRstatistic(fileName, prob, IQD_SQ_DIAG_SDP, &mosekParams,outFile2, 1.0e-5);
		IQD_generateRstatistic(fileName, prob, IQD_SQ_DIAG_DOM, &mosekParams,outFile2, 1.0e-5);
		IQD_generateRstatistic(fileName, prob, IQD_SQ_IDENTITY, &mosekParams,outFile2, 1.0e-5);
		
		//goto desallocate_memory;
    }*/
    
    /*
    IQD_writeEigInfo(fileName, prob.Q, "output_eigen2.txt", 1.0e-5);
    
    goto desallocate_memory;
	*/
	
	
	IQD_runProblem( fileName, outDirectory, prob, bb, &sParams, &mosekSDPParams );
	
	
	returnValue = bb.out_return_code;
	
	
	
desallocate_memory:
	
	
	if(file)	fclose(file);
	if(rows)	free(rows);
	if(cols)	free(cols);
	if(vals)	free(vals);
	
	return returnValue;
}
