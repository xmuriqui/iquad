/*
 * functions to read quadratic programs from:
 * 
 * 
 * 
 * 
 */

#include <cstdio>

#include <new>

#include "iquad.hpp"
#include "IQD_tools.hpp"






using namespace iquad;
using namespace std;





int iquad::IQD_readQuadProg(const char *fileName, const char *outDirectory)
{
	int returnValue, i, j, m, n;
	const int lineSize = 100;
	char *line = NULL;
	double *M = NULL;
	const double lbx = -1000, ubx = 1000;
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
    
    
    line = (char *) malloc(lineSize * sizeof(char));
	if(!lineSize)
	{
		returnValue = IQD_MEMORY_ERROR;
		goto desallocate_memory;
	}
    
    do
	{
	    fgets(line, lineSize, file);
		
		if(feof(file))
		{
			returnValue = IQD_BAD_DEFINITIONS;
			goto desallocate_memory;
		}
		
	}while( fscanf(file, "%d%d", &n, &m) == 0);
	
    printf("n: %d m: %d\n", n,m);
	
    j = IQD_max(n, m);
	M = (double *) malloc( (j*n + 1) * sizeof(double));
	
	if( !M )
	{
		returnValue = IQD_MEMORY_ERROR;
		goto desallocate_memory;
	}
	
	i = prob.setParametersAndAllocate(n, m);
	if(i != 0)
	{
		returnValue = IQD_MEMORY_ERROR;
		goto desallocate_memory;
	}
	
    
    //reading c
    do
	{
	    fgets(line, lineSize, file);
		
		if(feof(file))
		{
			returnValue = IQD_BAD_DEFINITIONS;
			goto desallocate_memory;
		}
		
	}while(fscanf(file, "%lf", &M[0]) == 0); 
	
	
	
	for(i = 1; i < n; i++)
		fscanf(file, "%lf", &M[i]);
	
	prob.setObjLinearCoefficients(M);
	
	
	//reading Q
	do
	{
		fgets(line, lineSize, file);
		
		if(feof(file))
		{
			returnValue = IQD_BAD_DEFINITIONS;
			goto desallocate_memory;
		}
	}while(fscanf(file, "%lf", &M[0]) == 0);
	
    
	for(j = 1; j < n; j++)
		fscanf(file, "%lf", &M[j]);
	
	for(i = 1; i < n; i++)
	{
		for(j = 0; j < n; j++)
			fscanf(file, "%lf", &M[i*n + j]);
	}
	
	i = prob.setObjQuadCoefsMatrix(M);
	if(i)
	{
		returnValue = IQD_MEMORY_ERROR;
		goto desallocate_memory;
	}
    
    if(m > 0)
	{
		//reading A
	    do
		{
			fgets(line, lineSize, file);
			
			if(feof(file))
			{
				returnValue = IQD_BAD_DEFINITIONS;
				goto desallocate_memory;
			}
		}while(fscanf(file, "%lf", &M[0]) == 0);
	    
	    
	    for(j = 1; j < n; j++)
			fscanf(file, "%lf", &M[j]);
		
		i = prob.setConstraintLinearPart(0, M);
		if(i)
		{
			returnValue = IQD_MEMORY_ERROR;
			goto desallocate_memory;
		}
		
		for(i = 1; i < m; i++)
		{
			for(j = 0; j < n; j++)
				fscanf(file, "%lf", &M[j]);
			
			j = prob.setConstraintLinearPart(i,M);
			if(j)
			{
				returnValue = IQD_MEMORY_ERROR;
				goto desallocate_memory;
			}
		}
		
		
		//reading b
		do
		{
			fgets(line, lineSize, file);
			
			if(feof(file))
			{
				returnValue = IQD_BAD_DEFINITIONS;
				goto desallocate_memory;
			}
		}while(fscanf(file, "%lf", &M[0]) == 0);
		
		for(i = 1; i < m; i++)
			fscanf(file, "%lf", &M[i]);
		
		i = prob.setConstraintUpperBounds(M);
		if(i)
		{
			returnValue = IQD_MEMORY_ERROR;
			goto desallocate_memory;
		}
		
		
		//for(i = 0; i < m; i++)
			//prob.setQuadConstraintType(i, IQD_CT_LEEQ);
		
	}
    
    
    for(i = 0; i < n; i++)
	{
		prob.setVariableLowerBound(i, 0.0); //VarBounds(i, 0, ubx);
		prob.setVariableUpperBound(i, ubx);
	}
    
    
    free(M); M = NULL;
	free(line); line = NULL;
	
	
	//prob.print();
	
	
	IQD_runProblem( fileName, outDirectory, prob, bb, &sParams, &mosekSDPParams );
	
	returnValue = bb.out_return_code;
	
	
    
desallocate_memory:
	
	
	
	if(file) fclose(file);
	
	if(M) free(M);
	if(line) free(line);

	return returnValue;
}


