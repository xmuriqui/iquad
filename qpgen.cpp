/*
 * Implementation of quadratic program generator from:
 * 
 * P. H. Calamai, L. N. Vicente and J. J. Judice, A new technique for generating quadratic programming test programming test problems
 * Mathematical Programming 61 (1993), 215-231
 * 
 * 
 * Author: Wendel Melo
 * 
 * Date: 21-May-2013
 * 
*/

# include <cstdio>
# include <cmath>
# include <ctime>
#if IQD_HAVE_RANDOMA
# include "sfmt.h"
#endif
# include "iquad.hpp"
# include "IQD_qpGen.hpp"
# include "IQD_tools.hpp"



using namespace iquad;
    

   
    
    

IQD_QPGenerator::IQD_QPGenerator()
{
	#if IQD_HAVE_RANDOMA
		int seed = 1986;
		initialize(seed);
	#endif
}


IQD_QPGenerator::IQD_QPGenerator(const int seed)
{
	#if IQD_HAVE_RANDOMA
		initialize(seed);
	#endif
}

int IQD_QPGenerator::initialize(const int seed)
#if IQD_HAVE_RANDOMA
{
	L = 4;
    condd = 1;
	seedToRandomNumbers = seed; 
    denv = 0.5;
	factor = 10.0;
	lowerBoundToGenerateD = 10;
	
	setSeed(&seedToRandomNumbers);
	
	return 0;
}
#else
{
	return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif





int IQD_QPGenerator::setSeed(const int *seed)
#if IQD_HAVE_RANDOMA
{
	return random.setSeed(seed);
}
#else
{
	return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif




	
int IQD_QPGenerator::generateProblem(const int nConcaveSubProbs, const int nBilinearSubProbs, const int nConvexSubProbs, IQD_IQuadProb &prob)
#if IQD_HAVE_RANDOMA
{
    const int m = nConcaveSubProbs + nBilinearSubProbs + nConvexSubProbs;
    const int n = 2*m;
    const int cons = 3*m;
    int l, i, j, theta, numNzv;
    int code;
    //int *type = NULL;
    double *Q = NULL, *s = NULL, *A = NULL, *c = NULL, *M = NULL, *v = NULL, *d = NULL;
    double *AM = NULL, *MQM = NULL, *Ms = NULL; 
    //IQD_Random random(seedToRandomNumbers);
    
    double sl, q;
    //double objfactor;
	double alpha, beta, xi = 0;
    double ro, omega;
    
    A = (double *) calloc( cons*n, sizeof(double) );
    c = (double *) malloc( cons * sizeof(double) );
    Q = (double *) calloc( n*n, sizeof(double) );
    s = (double *) calloc( n, sizeof(double) );
    
    if( !A || !c || !Q || !s )
    {
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
    
    
    
    for(l = 0; l < nConcaveSubProbs; l++)
    {
		theta = random.randInt(0, 1);
	
		if(theta == 0)
		{
		    q = -pow(4, l-L);
		    sl = q;
		    xi += q;
		}
		else if(theta == 1)
		{
		    q = -1;
		    sl = 4*q;
		    xi += 16*q;
		}
		else
		{
		    printf("Erro bizarro! theta: %d", theta);
		}
	
	
		if( random.randBool(0.5) )
		{
		    alpha = 1.5;
		    beta = 2;
		}
		else
		{
		    alpha = 2;
		    beta = 1.5;
		}
	
		A[ 3*l*n + l ] = alpha;
		A[ 3*l*n + l + m ] = beta;
		
		A[ (3*l+1)*n + l ] = 1.0;
		A[ (3*l+1)*n + l + m ] = -(beta + 1);
		
		A[ (3*l+2)*n + l ] = -(alpha + 1);
		A[ (3*l+2)*n + l + m ] = 1.0;
		
		c[3*l] = alpha + beta + alpha*beta;
		c[3*l + 1] = 0.0;
		c[3*l + 2] = 0.0;
		
		
		Q[ l*n + l ] += q;
		Q[ (l+m)*n + l + m ] += q;
		
		s[l] += sl;
		s[l+m] += sl;
    }
    
    
    for(l = nConcaveSubProbs; l < nConcaveSubProbs + nBilinearSubProbs; l++)
    {
		alpha = random.random();
		
		
		A[ 3*l*n + l ] = alpha;
		A[ 3*l*n + l + m ] = alpha + 1;
		
		A[ (3*l+1)*n + l ] = -(alpha + 1);
		A[ (3*l+1)*n + l + m ] = -alpha;
		
		A[ (3*l+2)*n + l ] = 1.0;
		A[ (3*l+2)*n + l + m ] = -1.0;
		
		c[3*l] = 3*alpha + 1;
		c[3*l + 1] = -(alpha + 1);
		c[3*l + 2] = 1.0;
		
		Q[l*n + l+m] = Q[ (l+m)*n + l] = 1.0; //it is not a diagonal position, but we have 0.5 multiplying the matrix.
		s[l] -= 1;
		s[l+m] -= 1;
		
		xi += 1;
    }
    
    
    for(l = nConcaveSubProbs + nBilinearSubProbs; l < m; l++)
    {
		alpha = random.random(5.0, 7.5);
		ro = random.randInt(0, 1);
		omega = random.randInt(0, 1);
		
		theta = 1 - ro*omega;
		
		A[ 3*l*n + l ] = -3.0;
		A[ 3*l*n + l + m ] = -2.0;
		
		A[ (3*l+1)*n + l ] = -2.0;
		A[ (3*l+1)*n + l + m ] = -3.0;
		
		A[ (3*l+2)*n + l ] = A[ (3*l+2)*n + l + m] = 1.0;
		
		c[3*l] = -alpha;
		c[3*l + 1] = -alpha;
		c[3*l + 2] = 3;
		
		
		Q[ l*n + l ] += 1.0;
		Q[ (l+m)*n + l + m ] += ro;
	
	
		if(theta == 0)
		{
		    s[l] = 1;
		    s[l+m] = ro;
		    xi += 0.5*(1 + ro);
		}
		else
		{
		    s[l] = 3;
		    s[l+m] = 3*ro;
		    xi += 0.5*(1 + ro)*9;
		}
    }
    
    
    // ********************************************************
    //Generating M
    
    
    M = (double *) malloc( n*n * sizeof(double) );
    d = (double *) malloc( n * sizeof(double) );
    v = (double *) calloc(n, sizeof(double));
    if(!v || !M)
    {
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
    
    numNzv = denv*n;
    
    for(l = 0; l < numNzv; l++)
    {
		do
		{
		    theta = random.randInt(0, n-1);
		}while( v[theta] != 0.0);
		
		v[theta] = random.random();
    }
    
    alpha = IQD_norm2(n, v);
    
    for(l = 0; l < n; l++)
		v[l] = (1.0/alpha) * v[l] ; //should be a unitary vector
    
    
    for(l = 0; l < n; l++)
		d[l] = pow(  10, -condd * (l/(n-1.0) )      );
	
	i = j = 0;
	theta = lowerBoundToGenerateD * ((int) pow(10, condd) );
	for(l = 0; l < n; l++)
	{
		d[l] = random.randInt(lowerBoundToGenerateD, theta);
		
		if(d[l] == lowerBoundToGenerateD)
			i++;
		else if(d[l] == theta)
			j++;
	}
	
	if(i == 0)
	{
		//we choose a coordinate to set at lowerBoundToGenerateD
		d[ random.randInt(0, n-1) ] = lowerBoundToGenerateD;
	}
	
	if(j == 0)
	{
		//we choose a coordinate to set at lowerBoundToGenerateD
		d[ random.randInt(0, n-1) ] = theta;
	}
	
	//now, we know the 2-norm condidition number is 10**condd
    
    
    //calculating 2*vv'
    for(i = 0; i < n; i++)
    {
		for(l = 0; l < n; l++)
		    M[i*n + l] = -2.0*v[i]*v[l];
    }
    
    //adding the eye matrix
    for(i = 0; i < n; i++)
		M[i*n + i] += 1.0; 
    
    
    //multiplying by D
    for(i = 0; i < n; i++)
    {
		for(l = 0; l < n; l++)
		{
		    M[i*n + l] = M[i*n + l]*d[i];
		}
    }
    
    
	printf("Q:\n");
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			printf("%f ", Q[i*n + j]);
		printf(";\n");
	}
	printf("\n");
	
	
	printf("A:\n");
	for(i = 0; i < cons; i++)
	{
		for(j = 0; j < n; j++)
			printf("%f ", A[i*n + j]);
		printf(";\n");
	}
	printf("\n");
	
	printf("c:\n");
	for(i = 0; i < n; i++)
		printf("%f ;", s[i]);
	printf("\n");
	
	printf("b:\n");
	for(i = 0; i < cons; i++)
		printf("%f ;", c[i]);
	printf("\n");
	
	printf("d:\n");
	for(i = 0; i < n; i++)
		printf("%f ;", d[i]);
	printf("\n");
	
	printf("v:\n");
	for(i = 0; i < n; i++)
		printf("%f ;", v[i]);
	printf("\n");
	
	printf("M:\n");
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			printf("%f ", M[i*n + j]);
		printf(";\n");
	}
	printf("\n");
    
	
	free(d); d = NULL;
    free(v); v = NULL;
	
	
    // ***********************************************************************************************************
    //Multiplying our matrices by M
    
    
    Ms = (double *) malloc( n*sizeof(double) );
    AM = (double *) malloc( cons*n*sizeof(double) );
    MQM= (double *) malloc( n*n*sizeof(double) );
    
    
    if(!Ms || !AM || !MQM )
    {
		code = IQD_MEMORY_ERROR;
		goto desallocate_memory;
    }
    
    
    
    //M'*s
    for(i = 0; i < n; i++)
    {
		alpha = 0.0;
		
		for(j = 0; j < n; j++)
		    alpha += M[j*n + i] * s[j];// M' * s
		    
		Ms[i] = alpha;
    }
    
    
    //A*M
    for(i = 0; i < cons; i++)
    {
		for(j = 0; j < n; j++)
		{
		    alpha = 0.0;
		    
		    for(l = 0; l < n; l++)
				alpha += A[i*n + l] * M[l*n + j];
		    
		    AM[i*n + j] = alpha;
		}
    }
    
    //M'*Q*M
    
    //first step: M'Q
    
    for(i = 0; i < n; i++)
    {
		for(j = 0; j < n; j++)
		{
		    alpha = 0.0;
		    for(l = 0; l < n; l++)
				alpha += M[l*n + i]* Q[l*n + j]; //M should be transposed
		    
		    MQM[i*n + j] = alpha;
		}
    }
    
    
    //now, we take Q to store M'Q * M
    
    for(i = 0; i < n; i++)
    {
		for(j = 0; j < n; j++)
		{
		    alpha = 0.0;
		    
		    for(l = 0; l < n; l++)
				alpha += MQM[i*n + l]*M[l*n + j];
		    
		    Q[i*n + j] = alpha;
		}
    }
    
    free(MQM);
    
    MQM = Q;
    Q = NULL;
    
    free(M);	M = NULL;
    
    free(A);	A = NULL;
    free(s);	s = NULL;
	
	if(factor != 1.0)
	{
		for(i = 0; i < n*n; i++)
			MQM[i] *= factor;
		
		for(i = 0; i < cons*n; i++)
			AM[i] *= factor;
		
		for(i = 0; i < n; i++)
			Ms[i] *= factor;
		
		for(i = 0; i < cons; i++)
			c[i] *= factor;
	}
	
	
	printf("MQM:\n");
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
			printf("%f ", MQM[i*n + j]);
		printf(";\n");
	}
	
	printf("\nAM:\n");
	for(i = 0; i < cons; i++)
	{
		for(j = 0; j < n; j++)
			printf("%f ", AM[i*n + j] );
		printf(";\n");
	}
	
	printf("\nMc:\n");
	for(i = 0; i < n; i++)
		printf("%f ;", Ms[i] );
	printf("\n");
	
	
	printf("\nb:\n");
	for(i = 0; i < cons; i++)
		printf("%f ;", c[i] );
	printf("\n\n");
    
	
	prob.setParametersAndAllocate(n, cons);
	
	prob.setObjQuadCoefsMatrix(MQM);
	
	prob.setObjLinearCoefficients(Ms);
	
	for(i = 0; i < cons; i++)
		prob.setConstraintLinearPart(i, &AM[i*n]);
	
	prob.setConstraintUpperBounds(c); //setQuadConstraintsRHS(c);
	
	
	//prob.printProblem();
    
    
desallocate_memory:
    
    if(A) 	free(A);
    if(Q)	free(Q);
    if(c)	free(c);
    if(s)	free(s);
    
    if(d)	free(d);
    if(v)	free(v);
    if(M)	free(M);
    
    if(Ms)	free(Ms);
    if(AM)	free(AM);
    if(MQM)	free(MQM);
    
    
    
    return code;
}
#else
{
	return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif



 
