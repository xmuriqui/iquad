/** That file contains a interface to AMPL modelig system. To work, is necessary
 * download and install the AMPL Solver Library (AMPL).
 * 
 * Unfourtunatelly, at the present moment, we only deal with the problem coming from AMPL:
 * 
 * Min c'x + 0.5 x'Qx
 * s. t.:
 *	A x + 0.5 x'QQ x {<=, =, >=} b
 *	g(x) <= 0
 *	
 *	l <= x <= u
 * 
 *	x_i is integer for i in I
 *
 * 
 * Note that we cannot consider convex nonlinear functions in the objetive function because we do not know how to handle separable functions using ASL and we have no information in the offical manual.
 * 
 * 
 * 
 * Author: Wendel Melo
 * 
 * Date: 18-jan-2013
 * 
 */


#include <cstdio>
#include <cstdlib>

#include <new>
#include <iostream>

#include "iquad.hpp"
#include "IQD_tools.hpp"
#include "MIP_ampl.hpp"


using namespace iquad;
using namespace minlpproblem;
using namespace std;


#define IQD_READ_AMPL_PARAMS 1


#if IQD_HAVE_ASL

	#include "asl.h"
	#include "getstub.h"

	#define asl cur_ASL
	
	class IQD_AMPLParamData
	{
	public:
		IQD_BranchAndBound *bb;
		
		IQD_GeneralSolverParams *solverParams;
		IQD_GeneralSolverParams *sdpParams;
	};
	
	
	
	char * IQD_readIquadAMPLParams(Option_Info *oi, keyword *kw, char *value, IQD_BranchAndBound *bb)
	{
		const char *type = kw->name;
		char subname[100], strValue[100];
		long int iValue;
		int aux, ret;
		double dValue;
		
		aux = sscanf(value, "%s %s", subname, strValue);
		
		if(aux == 2)
		{
			if(type[0] == 'i') //if( strcmp(type, "int") == 0 )
			{
				ret = sscanf(strValue, "%ld", &iValue);
				
				if(ret == 1)
					ret = bb->setIntegerParameter(subname, iValue);
				else
					ret = -1;
			}
			else if(type[0] == 'd') //if( strcmp(type, "dbl") == 0 )
			{
				ret = sscanf(strValue, "%lf", &dValue);
				
				if(ret == 1)
					ret = bb->setDoubleParameter(subname, dValue);
				else
					ret = -1;
			}
			else if(type[0] == 's')
			{
				ret = bb->setStringParameter(subname, strValue);
			}
			else
			{
				ret = -1;
				std::cerr << IQD_PREPRINT << "Invalid type of solver parameter specification: " << type << " It should be [int, dbl, str]" << std::endl;
			}
			
			
			if(ret == 0)
			{
				std::cout << IQD_PREPRINT << type << " iquad parameter " << subname << " to " << strValue << std::endl;
			}
			else //if(ret < 0) //if ret == 1, we should not show error messge
			{
				//fprintf(stderr, "Error at setting %s iquad parameter %s to %s.\n", type, subname, strValue );
				
				std::cerr << IQD_PREPRINT <<  (ret == IQD_NAME_ERROR ? "Name" : "Value") << "error at setting " << type << " iquad parameter " << subname << " to " << strValue << std::endl ;
			}
			
		}
		else
		{
			//fprintf(stderr, "Error at reading value of parameter: %s ", kw->name);
			std::cerr << IQD_PREPRINT << "Error at reading value of parameter: " << kw->name;
			if(aux >= 1)
			{
				//fprintf(stderr, "%s", subname);
				std::cerr << " " << subname;
			}
			
			std::cerr << std::endl;
		}
		
		
		
		if(aux >= 1)
		{
			//value should point to rest of string parameter
			value = value + strlen(subname); //pointer Arithmetic
			
			//avoiding white space between subname and pValue
			while( *value != '\0' && (*value == ' ' || *value == '\t') )
				value++; //pointer Arithmetic
		}
		
		
		if(aux >= 2)
		{
			value = value + strlen(strValue); //pointer Arithmetic
		}
		
		
		return value;
		
	}
	
	
	char * IQD_readSolverAMPLParams(Option_Info *oi, keyword *kw, char *value, IQD_GeneralSolverParams *sparams)
	{
		const char *type = kw->name;
		char subname[100], strValue[100];
		long int iValue;
		int aux, ret;
		double dValue;
		
		
		
		aux = sscanf(value, "%s %s", subname, strValue);
		
		
		if(aux == 2)
		{
			
			if(type[0] == 'i') //if( strcmp(type, "int") == 0 )
			{
				ret = sscanf(strValue, "%ld", &iValue);
				
				if( ret == 1 )
					ret = sparams->storeIntegerParameter(subname, iValue);
				else
					ret = -1;
			}
			else if(type[0] == 'd') //if( strcmp(type, "dbl") == 0 )
			{
				ret = sscanf(strValue, "%lf", &dValue);
				
				if( ret == 1 )
					ret = sparams->storeDoubleParameter(subname, dValue);
				else
					ret = -1;
			}
			else if(type[0] == 's') //if( strcmp(type, "str") == 0 ) //could be only else, but ok...
			{
				ret = sparams->storeStringParameter(subname, strValue);
			}
			else
			{
				ret = -1;
				std::cerr << IQD_PREPRINT <<  "Invalid type of solver parameter specification: " << type << ". It should be [int, dbl, str]." << std::endl;
			}
			
			if(ret == 0)
			{
				std::cout << IQD_PREPRINT << "Operation of setting parameter " << subname << " to " << strValue << " stored." << std::endl;
			}
			else
			{
				std::cerr << IQD_PREPRINT << "Error at setting " << type << " parameter " << subname << " to " << strValue << std::endl; 
			}
			
			
		}
		else
		{
			std::cerr << IQD_PREPRINT << "Error at reading value of parameter: " << kw->name;
			
			if(aux >= 1)
				std::cerr << " " << subname; //fprintf(stderr, "%s", subname);
			
			std::cerr << std::endl;
			//IQD_getchar();
		}
		
		
		
		if(aux >= 1)
		{
			//value should point to rest of string parameter
			value = value + strlen(subname); //pointer Arithmetic
			
			//avoiding white space between subname and pValue
			while( *value != '\0' && (*value == ' ' || *value == '\t') )
				value++; //pointer Arithmetic
		}
		
		
		if(aux >= 2)
		{
			value = value + strlen(strValue); //pointer Arithmetic
		}
		
		
		return value;
	}
	
	
	char * IQD_readAMPLParams(Option_Info *oi, keyword *kw, char *value)
	{
		char firstEnv = oi->opname[6];
		IQD_AMPLParamData *data = (IQD_AMPLParamData *) kw->info;
		
		
		//do not echo parameter name and value...
		oi->option_echo &= ~ASL_OI_echo;
		
		if( firstEnv == 'o' ) //( strcmp(&(oi->opname[6]), "options" )
		{
			return IQD_readIquadAMPLParams(oi, kw, value, data->bb);
		}
		else if( oi->opname[7] == 'u' ) //( strcmp(&(oi->opname[6]), "subsolver_options" )
		{
			return IQD_readSolverAMPLParams(oi, kw, value, data->solverParams);
		}
		else if( oi->opname[7] == 'd' ) //( strcmp(&(oi->opname[6]), "sdp_options" )
		{
			return IQD_readSolverAMPLParams(oi, kw, value, data->sdpParams);
		}
		
		
		IQD_PRINTERRORMSG( "Error! Unknow solver environment at IQD_readAMPLParam!");
		//IQD_getchar();
		return NULL;
		
	}
	
	
#endif






class IQD_ReadAmplModel : public MIP_ReadAmplModel
{
public:
	
	void readIquadParameters(IQD_AMPLParamData *paramData)
	{
	#if IQD_HAVE_ASL
		
		char **argv = NULL;
		
		static keyword keywords[3] = {
		{(char *)"dbl", IQD_readAMPLParams, paramData, (char *) "double float (real) parameters at iquad"}, 
		{(char *)"int", IQD_readAMPLParams, paramData, (char *) "integer parameters at iquad"}, 
		{(char *)"str", IQD_readAMPLParams, paramData, (char *) "string  parameters at iquad"}
		};
		
		
		static Option_Info opInfo = {(char *)"iquad", (char *) "iquad", (char *)"iquad_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 };
		
		
		static Option_Info opInfoSub = {(char *)"iquad", (char *) "iquad", (char *)"iquad_subsolver_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 };
		
		
		static Option_Info opInfoSdp = {(char *)"iquad", (char *) "iquad", (char *)"iquad_sdp_options", keywords, 3, 0, 0, 0, 0, 0, 0, 0 };
		
		
		argv = (char **) malloc( 2*sizeof(char*) );
		if(!argv)
			goto desallocate_memory;
		
		argv[0] = (char *) malloc( sizeof(char) );
		if(!argv)
			goto desallocate_memory;
			
		argv[1] = NULL;
		
		argv[0][0] = '\0';
		
		
		readParameters(argv, &opInfo);
		readParameters(argv, &opInfoSub);
		readParameters(argv, &opInfoSdp);
		
		
		
		
	desallocate_memory:
		
		if(argv)
		{
			if( argv[0] )	free(argv[0]);
			free(argv);
		}
		
	#endif
	}
	
	
};









int iquad::IQD_ampl(char *stub, char **argvo, const char *outDirectory)
#if IQD_HAVE_ASL
{
	//int i;
	int code, ret;
	
	char myMsg[200];
	double aux;
	double *psol, *dsol;
    
    MIP_NonLinearEval *myEval = NULL;//MIP_NLEvalAmlp *myEval = NULL;
	IQD_ReadAmplModel reader;
	
    IQD_IQuadProb prob;
    //IQD_RefProb refProb;
    IQD_GeneralSolverParams solverParams, sdpParams;
    IQD_BranchAndBound bb;
	//IQD_DummyNLEvalAmpl dummyEval;
	
	IQD_AMPLParamData paramData;
	
	paramData.bb = &bb;
	paramData.solverParams = &solverParams;
	paramData.sdpParams = &sdpParams;
	
	
	
    
	IQD_welcome();
	
	cout << IQD_PREPRINT << "Reading user model:" << endl;
	ret = reader.readProblem(stub, prob, myEval);
	
	if(ret != 0)
	{
		cerr << IQD_PREPRINT << "Error " << ret << " at reading ampl model!" << endl;
		
		goto desallocate_memory;
		code = IQD_UNDEFINED_ERROR;
	}
	
	cout << IQD_PREPRINT << "Reading user parameters:" << endl;
	reader.readIquadParameters(&paramData);
	
	
	if( reader.isMaximizationProblem() )
	{
		//we have to reverse possible lower and upper bounds given by the user
		
		aux = bb.in_upper_bound;
		
		if( bb.in_lower_bound > -IQD_INFINITY )
			bb.in_upper_bound = -bb.in_lower_bound;
		
		if( aux < IQD_INFINITY )
			bb.in_lower_bound = -aux;
	}
	
	
	
	//prob.print();
	
	
	if( reader.isMaximizationProblem() )
		cout << IQD_PREPRINT << "Maximization problem addressed as minimization" << endl;
	
	
	
	/*{
		IQD_RefProb2 refProb;
		
		
		refProb.oprob = &prob;
		refProb.zero_tol = bb.in_zero_tol_to_split;
		refProb.splitWay = IQD_SQ_MIX_SCHUR_DIAG_SDP;
		
		cout << "refProb.zero_tol: " << refProb.zero_tol << endl;
		
		int r = refProb.buildRefProb( bb.in_sdp_solver, bb.in_qp_solver, bb.in_qcp_solver, bb.in_nlp_solver, NULL, NULL );
		
		cout << "ret: " << r << endl;
		
		refProb.rprob.writeMIQCPModelInAMPLFile( "refprob2.mod" );
		
		cout << "Escrevi arquivo refprob2.mod" << endl;
		
		cout << endl << endl;
		cout << "--------------------------------------------------------------------------------------------------" << endl;
		
		//printing v as colunm order
		cout << "dimv: " << refProb.dimv << endl;
		refProb.v.printSparseMatrix();
		
		cout << "v: " << endl;
		refProb.v.printAllMatrix();
		//for(int i = 0; i < prob.n; i++)
		{
			//for(int j = 0; j < refProb.dimv; j++)
				//cout << refProb.v[j][i] << " ";
			
			//cout << ";" << endl;
		//} 
		
		cout << "lambdaObj: " << endl;
		for(int i = 0; i < refProb.dimv; i++)
			cout <<  refProb.lambdaObj[i] << " ";
		cout << endl;
		
		
		cout << "Q: " << endl;
		prob.Q.printAllMatrix();
		
		cout << "P: " << endl;
		refProb.rprob.Q.printAllMatrix();
		
		
		for(int i = 0; i < refProb.rprob.m; i++ )
		{
			if( refProb.rprob.QC[i].getNumberOfElements() > 0 )
			{
				cout << endl << "lambdaConstr " << i << ": ";
				for(int j = 0; j < refProb.dimv; j++)
					cout << refProb.lambdaConstr[i][j] << " " ;
				cout << endl;
				
				cout << "QC " << i << ": " ;
				if( i >= prob.m )
					cout << "(-1.0)*" ;
				
				cout << endl;
				prob.QC[i < prob.m ? i : refProb.origConstr[i]].printAllMatrix();
				
				cout << "P: " << endl;
				refProb.rprob.QC[i].printAllMatrix();
			}
		}
		
		
		
		//cout << "v: " << refProb.v << endl;
		//IQD_printMatrix( refProb.dimv, prob.n, &refProb.v[0][0] );
		
		
		//getchar();
	} */
	
	
	#if IQD_FILTER_PROBLEMS
	{
		const int nmax = 200;
		bool convexObj, convexConstr;
		char cmd[500];
		
		const double zero_tol = 1.0e-6;
		
		
		if( prob.getNumberOfIntegerVars() > 0 )
		{
			sprintf(cmd, "mkdir integer; mv %s integer", stub);
			system(cmd);
		}
		
		goto desallocate_memory;
		
		
		if( prob.n > nmax )
		{
			sprintf( cmd, "mv %s ./huge", stub );
			system("mkdir huge");
			system(cmd);
			
			goto desallocate_memory;
		}
		
		
		IQD_isConvexProblemAboutQuadratics( prob, zero_tol, convexObj, convexConstr );
		
		if( convexConstr )
		{
			sprintf( cmd, "mv %s ./convex", stub );
			system("mkdir convex");
			system(cmd);
			
			goto desallocate_memory;
		}
	}
	#endif
	
	
	IQD_runProblem(stub, outDirectory, prob, bb, &solverParams, &sdpParams);
	
	
	#if IQD_FILTER_PROBLEMS
	//we filter here problems having no bounded feasilbe region and problems too huge...
	{
		char cmd[500];
		if( bb.out_return_code == IQD_SPLITING_ERROR && bb.out_return_subcode == IQD_UNBOUNDED_PROBLEM)
		{
			sprintf( cmd, "mv %s ./waste", stub );
			system("mkdir waste");
			system(cmd);
		}
		else if( bb.out_return_code == IQD_OPTIMAL_SOLUTION && bb.out_number_of_iterations == 1 )
		{
			sprintf( cmd, "mv %s ./convex", stub );
			system("mkdir convex");
			system(cmd);
		}
		else
		{
			sprintf( cmd, "mv %s ./good", stub );
			system("mkdir good");
			system(cmd);
		}
	}
	#endif
	
	
	switch( bb.out_return_code )
	{
		case IQD_OPTIMAL_SOLUTION:
		    ret = 0;
			break;
		
		case IQD_UNBOUNDED_PROBLEM:
			ret = 300;
			break;
			
		case IQD_MAX_ITERATIONS_STOP:
			ret = 400;
			break;
		
		case IQD_MAX_TIME_STOP:
			ret = 401;
			break;
		
		case IQD_INFEASIBLE_PROBLEM:
			ret = 200;
			break;
		
		case IQD_LIBRARY_NOT_AVAILABLE:
			ret = 501;
			break;
			
		case IQD_CALLBACK_FUNCTION_ERROR:
			ret = 502;
			break;
		
		case IQD_SDP_SOLVING_ERROR:
			ret = 503;
			break;
			
		case IQD_QCP_SOLVER_ERROR:
			ret = 504;
			break;
		
		case IQD_NLP_SOLVER_ERROR:
			ret = 505;
			break;
			
		default:
			ret = 500;
	}
	
	
	if( bb.out_best_obj < IQD_INFINITY  )
    {
		if( reader.isMaximizationProblem() )
			bb.out_best_obj = -bb.out_best_obj;
		
		
		if(bb.out_return_code == IQD_OPTIMAL_SOLUTION)
		    sprintf( myMsg, "An optimal solution was found! Obj function: %f. CPU Time: %f", bb.out_best_obj, bb.out_cpu_time );
		else
		    sprintf( myMsg, "A feasible solution was found. Obj function: %f. CPU Time: %f", bb.out_best_obj, bb.out_cpu_time );
		
		psol= bb.out_best_sol;
		dsol = NULL;
    }
    else
    {
		sprintf( myMsg, "No feasible solution was found. CPU Time: %f", bb.out_cpu_time );
		
		psol = dsol = NULL;
    }
	
	
	reader.putInfoToAmpl(myMsg, ret, psol, dsol);
	
	
	
	
	
	code = 0;
	
desallocate_memory:
	
	if(myEval) delete myEval;
	
	return code;
}


#else
{
    return IQD_LIBRARY_NOT_AVAILABLE;
}
#endif






