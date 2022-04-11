/*
 * File that implement a function to run iquad branch and and bound on a MINLPProblem
 */

#include <new>

#include "iquad.hpp"
#include "IQD_tools.hpp"

#if IQD_USE_MURIQUI_ON_INT_ORIGINAL_PROBLEM
    #include "muriqui.hpp"
#endif





using namespace std;
using namespace iquad;




#if IQD_SAVE_OUTPUT_FILE_BY_HOUR
	extern iquad::IQD_OutFile *IQD_outFileHour;
#endif




int iquad::IQD_runProblem( const char *fileName, const char *outDirectory, IQD_IQuadProb &prob, IQD_BranchAndBound &bb, IQD_GeneralSolverParams *solverParams, IQD_GeneralSolverParams *sdpSolverParams )
{
	int i;
	int code;
	IQD_OutFile *outFile = NULL;
	
	
	#if IQD_SAVE_OUTPUT_FILE
		char *nameOut = new (nothrow) char[ strlen(outDirectory) + strlen(IQD_OUT_FILE_NAME) + 2 ];
		
		strcpy( nameOut, outDirectory );
		nameOut[ strlen(outDirectory) ] = '/';
		nameOut[ strlen(outDirectory) + 1] = '\0'; 
		strcat( nameOut, IQD_OUT_FILE_NAME );
	#endif
	
	#if IQD_SAVE_OUTPUT_FILE_BY_HOUR
		char *nameOutHour = new (nothrow) char[ strlen(outDirectory) + strlen(IQD_OUT_FILE_BY_HOUR_NAME) + 2 ];
		
		strcpy( nameOutHour, outDirectory );
		nameOutHour[ strlen(outDirectory) ] = '/';
		nameOutHour[ strlen(outDirectory) + 1] = '\0';
		strcat( nameOutHour, IQD_OUT_FILE_BY_HOUR_NAME );
	#endif
	
	
	#if IQD_APPLY_JON_TRANSFORMATION
	{
		char line[500];
		
// 		
		IQD_scaleProblemAndTurnIntegerVars(prob, 10.0);
		
		/*system("mkdir intboxqpmod");
		
		sprintf(line, "intboxqpmod/%s", fileName );
		
		prob.writeMIQCPModelInAMPLFile( line);
		
		goto termination; */
	}
	#endif
	
	
	
	#if 0
	{
		#define IQD_AMPL_PATH	"/home/wendel/ampl/ampl"
		char name[500], cmd[500];
		
		
		sprintf( name, "%s.mod", fileName );
		sprintf( cmd, "iaia; write g%s", name); //using option solver to generate nl file
		
		prob.writeMIQCPModelInAMPLFile(name, cmd );
		
		sprintf(cmd, IQD_AMPL_PATH " %s", name);
		
		system(cmd);
		
		goto termination;
	}
	#endif
	
	
	#if 0
	{
		char name[500], cmd[500];
		
		
		sprintf( name, "%s.mod", fileName );
		prob.writeMIQCPModelInAMPLFile(name);
		
		goto termination;
	}
	#endif
	
	
	#if IQD_SAVE_OUTPUT_FILE
		outFile = new (nothrow) IQD_OutFile;
		
		if(outFile)
		{
			printf("nameOut: %s\n", nameOut);
			
			i = outFile->open( nameOut );
			if(i != 0)
			{
				delete outFile;
				outFile = NULL;
			}
			else
			{
				outFile->writeBreakLine();
				outFile->writeProblemInfo(fileName, prob);
				outFile->flush();
			}
			
			//goto termination;
		}
	#endif
	
	#if IQD_SAVE_OUTPUT_FILE_BY_HOUR
	
		IQD_outFileHour = new (nothrow) IQD_OutFile;
		
		if(IQD_outFileHour)
		{
			i = IQD_outFileHour->open( nameOutHour );
		    if(i != 0)
		    {
				delete IQD_outFileHour;
				IQD_outFileHour = NULL;
		    }
		    else
		    {
				IQD_outFileHour->writeBreakLine();
				IQD_outFileHour->writeProblemInfo(fileName, prob);
				IQD_outFileHour->flush();
		    }
		}
	
	#endif
	
	if(prob.getNumberOfVars() <= 50)
		prob.print();
	
	
	
	#if IQD_USE_CPLEX_ON_NONCONVEX_PROBLEM
	{
		const double MAXCPUTIME = 7200;
		int r, niter;
		double lb, cputime, realTime;
		clock_t clockStart;
		optsolvers::OPT_Cplex cplex;
		
		
		
		if( prob.hasQuadMatrixInSomeConstraint() )
		{
			r = optsolvers::OPT_setQCPProblemOnCplex(prob, &cplex, true, true, true, true);
		}
		else
		{
			r = optsolvers::OPT_setMINLPProblem(prob, &cplex, true, true, true, true);
		}
		
		#if IQD_DEBUG_MODE
			if( r != 0 )
			{
				IQD_PRINTERRORNUMBER(r);
				IQD_getchar();
			}
		#endif
		
		
		
		r = cplex.setNumberOfThreads(1);
		
		#if IQD_DEBUG_MODE
			if( r != 0 )
			{
				IQD_PRINTERRORNUMBER(r);
				IQD_getchar();
			}
		#endif
		
		r = cplex.setIntegerParameter( "CPX_PARAM_SOLUTIONTARGET", 3 );
		
		if( r != 0 )
		{
			IQD_PRINTERRORNUMBER(r);
			std::cout << IQD_PREPRINT << "Error to set parameter about nonconvex optimization in cplex" << std::endl;
			IQD_getchar();
		}
		
		r = cplex.setMaxCPUTime(MAXCPUTIME);
		
		#if IQD_DEBUG_MODE
			if( r != 0 )
			{
				IQD_PRINTERRORNUMBER(r);
				IQD_getchar();
			}
		#endif
		
		
		realTime = IQD_getTime();
		clockStart = clock();
		
		cplex.generateModelFile("cplex.lp");
		IQD_getchar();
		
		cplex.setOutputLevel(10);
		
		
		r = cplex.solve();
		
		
		
		cputime = (double( clock() - clockStart ))/CLOCKS_PER_SEC;
		realTime = IQD_getTime() - realTime;
		
		
		//cout << "CPXgetitcnt: " << CPXgetitcnt( cplex.env, cplex.prob ) << std::endl;
		//cout << "CPXgetbaritcnt: " << CPXgetbaritcnt( cplex.env, cplex.prob ) << std::endl;
		//cout << "CPXgetmipitcnt: " << CPXgetmipitcnt( cplex.env, cplex.prob ) << std::endl;
		//cout << "CPXgetnodecnt: " << CPXgetnodecnt( cplex.env, cplex.prob ) << std::endl;
		
		if( prob.getNumberOfIntegerVars() > 0 )
		{
			niter = CPXgetnodecnt( cplex.env, cplex.prob );
			/*r = CPXgetmiprelgap( cplex.env, cplex.prob, &gap );
			if( r != 0 )
			{
				IQD_PRINTERRORNUMBER(r);
				IQD_getchar();
			} */
			
			/*r = CPXgetbestobjval( cplex.env, cplex.prob, &lb );
			if( r != 0 )
			{
				IQD_PRINTERRORNUMBER(r);
				IQD_getchar();
			} */
		}
		else
		{
			niter = CPXgetitcnt( cplex.env, cplex.prob );
		}
		
		
		
		lb = cplex.dualObjValue;
		
		
		std::cout << "____________________________________________________________________________" << std::endl;
		
		std::cout << "cplex Problem: " << fileName << " return code: " << cplex.retCode << " obj function: " << cplex.objValue << " lower bound: " << lb << " time: " << realTime << " cpu time: " << cputime << " iters: " << niter <<  std::endl;
		
		std::cout << "____________________________________________________________________________" << std::endl;
		
		
		#if IQD_SAVE_OUTPUT_FILE
			if(outFile)
			{
				outFile->writeString("cplex");
				outFile->writeInt(cplex.retCode);
				outFile->writeSolInfo(cplex.origSolverRetCode, lb, cplex.objValue, cputime, realTime, niter, NAN, NAN, NAN);
				outFile->flush();
			}
		#endif
		
		
		
		
		goto termination;
		
	}
	#endif
	
	
	
	
	
	
	
	/* sParams.addDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_REL_GAP", 1.0e-4); //default: 1.0e-6
	sParams.addDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_MU_RED", 1.0e-6); //default: 1.0e-12
	sParams.addDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_PFEAS", 1.0e-6); //default: 1.0e-8
	sParams.addDoubleParameter( "MSK_DPAR_INTPNT_NL_TOL_DFEAS", 1.0e-6); //default: 1.0e-8
	sParams.addDoubleParameter( "MSK_DPAR_INTPNT_TOL_INFEAS", 1.0e-6);*/
	
	solverParams->storeIntegerParameter( "MSK_IPAR_NUM_THREADS", 1 );
	
	solverParams->storeDoubleParameter( "MSK_DPAR_OPTIMIZER_MAX_TIME", 2*60*60);
	
	solverParams->storeIntegerParameter("MSK_IPAR_INTPNT_MAX_ITERATIONS", 2000);
	//solverParams->addIntegerParameter("max_iter", 5000);
	//solverParams->addIntegerParameter("MaxIter", 1500);
	

    sdpSolverParams->storeIntegerParameter("MSK_IPAR_NUM_THREADS", 1 );
    sdpSolverParams->storeDoubleParameter("MSK_DPAR_OPTIMIZER_MAX_TIME", 2*60*60);
	
	//sParams.addIntegerParameter("print_level", 4);
	//sParams.addStringParameter("derivative_test", "second-order");
    
    //mosekParams.addDoubleParameter("MSK_DPAR_INTPNT_CO_TOL_REL_GAP", 1.0e-10);
    //mosekParams.addDoubleParameter("", );
    //mosekParams.addDoubleParameter("", );
    
	
	
	
	#if IQD_TEST_CONVEX_PROBLEM
	{
		bool convex, cobj, cconstr;
		int r;
		const double zeroTol = 1.0e-8;
		
		r = IQD_isConvexProblemAboutQuadratics( prob, zeroTol, cobj, cconstr );
		
		if( r != 0 )
		{
			IQD_PRINTERRORNUMBER( r );
		}
		else
		{
			convex = cobj && cconstr;
			
			if( convex )
			{
				char cmd[ IQD_max(strlen(fileName), strlen(IQD_CONVEX_DIR_MOVING)) + 50];
				
				std::cout << "Convex problem" << std::endl;
				
				sprintf(cmd, "mkdir %s", IQD_CONVEX_DIR_MOVING);
				system( cmd ); //if folder already exist, there is no difference
				
				sprintf( cmd, "mv %s %s", fileName, IQD_CONVEX_DIR_MOVING );
				system( cmd );
			}
			else
			{
				std::cout << "Nonconvex problem" << endl;
			}
		}
		
		goto termination;
	}
	#endif
	
	
	
	#if IQD_USE_MURIQUI_ON_INT_ORIGINAL_PROBLEM
		if( prob.getNumberOfIntegerVars() > 0 )
		{
			int r;
			double obj;
			const double MAXHEURTIME = 60;
			muriqui::MRQ_HeuristicExecutor heurExec;
			muriqui::MRQ_GeneralSolverParams nlpParams;
			
			heurExec.in_use_diving = false;
			heurExec.in_use_igma2 = false;
			heurExec.in_use_oafp = false;
			
			heurExec.setSolvers(muriqui::MRQ_CPLEX, muriqui::MRQ_IPOPT);
			
			heurExec.setMaxTimes(MAXHEURTIME);
			heurExec.in_max_total_cpu_time = MAXHEURTIME;
			
			nlpParams.addDoubleParameter( "max_cpu_time", MAXHEURTIME );
			
			std::cout << IQD_PREPRINT << "Aplying Muriqui heuristcs..." << std::endl;
			
			r = heurExec.run(prob, NULL, &nlpParams, MRQ_INFINITY, obj, NULL, NULL );
			
			if( r == muriqui::MRQ_HEURISTIC_SUCCESS || r == muriqui::MRQ_OPTIMAL_SOLUTION )
			{
				std::cout << IQD_PREPRINT << "Muriqui found a feasible solution: " << obj << std::endl;
				bb.in_upper_bound = obj; // + 0.01 + 0.01*obj;
				
				
				for(int i = 0; i < prob.n; i++)
					std::cout << "x[" << i << "]: " << prob.x[i] << endl ;
				
				
				double myobj;
				prob.objEval(1, true, prob.x, myobj);
				
				std::cout << "myobj: " << myobj << std::endl;
			}
			
			std::cout << "iters: " << heurExec.fp.out_number_of_iterations << std::endl;
		}
	#endif
	
	
	
	#if IQD_FILTER_PROBLEMS
	{
		
		bb.in_max_cpu_time = 60*60;
		bb.in_max_iters = 10;
		bb.in_split_way = IQD_SQ_MIX_SCHUR_DIAG_SDP;
		bb.in_set_special_solver_parameters = false;
		bb.in_abs_feas_tol = 1.0e-5;
		bb.in_printing_frequency = 1;
		bb.in_qcp_solver = IQD_QCP_IPOPT;
		
		//bb.in_sdp_solver = IQD_SS_CSDP;
		
		bb.in_number_of_threads = 1;
		
		bb.run(prob, solverParams, sdpSolverParams);
		printf("____________________________________________________________________________\n");
		printf("Problem: %s Return code: %d obj function: %f time: %0.2f cpu time: %0.2f splitting time: %0.2f iters: %d\n", fileName, bb.out_return_code, bb.out_best_obj, bb.out_clock_time, bb.out_cpu_time, bb.out_splitting_cpu_time, (int) bb.out_number_of_iterations);
		//for(i = 0; i < n; i++)
		//printf("x[%d]: %f\n", i, bb.out_best_sol[i]);
		printf("____________________________________________________________________________\n");
		
		
		goto termination;
	}
	#endif
	
	
	
	
	
    bb.in_printing_frequency = 1;//10000;
    bb.in_max_cpu_time = 2*60*60;
	bb.in_print_level = 100;
	
	//bb.in_relative_convergence_tol = 1.0e-6;
	//bb.in_recalculate_aux_bounds_on_zu_updt = true;
	//bb.in_set_special_solver_parameters = false;
	
	//bb.in_relative_eps_added_to_diag_P = 1.0e-6;
	//bb.in_absolute_eps_added_to_diag_P = 1.0e-5;
	
    bb.in_number_of_threads = 1;
	bb.in_integrality_strategy = IQD_IS_ON_BB_ENUMERATION;
	bb.in_set_special_solver_parameters = false;
	
	bb.in_use_heuristics_on_integer_probs = false;
	
	//bb.in_qp_solver = IQD_QP_CPLEX;
	
    bb.in_split_way = IQD_SQ_SCHUR;
    //bb.in_split_way = IQD_SQ_DIAG_SDP;
    //bb.in_split_way = IQD_SQ_DIAG_DOM;
    //bb.in_split_way = IQD_SQ_IDENTITY;
	//bb.in_split_way = IQD_SQ_MILP_BASIS;
	bb.in_split_way = IQD_SQ_MIX_SCHUR_DIAG_SDP;
	
	//prob.setNLEvalObject(&eval, true);
	
    bb.run(prob, solverParams, sdpSolverParams);
    
    printf("____________________________________________________________________________\n");
    printf("Problem: %s Return code: %d obj function: %f time: %0.2f cpu time: %0.2f splitting time: %0.2f iters: %d\n", fileName, bb.out_return_code, bb.out_best_obj, bb.out_clock_time, bb.out_cpu_time, bb.out_splitting_cpu_time, (int) bb.out_number_of_iterations);
    //for(i = 0; i < n; i++)
	//printf("x[%d]: %f\n", i, bb.out_best_sol[i]);
    printf("____________________________________________________________________________\n");
	
	
	if( bb.out_return_code == IQD_INFEASIBLE_PROBLEM && bb.in_upper_bound < IQD_INFINITY )
		bb.out_return_code = IQD_OPTIMAL_SOLUTION;
    
    
    #if IQD_SAVE_OUTPUT_FILE
		if(outFile)
		{
		    outFile->writeSolInfo(bb.out_return_code, bb.out_lower_bound, bb.out_upper_bound, bb.out_cpu_time, bb.out_clock_time, bb.out_number_of_iterations, bb.out_splitting_cpu_time, bb.out_root_node_lower_bound, bb.out_root_node_upper_bound);
		    outFile->flush();
		}
    #endif
	
	
	/*bb.in_split_way = IQD_SQ_DIAG_SDP;
	bb.run(prob, solverParams, sdpSolverParams);
    printf("____________________________________________________________________________\n");
    printf("Problem: %s Return code: %d obj function: %f time: %0.2f cpu time: %0.2f splitting time: %0.2f iters: %d\n", fileName, bb.out_return_code, bb.out_best_obj, bb.out_clock_time, bb.out_cpu_time, bb.out_splitting_cpu_time, (int) bb.out_number_of_iterations);
    //for(i = 0; i < n; i++)
	//printf("x[%d]: %f\n", i, bb.out_best_sol[i]);
    printf("____________________________________________________________________________\n");
	
	if( bb.out_return_code == IQD_INFEASIBLE_PROBLEM && bb.in_upper_bound < IQD_INFINITY )
		bb.out_return_code = IQD_OPTIMAL_SOLUTION;
    
    
    #if IQD_SAVE_OUTPUT_FILE
		if(outFile)
		{
		    outFile->writeSolInfo(bb.out_return_code, bb.out_lower_bound, bb.out_upper_bound, bb.out_cpu_time, bb.out_clock_time, bb.out_number_of_iterations, bb.out_splitting_cpu_time, bb.out_root_node_lower_bound, bb.out_root_node_upper_bound);
		    outFile->flush();
		}
    #endif */
	
	
	
	
	/*bb.in_split_way = IQD_SQ_DIAG_DOM;
	bb.run(prob, solverParams, sdpSolverParams);
    printf("____________________________________________________________________________\n");
    printf("Problem: %s Return code: %d obj function: %f time: %0.2f cpu time: %0.2f splitting time: %0.2f iters: %d\n", fileName, bb.out_return_code, bb.out_best_obj, bb.out_clock_time, bb.out_cpu_time, bb.out_splitting_cpu_time, (int) bb.out_number_of_iterations);
    //for(i = 0; i < n; i++)
	//printf("x[%d]: %f\n", i, bb.out_best_sol[i]);
    printf("____________________________________________________________________________\n");
	
	
	if( bb.out_return_code == IQD_INFEASIBLE_PROBLEM && bb.in_upper_bound < IQD_INFINITY )
		bb.out_return_code = IQD_OPTIMAL_SOLUTION;
    
    
    #if IQD_SAVE_OUTPUT_FILE
		if(outFile)
		{
		    outFile->writeSolInfo(bb.out_return_code, bb.out_lower_bound, bb.out_upper_bound, bb.out_cpu_time, bb.out_clock_time, bb.out_number_of_iterations, bb.out_splitting_cpu_time, bb.out_root_node_lower_bound, bb.out_root_node_upper_bound);
		    outFile->flush();
		}
    #endif */
	
	
	
	
	
	
	/*bb.in_split_way = IQD_SQ_IDENTITY;
	bb.run(prob, solverParams, sdpSolverParams);
    printf("____________________________________________________________________________\n");
    printf("Problem: %s Return code: %d obj function: %f time: %0.2f cpu time: %0.2f splitting time: %0.2f iters: %d\n", fileName, bb.out_return_code, bb.out_best_obj, bb.out_clock_time, bb.out_cpu_time, bb.out_splitting_cpu_time, (int) bb.out_number_of_iterations);
    //for(i = 0; i < n; i++)
	//printf("x[%d]: %f\n", i, bb.out_best_sol[i]);
    printf("____________________________________________________________________________\n");
	
	
	if( bb.out_return_code == IQD_INFEASIBLE_PROBLEM && bb.in_upper_bound < IQD_INFINITY )
		bb.out_return_code = IQD_OPTIMAL_SOLUTION;
    
    
    #if IQD_SAVE_OUTPUT_FILE
		if(outFile)
		{
		    outFile->writeSolInfo(bb.out_return_code, bb.out_lower_bound, bb.out_upper_bound, bb.out_cpu_time, bb.out_clock_time, bb.out_number_of_iterations, bb.out_splitting_cpu_time, bb.out_root_node_lower_bound, bb.out_root_node_upper_bound);
		    outFile->flush();
		}
    #endif */
	
	
	
	/*bb.in_split_way = IQD_SQ_MIX_SCHUR_DIAG_SDP;
	bb.run(prob, solverParams, sdpSolverParams);
    printf("____________________________________________________________________________\n");
    printf("Problem: %s Return code: %d obj function: %f time: %0.2f cpu time: %0.2f splitting time: %0.2f iters: %d\n", fileName, bb.out_return_code, bb.out_best_obj, bb.out_clock_time, bb.out_cpu_time, bb.out_splitting_cpu_time, (int) bb.out_number_of_iterations);
    //for(i = 0; i < n; i++)
	//printf("x[%d]: %f\n", i, bb.out_best_sol[i]);
    printf("____________________________________________________________________________\n");
	
	
	if( bb.out_return_code == IQD_INFEASIBLE_PROBLEM && bb.in_upper_bound < IQD_INFINITY )
		bb.out_return_code = IQD_OPTIMAL_SOLUTION;
    
    
    #if IQD_SAVE_OUTPUT_FILE
		if(outFile)
		{
		    outFile->writeSolInfo(bb.out_return_code, bb.out_lower_bound, bb.out_upper_bound, bb.out_cpu_time, bb.out_clock_time, bb.out_number_of_iterations, bb.out_splitting_cpu_time, bb.out_root_node_lower_bound, bb.out_root_node_upper_bound);
		    outFile->flush();
		}
    #endif */
	
	
	
	
	
	
	
	code = bb.out_return_code;
	
	
termination:
	
	#if IQD_SAVE_OUTPUT_FILE
		if(outFile)
		    outFile->close();
	#endif
	
	
	#if IQD_SAVE_OUTPUT_FILE_BY_HOUR
		if(IQD_outFileHour)
		{
			IQD_outFileHour->close();
			delete IQD_outFileHour;
		}
	#endif
	
	#if IQD_SAVE_OUTPUT_FILE
		delete[] nameOut;
	#endif
	
	#if IQD_SAVE_OUTPUT_FILE_BY_HOUR
		delete[] nameOutHour; 
	#endif
	
	
	if(outFile)	delete outFile;
	
	return code;
}

